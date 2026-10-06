# SpectralMatching_v2.R  --  rewritten search engine for the DynamicTol package with:
#   S1  Library meta + peaks are read ONCE, not once per query.
#   S2  Candidate lookup by binary search on sorted precursor m/z.
#   S3  Peak matching by index intervals (findInterval + sequence) instead of
#       outer() + reshape2::melt + two full_join()s. No dense matrix, no join keys
#       built from rounded floats, no reshape2 (retired) dependency.
#   S4  Decoy spectra are built ONCE PER QUERY, not once per library candidate.
#   S5  Every shifted copy is tagged with tau, so the mean-over-shifts decoy score
#       required by Eq. (6)-(7) comes out of the same pass as the pooled one.
#   S6  F_CalPMMT() is the single source of truth for the tolerance
#   S7  Scores are returned numeric, not character.
#   S8  on.exit(add = TRUE); connections actually close.
#   S9  rttol reads dbDetails$rttol, not a global.
#   TolFactor                 = 1      matching window is +/- TolFactor * PMMT
#   Round               = TRUE   round m/z to 5 dp before comparing
#   MinMatch                  = 2      score forced to 0 below this many matched peaks

suppressPackageStartupMessages({
  library(data.table)
  library(foreach)
})
#' Global, m/z-sorted peak index of a store: built once, used when usePrecursors = FALSE.
F_BuildPeakIndex <- function(store, Round = TRUE) {
  len   <- store$meta$end - store$meta$start + 1L
  owner <- integer(length(store$mz))                 # peak -> row of store$meta
  owner[sequence(len, from = store$meta$start)] <- rep.int(seq_len(store$n), len)
  mz <- if (Round) round(store$mz, 5) else store$mz
  o  <- order(mz)
  list(mz = mz[o], row = owner[o])
}

#  Tolerance model -- the ONLY place the PMMT equations are written down
#' Peak matching mass tolerance (Da), vectorised over mz.
#'
#' Returns PMMT as defined in Eq. (1a) (Orbitrap) / Eq. (2a) (Q-TOF) of the
#' manuscript. NOTE: PMMT is a standard deviation. The window actually applied
#' during matching is +/- TolFactor * PMMT (see F_MatchPeaks).
#'
#' @param mz          numeric vector of m/z values
#' @param instrument  'Orbitrap' or 'Qtof'
#' @param MRP         reference resolving power at RefMZ (MS/MS)
#' @param RefMZ       reference m/z (Orbitrap only; 200 in the manuscript)
#' @param MF          mass fluctuation / mass accuracy in ppm.
#'                    Defaults follow the manuscript: 1 ppm Orbitrap, 2 ppm Q-TOF.
#' @param SDRatio     multiplier on the resolution term only (not on MF)
#' @param mztol       NA -> dynamic; > 1 -> fixed ppm; <= 1 -> fixed Da
F_CalPMMT <- function(mz,
                      instrument = c("Orbitrap", "Qtof"),
                      MRP = 17500, RefMZ = 200, MF = 1,
                      SDRatio = 1, mztol = NA) {
  
  instrument <- match.arg(instrument)
  if (is.null(MF)) MF <- if (instrument == "Orbitrap") 1 else 2
  if (is.character(mztol)) mztol <- suppressWarnings(as.numeric(mztol))
  if (length(mztol) != 1L) stop("mztol must be a single value")
  if (is.na(mztol)) {
    if (instrument == "Orbitrap") {
      # sigma = FWHM / 2.35482, FWHM = mz^1.5 / (MRP * sqrt(RefMZ))
      B <- SDRatio / (2.35482 * MRP * sqrt(RefMZ))
      return(B * mz^1.5 + MF * mz * 1e-6)
    } else {
      # For a Q-TOF this is a CONSTANT in ppm: 1e6/(2.35482*MRP) + MF
      return(SDRatio * mz / (2.35482 * MRP) + MF * mz * 1e-6)
    }
  } else if (mztol > 1) {
    return(mz * mztol * 1e-6)          # fixed ppm
  } else {
    return(rep_len(mztol, length(mz)))  # fixed Da
  }
}
#' Reverse ("library-based") cosine: every library peak is kept; the intensities of
#' all query peaks matching one library peak are summed. Matches the IUPAC reverse
#' search definition quoted in M&M and the legacy Merged/F_Dpc code path.
F_CosReverse <- function(it, ib, int_t, int_b) {
  nb <- length(int_b)
  if (!length(ib)) return(0)
  u <- numeric(nb)
  agg <- rowsum(int_t[it], group = ib, reorder = FALSE)
  u[as.integer(rownames(agg))] <- agg[, 1L]
  den <- sqrt(sum(u^2)) * sqrt(sum(int_b^2))
  if (den == 0) 0 else sum(u * int_b) / den
}

#' Forward / union cosine:
#' matched pairs contribute once per pair (so a query peak matching two library
#' peaks is counted twice), unmatched peaks on either side are padded with 0.
F_CosUnion <- function(it, ib, int_t, int_b) {
  num <- if (length(it)) sum(int_t[it] * int_b[ib]) else 0
  ss_t <- if (length(it)) sum(int_t[it]^2) else 0
  ss_b <- if (length(ib)) sum(int_b[ib]^2) else 0
  un_t <- if (length(it)) int_t[-unique(it)] else int_t
  un_b <- if (length(ib)) int_b[-unique(ib)] else int_b
  den <- sqrt(ss_t + sum(un_t^2)) * sqrt(ss_b + sum(un_b^2))
  if (den == 0) 0 else num / den
}
#' Candidate library spectra whose precursor is within +/- tol of q_precMZ. O(log n).
F_Candidates <- function(store, q_precMZ, tol) {
  if (is.na(q_precMZ)) return(seq_len(store$n))
  p  <- store$meta$precursor_mz
  i1 <- findInterval(q_precMZ - tol, p, left.open = TRUE) + 1L
  i2 <- findInterval(q_precMZ + tol, p)
  if (i2 < i1) integer(0) else i1:i2
}
#' Zero-score row for a query without any precursor-matched library candidate.
F_EmptyRow <- function(qm, par) {
  s0 <- list(dpc = 0, rdpc = 0, Match = 0L, MatchLib = 0L, decoy.mean = 0, decoy.mean.dpc = 0, xcorr = 0, Match.decoy = 0L)
  F_MakeRow(qm, NULL, NA, s0, par)
}

#  Spectrum stores -- read each database once, keep peaks as a ragged array
.getTbl <- function(con, candidates) {
  for (nm in candidates) if (DBI::dbExistsTable(con, nm)) return(nm)
  stop("none of these tables exist: ", paste(candidates, collapse = ", "))
}
#' Read a spectral SQLite database into memory.
#' Peaks are stored as three parallel vectors plus a per-spectrum index range,
#' which is far cheaper than a list of data.tables and cheap to export to workers.
F_LoadStore <- function(dbPth, pids = NA, pol = NA, spectraFilter = TRUE) {
  con <- DBI::dbConnect(RSQLite::SQLite(), dbPth)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  meta_tbl <- .getTbl(con, c("s_peak_meta", "library_spectra_meta"))
  peak_tbl <- .getTbl(con, c("s_peaks", "library_spectra"))
  meta <- as.data.table(DBI::dbGetQuery(con, sprintf("SELECT * FROM %s", meta_tbl)))
  peak <- as.data.table(DBI::dbGetQuery(con, sprintf("SELECT * FROM %s", peak_tbl)))
  
  id_col  <- if ("pid" %in% names(meta)) "pid" else "id"
  fk_col  <- if ("pid" %in% names(peak)) "pid" else "library_spectra_meta_id"
  int_col <- if ("i"   %in% names(peak)) "i"   else "intensity"
  
  if (!anyNA(pids))  meta <- meta[get(id_col) %in% pids]
  if (!is.na(pol) && "polarity" %in% names(meta))
    meta <- meta[tolower(polarity) == tolower(pol)]
  if (spectraFilter && "pass_flag" %in% names(peak)) peak <- peak[pass_flag == 1L | pass_flag == TRUE]
  
  setnames(meta, id_col, ".sid")
  setnames(peak, c(fk_col, int_col), c(".sid", ".int"))
  peak <- peak[.sid %in% meta$.sid]
  setorder(peak, .sid, mz)                    # sorted m/z within each spectrum: required
  
  # per-spectrum index range into the peak vectors
  idx <- peak[, .(start = .I[1], end = .I[.N]), by = .sid]
  meta <- merge(meta, idx, by = ".sid", all.x = FALSE, sort = FALSE)
  
  if (!"precursor_mz" %in% names(meta)) meta[, precursor_mz := NA_real_]
  meta[, precursor_mz := as.numeric(precursor_mz)]
  setorder(meta, precursor_mz)                # required for the binary search
  
  list(meta = meta,
       mz  = peak$mz,
       int = peak$.int,
       n   = nrow(meta))
}
#  Peak matching -- index based, no dense matrix
#' Match query peaks to library peaks within +/- w of each query peak.
#'
#' @param mz_t  query m/z (any order)
#' @param w_t   half-window per query peak (same length as mz_t)
#' @param mz_b  library m/z, MUST be sorted increasing
#' @return list(it, ib) of 1-based indices into mz_t / mz_b, one entry per matched pair
F_MatchPeaks <- function(mz_t, w_t, mz_b) {
  nb <- length(mz_b)
  if (nb == 0L || length(mz_t) == 0L) return(list(it = integer(0), ib = integer(0)))
  lo <- mz_t - w_t
  hi <- mz_t + w_t
  i1 <- findInterval(lo, mz_b, left.open = TRUE) + 1L   # first index with mz_b >= lo
  i2 <- findInterval(hi, mz_b)                          # last  index with mz_b <= hi
  n  <- i2 - i1 + 1L
  n[n < 0L] <- 0L
  if (!any(n > 0L)) return(list(it = integer(0), ib = integer(0)))
  keep <- which(n > 0L)
  ib   <- sequence(n[keep], from = i1[keep])            # base R >= 4.0
  it   <- rep.int(keep, n[keep])
  ok <- abs(mz_t[it] - mz_b[ib]) < w_t[it]
  list(it = it[ok], ib = ib[ok])
}



#' One output row for query qm against library row j (j = NA: no candidate).
F_MakeRow <- function(qm, lmeta, j, s, par) {
  lj <- function(col) if (anyNA(j) || is.null(lmeta[[col]])) NA else lmeta[[col]][j]
  data.table(
    dpc = s$dpc, rdpc = s$rdpc, Match = s$Match,
    decoy.mean = s$decoy.mean,
    xcorr = s$xcorr, Match.decoy = s$Match.decoy,
    rt_q = as.numeric(qm$retention_time %||% NA),
    rt_l = as.numeric(lj("retention_time")),
    pid_q = qm$.sid, pid_l = lj(".sid"),
    accession_q = qm$accession %||% NA, accession_l = lj("accession"),
    precursor_mz_q = as.numeric(qm$precursor_mz),
    precursor_mz_l = as.numeric(lj("precursor_mz")),
    inchikey_q = qm$inchikey_id %||% NA, inchikey_l = lj("inchikey_id"),
    inchi_q = NA, inchi_l = NA, computed_inchi_q = NA, computed_inchi_l = NA,
    smiles_q = qm$smiles %||% NA, smiles_l = lj("smiles"),
    computed_smiles_q = NA, computed_smiles_l = NA,
    splash_q = qm$splash %||% NA, splash_l = lj("splash"),
    entry_name_q = qm$name %||% NA, entry_name_l = lj("name"),
    precursor_type_q = qm$precursor_type %||% NA, precursor_type_l = lj("precursor_type"),
    instrument_type_q = qm$instrument_type %||% NA, instrument_type_l = lj("instrument_type"),
    instrument_q = qm$instrument %||% NA, instrument_l = lj("instrument"),
    collision_energy_q = qm$collision_energy %||% NA, collision_energy_l = lj("collision_energy"),
    resolution_q = qm$resolution %||% NA, resolution_l = lj("resolution"),
    mztol = par$mztol_label,
    MatchLib = s$MatchLib, SDRatio = par$SDRatio, DecoyCutoff = par$DecoyCutoff
  )
}

#  One query against the whole library
F_QueryOne <- function(qi, qstore, lstore, par) {
  
  qm <- qstore$meta[qi]
  ii <- qstore$int[qm$start:qm$end]
  mz <- qstore$mz [qm$start:qm$end]
  if (length(mz) == 0L || max(ii) <= 0) return(NULL)
  
  if (par$Round) mz <- round(mz, 5)
  ra <- ii / max(ii) * 100
  q  <- list(mz = mz, w = (mz^par$mzW) * (ra^par$raW))
  
  cand <- if (par$usePrecursors)
    F_Candidates(lstore, qm$precursor_mz, par$MS1Tol)
  else seq_len(lstore$n)
  
  if (length(cand) && !is.na(par$rttol) && "retention_time" %in% names(lstore$meta)) {
    keep <- abs(as.numeric(lstore$meta$retention_time[cand]) -
                  as.numeric(qm$retention_time)) < par$rttol
    cand <- cand[which(keep)]
  }
  # As "a query spectrum with no consensus/raw spectra
  # that match its precursor will obtain a cosine similarity score of 0 for both
  # target and decoy" and is labelled as an incorrect identification. Previously
  # such queries were silently dropped (return(NULL)), so they were missing from
  # the ROC/HOP curves and from the target/decoy score sets used for PEP/q-values.
  if (!length(cand)) {
    if (!isTRUE(par$keepUnmatched)) return(NULL)
    return(F_EmptyRow(qm, par))
  }
  # decoy spectra: built ONCE per query, not once per candidate
  d <- NULL
  if (par$decoy) {
    pmn <- F_CalPMMT(min(q$mz), instrument = par$instrument, MRP = par$MRPValue,
                     RefMZ = par$RefMZValue, MF = par$MF, SDRatio = par$SDRatio,
                     mztol = par$mztol)
    pmx <- F_CalPMMT(max(q$mz), instrument = par$instrument, MRP = par$MRPValue,
                     RefMZ = par$RefMZValue, MF = par$MF, SDRatio = par$SDRatio,
                     mztol = par$mztol)
    k   <- par$TolFactor
    tau <- c(seq(-par$HalfLength, -1L), seq(1L, par$HalfLength))
    off <- par$ShiftFactor * sign(tau) * k * pmx + tau * k * pmn
    dmz <- as.vector(outer(q$mz, off, "+"))
    dw  <- rep.int(q$w, length(off))
    dtau<- rep(tau, each = length(q$mz))
    ok  <- dmz > 0
    o   <- order(dmz[ok])                      # F_MatchPeaks needs nothing sorted on
    d   <- list(mz = dmz[ok][o], w = dw[ok][o], tau = dtau[ok][o])
  }
  # shared-peak prefilter replaces the precursor filter
  if (!par$usePrecursors && !is.null(lstore$pidx)) {
    cand <- cand[F_SharedPeakCandidates(lstore, q$mz, d, par)[cand]]
    if (!length(cand)) return(if (isTRUE(par$keepUnmatched)) F_EmptyRow(qm, par) else NULL)
  }
  lmeta <- lstore$meta
  sc <- vector("list", length(cand))
  for (k in seq_along(cand)) {
    j <- cand[k]
    sl <- lmeta$start[j]:lmeta$end[j]
    l_mz <- lstore$mz[sl]; l_ii <- lstore$int[sl]
    if (!length(l_mz) || max(l_ii) <= 0) next
    if (par$Round) l_mz <- round(l_mz, 5)
    l_ra <- l_ii / max(l_ii) * 100
    l_w  <- (l_mz^par$mzW) * (l_ra^par$raW)
    sc[[k]] <- F_ScorePair(q, l_mz, l_w, d, par)
  }
  ok <- !vapply(sc, is.null, logical(1))
  if (!any(ok)) return(if (isTRUE(par$keepUnmatched)) F_EmptyRow(qm, par) else NULL)
  
  S  <- rbindlist(sc[ok])
  jj <- cand[ok]
  if (!is.na(par$topN) && nrow(S) > par$topN) {
    o  <- head(order(-S[[par$topBy]]), par$topN)
    S  <- S[o]; jj <- jj[o]
  }
  F_MakeRow(qm, lmeta, jj, as.list(S), par)   # one vectorised data.table per query
}
#' Logical vector (length lstore$n): library spectra sharing >= MinMatch peaks with
#' the query (or with its pooled decoy). Same window/rounding as F_ScorePair, so the
#' counts equal Match / Match.decoy exactly; everything dropped would score 0.
F_SharedPeakCandidates <- function(lstore, q_mz, d, par) {
  pidx <- lstore$pidx
  mm   <- max(1L, par$MinMatch)
  tolw <- function(mz) par$TolFactor * F_CalPMMT(mz, instrument = par$instrument,
                                                 MRP = par$MRPValue, RefMZ = par$RefMZValue, MF = par$MF,
                                                 SDRatio = par$SDRatio, mztol = par$mztol)
  m   <- F_MatchPeaks(q_mz, tolw(q_mz), pidx$mz)
  hit <- tabulate(pidx$row[m$ib], nbins = lstore$n) >= mm
  if (!is.null(d)) {
    md  <- F_MatchPeaks(d$mz, tolw(d$mz), pidx$mz)
    hit <- hit | tabulate(pidx$row[md$ib], nbins = lstore$n) >= mm
  }
  hit
}

#  Scoring one query/library spectrum pair
#' @param q  list(mz, w) normalised+weighted query spectrum (w = weighted intensity)
#' @param d  NULL, or list(mz, w, tau) pooled decoy spectrum for this query
F_ScorePair <- function(q, l_mz, l_w, d, par) {
  w_t <- par$TolFactor * F_CalPMMT(q$mz, instrument = par$instrument, MRP = par$MRPValue,
                                   RefMZ = par$RefMZValue, MF = par$MF,
                                   SDRatio = par$SDRatio, mztol = par$mztol)
  m <- F_MatchPeaks(q$mz, w_t, l_mz)
  nmatch <- length(m$it)
  if (nmatch >= par$MinMatch) {
    dpc  <- F_CosUnion(m$it, m$ib, q$w, l_w)
    rdpc <- F_CosReverse(m$it, m$ib, q$w, l_w)
  } else {
    dpc <- 0; rdpc <- 0
  }
  
  out <- list(dpc = dpc, rdpc = rdpc, Match = nmatch,
              MatchLib = length(unique(m$ib)),
              decoy.mean = 0, xcorr = 0, Match.decoy = 0L)
  
  if (!is.null(d)) {
    w_d <- par$TolFactor * F_CalPMMT(d$mz, instrument = par$instrument, MRP = par$MRPValue,
                                     RefMZ = par$RefMZValue, MF = par$MF,
                                     SDRatio = par$SDRatio, mztol = par$mztol)
    md <- F_MatchPeaks(d$mz, w_d, l_mz)
    out$Match.decoy <- length(md$it)
    
    if (out$Match.decoy >= par$MinMatch) {
      #  mean over the individual shifts -- what Eq. (6)-(7) actually specifies.
      #     One aggregation gives all 2*HalfLength cosines at once.
      normv <- sqrt(sum(l_w^2))
      if (normv > 0) {
        DT <- data.table(tau = d$tau[md$it], ib = md$ib, u = d$w[md$it])
        DT <- DT[, .(u = sum(u)), by = .(tau, ib)]
        DT[, v := l_w[ib]]
        per <- DT[, .(num = sum(u * v), ss = sum(u^2)), by = tau]
        per[, cs := fifelse(ss > 0, num / (sqrt(ss) * normv), 0)]
        # shifts that matched nothing contribute a cosine of 0
        out$decoy.mean <- sum(per$cs) / (2 * par$HalfLength)
        # out$xcorr <- rdpc - out$decoy.mean
      }
    }
  }
  out
}
#' Public entry point of Spectral matching for LC-MS/MS datasets.
#' Scores: dpc (forward/union cosine), rdpc (reverse cosine).
#' decoy.mean: mean reverse cosine over the individual mass shifts (Eq. 6-7);
#' xcorr: rdpc - decoy.mean (Eq. 7).
#'
#' @param TolFactor matching window is +/- TolFactor * PMMT. Default 1, thus +/-1 sigma criterion
#' @param MS1Tol    precursor window in Da (+/-). The default 0.06 Da.
#' @param MF        mass fluctuation (ppm). 1 ppm for Orbitrap, 2 ppm for Q-TOF.
#' @param DecoyCutoff relative intensity (% of base peak) below which query peaks are NOT used to build the shifted decoy spectra (XcorrCutoff = 1).
#'                  The target score always uses the complete query spectrum.
#' @param keepUnmatched TRUE -> queries without any precursor-matched candidate are
#'                  reported with target = decoy = 0 .

SpectralMatching2 <- function(q_dbPth, l_dbPth,
                              usePrecursors = FALSE, MS1Tol = 0.06,
                              mztol = NA, SDRatio = 1, MF = 1,TolFactor = 1,
                              MRPValue = 17500, RefMZValue = 200,
                              cores = 1, topN = 10, topBy = c("dpc", "rdpc"),
                              Round = TRUE, MinMatch = 2,
                              instrument = c("Orbitrap", "Qtof"),
                              decoy = FALSE, # TRUE to generate the decoy spectra
                              HalfLength = 75, ShiftFactor = 1.5, DecoyCutoff = 0,
                              keepUnmatched = TRUE,
                              q_pol = NA, l_pol = NA,
                              q_pids = NA, l_pids = NA,
                              raW = 1, mzW = 0, rttol = NA,
                              outDir = ".", tag = NULL, write = TRUE,
                              verbose = TRUE) {
  
  instrument <- match.arg(instrument)
  t0 <- proc.time()[["elapsed"]]
  topBy <- match.arg(topBy)
  if (verbose) message("Loading query database ...")
  qstore <- F_LoadStore(q_dbPth, pids = q_pids, pol = q_pol)
  if (verbose) message("Loading library database ...")
  lstore <- F_LoadStore(l_dbPth, pids = l_pids, pol = l_pol)
  if (!usePrecursors) {
    if (verbose) message("Building global peak index for precursor-free search ...")
    lstore$pidx <- F_BuildPeakIndex(lstore, Round)
    if (is.na(topN)) warning("usePrecursors = FALSE without topN can produce a very large result")
  }
  if (verbose) message(sprintf("  %d query spectra, %d library spectra, %d library peaks",
                               qstore$n, lstore$n, length(lstore$mz)))
  
  par <- list(MS1Tol = MS1Tol, SDRatio = SDRatio, usePrecursors = usePrecursors,
              DecoyCutoff = DecoyCutoff, keepUnmatched = keepUnmatched,
              raW = raW, mzW = mzW, rttol = rttol, mztol = mztol,
              mztol_label = if (is.na(mztol)) "NA" else as.character(mztol),
              MRPValue = MRPValue, RefMZValue = RefMZValue, MF = MF,
              instrument = instrument, decoy = decoy,
              HalfLength = HalfLength, TolFactor = TolFactor, ShiftFactor = ShiftFactor,
              Round = Round, MinMatch = MinMatch, topN = topN, topBy = topBy)
  
  qi_all <- seq_len(qstore$n)
  
  if (cores > 1) {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    doParallel::registerDoParallel(cl)
    parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(data.table)))
    matched <- foreach(qi = qi_all, .packages = "data.table",
                       .export = c("F_QueryOne", "F_ScorePair", "F_SharedPeakCandidates", "F_MatchPeaks",
                                   "F_CosReverse", "F_CosUnion", "F_CalPMMT",
                                   "F_Candidates", "F_MakeRow", "F_EmptyRow", "%||%"),
                       .noexport = character(0)) %dopar%
      F_QueryOne(qi, qstore, lstore, par)
  } else {
    matched <- vector("list", length(qi_all))
    for (k in qi_all) {
      matched[[k]] <- F_QueryOne(k, qstore, lstore, par)
      if (verbose && k %% 100 == 0) message(sprintf("  %d / %d queries", k, qstore$n))
    }
  }
  
  matched <- rbindlist(matched[!vapply(matched, is.null, logical(1))],
                       use.names = TRUE, fill = TRUE)
  if (!nrow(matched)) { message("No matches found"); return(NULL) }
  setorder(matched, -dpc)                      # numeric now, so this is a numeric sort
  matched$xcorr <- matched$rdpc - matched$decoy.mean
  elapsed <- proc.time()[["elapsed"]] - t0
  if (verbose) message(sprintf("Finished in %.1f s (%.3f s per query spectrum)",
                               elapsed, elapsed / qstore$n))
  
  if (write) {
    if (!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)
    nm <- sprintf("%s%s_usePrecursors%s_MS1Tol%s_MF%s_TolFactor%s_MinMatch%s_decoy%s_DecoyCutoff%d.csv",
                  if (is.null(tag)) "" else paste0(tag, "_"),
                  F_TolLabel(mztol, SDRatio, instrument),usePrecursors, MS1Tol,MF,TolFactor,MinMatch, decoy,  DecoyCutoff)
    data.table::fwrite(matched, file.path(outDir, nm))
  }
  matched[]
}
#' Human-readable label for a tolerance setting; used in filenames and output.
F_TolLabel <- function(mztol, SDRatio = 1, instrument = "Orbitrap") {
  if (is.character(mztol)) mztol <- suppressWarnings(as.numeric(mztol))
  if (is.na(mztol)) sprintf("dynamic_%s_SD%s", instrument, SDRatio)
  else if (mztol > 1) sprintf("ppm%s", mztol)
  else sprintf("Da%s", mztol)
}


`%||%` <- function(a, b) if (is.null(a)) b else a



