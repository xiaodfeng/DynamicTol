# This file contains functions used for spectral matching
# The functions below are modified based on the spectralMatching.R file in the msPurity package with the following changes:
# 1) Added the parameter mztol to enable the dynamic matching
# 2) Added the parameter HalfLenth to enable creating the decoy MS/MS based on the query MS/MS and thus the calculation of decoy score.
# 3) replace the align() function with F_DynamicMatching() or F_FixedMatching() to speed up the matching process between query and library
# 4) Use reverse dot product cosine only to keep equal vector length
# 5) Merge the results in case there are multiple hits per library by using the function of F_MergeTarget()

#' @title Spectral matching for LC-MS/MS datasets
#' @aliases spectralMatching
#' @import data.table
#' @description
#' **General**
#'
#' Perform spectral matching to spectral libraries for an LC-MS/MS dataset.
#'
#' The spectral matching is performed from a **Query** SQLite spectral-database against a **Library** SQLite spectral-database.
#'
#' The query spectral-database in most cases should contain be the "unknown" spectra database
#'
#' The library spectral-database in most cases should contain the "known" spectra from either public or user generated resources.
#' The library SQLite database by default contains data from MoNA including Massbank, HMDB, LipidBlast and GNPS.
#' A larger database can be downloaded from [here](https://github.com/computational-metabolomics/msp2db/releases).
#' To create a user generated library SQLite database the following tool can be used to generate a SQLite database
#' from a collection of MSP files: [msp2db](https://github.com/computational-metabolomics/msp2db/releases).
#' It should be noted though, that as long as the schema of the spectral-database is as described above, then any database can be used
#' for either the library or query -  even allowing for the same database to be used.
#'
#' The spectral matching functionality has four main components, spectral filtering, spectral alignment, spectral matching,
#' and summarising the results.
#'
#' Spectral filtering is simply filtering both the library and query spectra to be search against (e.g. choosing
#' the library source, instrument, retention time, precursor PPM tolerance etc).
#'
#' The spectral alignment stage involves aligning the query peaks to the library peaks. The approach used is similar
#' to modified pMatch algorithm described in Zhou et al 2015.
#'
#' The spectral matching of the aligned spectra is performed against a combined intensity and m/z weighted vector - created for both
#' the query and library spectra (wq and wl). See below:
#'
#' \deqn{w=intensity^x * mz^y}
#'
#' Where x and y represent weight factors, defaults to *x*=0.5 and *y*=2 as per MassBank. These can be adjusted by
#' the user though.
#'
#' The aligned weighted vectors are then matched using reverse dot product cosine with the dynamic matching tolerance.
#' See below for dot product cosine equation.
#'
#' \deqn{dpc =  wq * wl / \sqrt{\sum wq^2} * \sqrt{\sum wl^2}}
#'
#' See the vigenttes for more details regarding matching algorithms used.
#'
#' @param q_dbPth character; Path of the database of queries that will be searched against the library spectra. Generated from createDatabase
#' @param l_dbPth character; path to library spectral SQLite database. Defaults to msPurityData package data.
#' @param q_ppmPrec numeric; ppm tolerance for query precursor
#' @param q_pol character; Polarity of query spectra ('positive', 'negative', NA).
#' @param q_instrumentTypes vector; Instrument types for query spectra.
#' @param q_instruments vector; Instruments for query spectra (note that this is used in combination with q_instrumentTypes - any
#'                              spectra matching either q_instrumentTypes or q_instruments will be used).
#' @param q_sources vector; Sources of query spectra (e.g. massbank, hmdb).
#' @param q_pids vector; pids for query spectra (correspond to column 'pid; in s_peak_meta)
#' @param q_rtrange vector; retention time range (in secs) of query spectra, first value mininum time and second value max e.g. c(0, 10) is between 0 and 10 seconds
#' @param q_accessions vector; accession ids to filter query spectra
#' @param l_ppmPrec numeric; ppm tolerance for library precursor
#' @param l_pol character; Polarity of library spectra ('positive', 'negative', NA)
#' @param l_instrumentTypes vector; Instrument types for library spectra.
#' @param l_instruments vector; Instruments for library spectra (note that this is used in combination with q_instrumentTypes - any
#'                              spectra matching either q_instrumentTypes or q_instruments will be used).
#' @param l_sources vector; Sources of library spectra (e.g. massbank, hmdb).
#' @param l_pids vector; pids for library spectra (correspond to column 'pid; in s_peak_meta)
#' @param l_rtrange vector; retention time range (in secs) of library spectra, first value mininum time and second value max e.g. c(0, 10) is between 0 and 10 seconds
#' @param l_accessions vector; accession ids to filter library spectra
#'
#' @param usePrecursors boolean; If TRUE spectra will be filtered by similarity of precursors based on ppm range defined by l_ppmPrec and q_ppmPrec
#' @param raW numeric; Relative abundance weight for spectra (default to 0.5 as determined by massbank for ESI data)
#' @param mzW numeric; mz weight for spectra (default to 2 as determined by massbank for ESI data)
#' @param rttol numeric ; Tolerance in time range between the library and query spectra retention time
#'
#' @param cores numeric; Number of cores to use
#'
#' @param q_dbType character; Query database type for compound database can be either (sqlite, postgres or mysql)
#' @param l_dbType character; Library database type for compound database can be either (sqlite, postgres or mysql)
#'
#' @return Returns a list containing the following elements
#'
#' **q_dbPth**
#'
#' Path of the query database

#' **matchedResults**
#'
#' All matched results from the query spectra to the library spectra. Contains the same columns as above
#' but without the XCMS details. This table is useful to observe spectral matching results
#' for all MS/MS spectra irrespective of if they are linked to XCMS MS1 features.
#'
#' @return list of database details and dataframe summarising the results for the xcms features
#'
#' @examples
#' \dontrun{
#' #====== XCMS =================================
#' ## Read in MS data
#' #msmsPths <- list.files(system.file("extdata", "lcms", "mzML",
#' #           package="msPurityData"), full.names = TRUE, pattern = "MSMS")
#' #ms_data = readMSData(msmsPths, mode = 'onDisk', msLevel. = 1)
#'
#' ## Find peaks in each file
#' #cwp <- CentWaveParam(snthresh = 5, noise = 100, ppm = 10, peakwidth = c(3, 30))
#' #xcmsObj  <- xcms::findChromPeaks(ms_data, param = cwp)
#'
#' ## Optionally adjust retention time
#' #xcmsObj  <- adjustRtime(xcmsObj , param = ObiwarpParam(binSize = 0.6))
#'
#' ## Group features across samples
#' #pdp <- PeakDensityParam(sampleGroups = c(1, 1), minFraction = 0, bw = 30)
#' #xcmsObj <- groupChromPeaks(xcmsObj , param = pdp)
#'
#' #====== dynamic ============================
#' #pa  <- purityA(msmsPths)
#' #pa <- frag4feature(pa = pa, xcmsObj = xcmsObj)
#' #pa <- filterFragSpectra(pa, allfrag=TRUE)
#' #pa <- averageAllFragSpectra(pa)
#' #q_dbPth <- createDatabase(pa, xcmsObj, metadata=list('polarity'='positive','instrument'='Q-Exactive'))
#' #sm_result <- spectralMatching(q_dbPth, cores=4, l_pol='positive')
#'
#' td <- tempdir()
#' q_dbPth <- system.file("extdata", "tests", "db", "createDatabase_example.sqlite", package="dynamic")
#'
#' rid <- paste0(paste0(sample(LETTERS, 5, TRUE), collapse=""),  paste0(sample(9999, 1, TRUE), collapse=""), ".sqlite")
#' sm_out_pth <- file.path(td, rid)
#'
#' result <- spectralMatching(q_dbPth, cores=1, l_accessions = c('PR100407', 'ML005101', 'CCMSLIB00003740024'))
#'
#' }
#' @md
#' @export
SpectralMatching <- function(q_dbPth, l_dbPth,
                             q_ppmPrec = 20, l_ppmPrec = 20,
                             q_pol = NA, l_pol = NA,
                             q_instrumentTypes = NA, l_instrumentTypes = NA,
                             q_instruments = NA, l_instruments = NA,
                             q_sources = NA, l_sources = NA,
                             q_pids = NA, l_pids = NA,
                             q_rtrange = c(NA, NA), l_rtrange = c(NA, NA),
                             q_accessions = NA, l_accessions = NA,
                             q_dbType = "sqlite", l_dbType = "sqlite",
                             usePrecursors = TRUE,
                             raW = 1, mzW = 0,
                             rttol = NA, mztol = NA,
                             MRPValue = 17500, RefMZValue = 200,
                             libsch = 'reverse', # reverse for reverse library search, and forward for forward library search
                             decoy = FALSE, # TRUE to generate the decoy spectra
                             cores = 1, HalfLenth= 75) {

  tictoc::tic('library search total running time')
  if (is.na(l_dbPth)) {
    l_dbPth <- system.file("extdata", "library_spectra", "library_spectra.db", package = "msPurityData")
  }
  tictoc::tic("Filter query and library dataset")
  q_con <- DBI::dbConnect(RSQLite::SQLite(), q_dbPth)
  q_speakmeta <- filterSMeta(
    pol = q_pol,
    instrumentTypes = q_instrumentTypes,
    instruments = q_instruments,
    sources = q_sources,
    pids = q_pids,
    rtrange = q_rtrange,
    con = q_con,
    accessions = q_accessions
  )
  l_con <- DBI::dbConnect(RSQLite::SQLite(), l_dbPth)
  l_speakmeta <- filterSMeta(
    pol = l_pol,
    instrumentTypes = l_instrumentTypes,
    instruments = l_instruments,
    sources = l_sources,
    pids = l_pids,
    rtrange = l_rtrange,
    con = l_con,
    accessions = l_accessions
  )
  q_fpids <- pullPid(q_speakmeta) # extract the query ids
  l_fpids <- pullPid(l_speakmeta)
  # Loop through the query dataset and spectra match against the library spectra
  dbDetails <- list(
    "q" = list(pth = q_dbPth, type = q_dbType),
    "l" = list(pth = l_dbPth, type = l_dbType),
    q_ppmPrec=q_ppmPrec, l_ppmPrec=l_ppmPrec, HalfLenth=HalfLenth,
    usePrecursors=usePrecursors, mzW=mzW, raW=raW, rttol=rttol,
    mztol=mztol,MRPValue=MRPValue,RefMZValue=RefMZValue,libsch=libsch,decoy=decoy
  )
  tictoc::toc()
  tictoc::tic("aligning and matching")
  # Check cores and choose if parallel or not (do or dopar)
  if (cores > 1) {
    cl <- parallel::makeCluster(cores, type = "SOCK")
    doSNOW::registerDoSNOW(cl)
    parallel::clusterExport(cl, c('F_CalPMMT','alignAndMatch','F_DynamicMatching','filterSMeta','filterPrecursors','F_FixedMatchingPPM','F_FixedMatching','getScanPeaksSqlite',
                                  'getSmeta','getMetaCols','getPeakCols','MSsim','F_MergeDecoy','F_MergeTarget','pullPid','queryVlibrary'))
    matched <- foreach(q_fpid = q_fpids,.packages = c('data.table','dplyr','foreach')) %dopar% queryVlibrary(q_pid=q_fpid, l_pids=l_fpids, dbDetails=dbDetails)
    parallel::stopCluster(cl)
  } else {
    matched <- foreach(q_fpid = q_fpids,.packages = c('data.table','dplyr','foreach')) %do% queryVlibrary(q_pid=q_fpid, l_pids=l_fpids, dbDetails=dbDetails)
  }
  # View(matched)
  # print(matched)
  matched <- bind_rows(matched) # Combine the results from lists
  if (nrow(matched) == 0) {
    message("No matches found")
    return(NULL)
  }
  tictoc::toc()
  tictoc::tic("Output results")
  # matched <- matched[, !names(matched) == "X1"]   # remove the plyr id column
  # ensure numeric
  # nmCols <- c("dpc","dpc.decoy","dpc.xcorr","mcount", "allcount", "mpercent")
  # matched[,nmCols] <- as.numeric(as.character(unlist(matched[,nmCols])))
  # make sure all NA values are fully NA values
  # Filter out the candidates with negative Xcorr
  # matched <- matched[(dpD.xcorr>0)]
  # Filter out the candidates with low number of match
  # matched$Match <- as.numeric(matched$Match)
  # matched <- matched[(Match>=2)]
  # Order the results with the dpc decreasing, while dpc.decoy increasing
  # setorder(matched, -dpc, dpc.decoy)
  # matched <- matched[order(-as.numeric(matched$dpc), as.numeric(matched$dpc)), ]
  # Add index of matched candidates
  # matched$mid <- 1:nrow(matched)
  # Add rtdiff
  # matched$rtdiff <- as.numeric(matched$library_rt) - as.numeric(matched$query_rt)
  matched <- data.table(matched) # transfer into the data.table format
  name <- paste0('mztol_',mztol,'_libsch_',libsch,'_decoy_',decoy, '.csv')
  # write.csv(matched, name,row.names = FALSE) #
  data.table::fwrite(matched, name,row.names = FALSE) # faster method to write the csv files.
  tictoc::toc()
  on.exit(DBI::dbDisconnect(q_con))
  on.exit(DBI::dbDisconnect(l_con))
  tictoc::toc()
  return(matched)
}
#' @title F_CalPMMT
#' @description
#' Calculate the peak matching mass tolerance (PMMT)
#' according to the mass resolving power (MRP) and reference m/z (RefMZ)
#' @export
F_CalPMMT <- function(mz, MRP = 17500, RefMZ = 200) {
  A <- 1 / (MRP * (RefMZ^0.5))
  B <- A / 2.35482
  return(B * mz^1.5 + mz / 1000000)
}
alignAndMatch <- function(q_speaksi, l_speaksi, dbDetails) {
  # For testing purpose begin
  # top <-data.frame("mz" = c(86.09695,486.42896,610.41681,669.48749,
  #                           104.10786,114.41492,131.05865,146.98178,230.20305),
  #                   "intensity" = c(83,100,79,27,
  #                                   5,4,6,15,20))
  # bottom <-data.frame("mz" = c(86.09702,486.43140,610.42194,669.49377,
  #                              143.80835,190.69144,243.68671,318.04742,399.69229),
  #                     "intensity" = c(38,85,100,43,
  #                                     24,20,5,4,6))
  # raW <- 1
  # mzW <- 0
  # mztol <- NA
  # MS2_18374[,ra:=i]
  # MS2_18373[,ra:=i]
  # q_speaksi<-MS2_18373
  #
  # l_speaksi<-MS2_18374
  # For testing purpose end
  # message('Normalization and scaling')
  q_speaksi <- as.data.frame(q_speaksi) # transfer from tibble into data frame
  l_speaksi <- as.data.frame(l_speaksi)
  q_speaksi$ra <- (q_speaksi$i / max(q_speaksi$i)) * 100 # Normalization
  l_speaksi$ra <- (l_speaksi$i / max(l_speaksi$i)) * 100
  q_speaksi$w <- (q_speaksi$mz^dbDetails$mzW) * (q_speaksi$ra^dbDetails$raW)   # Scaling factor
  l_speaksi$w <- (l_speaksi$mz^dbDetails$mzW) * (l_speaksi$ra^dbDetails$raW)
  top <- data.frame('mz'=q_speaksi$mz,'intensity'=q_speaksi$w)
  bottom <- data.frame('mz'=l_speaksi$mz,'intensity'=l_speaksi$w)
  # message('Target search to obtain dpc')
  # tictoc::tic()
  if (is.na(dbDetails$mztol)) {     ## Dynamic matching
    Align <- F_DynamicMatching(top = top, bottom = bottom, MRP= dbDetails$MRPValue,RefMZ = dbDetails$RefMZValue) #
  } else if (dbDetails$mztol <= 1) {     ## Fixed matching with Da unit
    Align <- F_FixedMatching(top = top, bottom = bottom, mztol = dbDetails$mztol)
  } else {     ## Fixed matching with PPM unit
    Align <- F_FixedMatchingPPM(top = top, bottom = bottom, mztol = dbDetails$mztol)
  }
  # tictoc::toc()
  Match <- nrow(Align[(Align[, 3] > 0)&(Align[, 2] > 0), ]) # check the number of matches
  if (dbDetails$libsch == 'reverse') {Align <- Align[Align[, 3] > 0, ]}
  if (Match>0){
    MergeAlign <- F_MergeTarget(Align)  # Merge the results in case there are multiple hits per library
    dpc <- MSsim(MergeAlign)
  } else {
    dpc <- 0
  }
  # message('Decoy search to obtain dpc.decoy and dpc.xcorr')
  dpc.decoy <- 0
  Match.decoy <- 0
  if (dbDetails$decoy){
    PMMTmin<-F_CalPMMT(min(top$mz),MRP= dbDetails$MRPValue,RefMZ = dbDetails$RefMZValue) #
    PMMTmax<-F_CalPMMT(max(top$mz),MRP= dbDetails$MRPValue,RefMZ = dbDetails$RefMZValue) #
    # print(paste0('MRPValue ',MRPValue,' RefMZValue ',RefMZValue))
    ## create the decoy
    decoy <-  top
    decoy.rbind <- data.table()
    if (is.na(dbDetails$mztol)) {     ## Dynamic matching
      for (id in c(seq(-dbDetails$HalfLenth, -1), seq(1, dbDetails$HalfLenth))) {
        decoy$mz <- top$mz + 1.5 * (id / abs(id)) * PMMTmax + id * PMMTmin
        decoy.rbind <- rbind(decoy.rbind, decoy)
      }
      decoy.rbind <- decoy.rbind [mz>0] %>% setorder(.,mz) # exclude the entries with mz less than zero
      Align.decoy <- F_DynamicMatching(top = as.data.frame(decoy.rbind), bottom = bottom, MRP= dbDetails$MRPValue,RefMZ = dbDetails$RefMZValue)
    } else { ## Fixed matching
      for (id in c(seq(-dbDetails$HalfLenth, -1), seq(1, dbDetails$HalfLenth))) {
        decoy$mz <- top$mz + 1.5 * (id / abs(id)) * PMMTmax + id * PMMTmin
        decoy.rbind <- rbind(decoy.rbind, decoy)
      }
      decoy.rbind <- decoy.rbind [mz>0] %>% setorder(.,mz) # exclude the entries with mz less than zero
      Align.decoy <- F_FixedMatching(top = as.data.frame(decoy.rbind), bottom = bottom)
    }
    Match.decoy<-nrow(Align.decoy[(Align.decoy[, 3] > 0)&(Align.decoy[, 2] > 0), ])
    if (dbDetails$libsch == 'reverse') {Align.decoy <- Align.decoy[Align.decoy[, 3] > 0, ]}
    if (Match.decoy > 0){
      MergeAlign.decoy<-F_MergeDecoy(Align.decoy)
      dpc.decoy <- MSsim(MergeAlign.decoy)
    }
  }
  dpc.xcorr <- dpc-dpc.decoy
  ## Echo
  # print(paste('dpc',dpc,'dpc.decoy',dpc.decoy))
  # print(mztol)
  # print('top')
  # print(top)
  # print('bottom')
  # print(bottom)
  # print('Align')
  # print(Align)
  # print('MergeAlign')
  # print(MergeAlign)
  # print('MergeAlign.decoy')
  # print(MergeAlign.decoy)
  ## Return results
  return(c(
    "dpc" = dpc, "dpc.decoy" = dpc.decoy, "dpc.xcorr" = dpc.xcorr,
    "Match" = Match,  "Match.decoy" = Match.decoy
  ))
}


#' @title F_DynamicMatching
#' @description
#' Matching between query and library using the dynamic tolerance
#' @export
F_DynamicMatching <- function(top = top, bottom = bottom, RefMZ=200, MRP = 17500) {## Dynamic matching
  A <- 1 / (MRP * (RefMZ^0.5))
  B <- A / 2.35482
  for (i in 1:nrow(bottom)) {
    top[, 1][abs(bottom[, 1][i] - top[, 1]) < B * bottom[, 1][i]^1.5 + bottom[, 1][i] / 1000000] <- bottom[, 1][i]
  }
  alignment <- merge(top, bottom, by = 1, all = TRUE)
  alignment[, c(2, 3)][is.na(alignment[, c(2, 3)])] <- 0 # convert NAs to zero (R-Help, Sept. 15, 2004, John Fox)
  names(alignment) <- c("mz", "intensity.top", "intensity.bottom")
  # print(alignment)
  return(alignment)
}


filterSMeta <- function(raThres = 0, pol = NA, instrumentTypes = NA, instruments = NA, sources = NA, xcmsGroups = NA,
                        pids = NA, rtrange = c(NA, NA), spectraTypes = NA, accessions = NA, con) {
  # get column names
  # PRAGMA table_info();
  meta_cn <- getMetaCols(con)
  speakmeta <- getSmeta(con, pids)
  if ("accession" %in% meta_cn$name && !anyNA(accessions)) {
    speakmeta <- speakmeta %>% dplyr::filter(accession %in% accessions)
  }
  if ("polarity" %in% meta_cn$name && !is.na(pol)) {
    speakmeta <- speakmeta %>% dplyr::filter(lower(polarity) == lower(pol))
  }
  if ("instrument_type" %in% meta_cn$name && "instrument" %in% meta_cn$name && !anyNA(instrumentTypes) && !anyNA(instruments)) {
    speakmeta <- speakmeta %>% dplyr::filter(instrument_type %in% instrumentTypes || instrument %in% instruments)
  } else if ("instrument_type" %in% meta_cn$name && !anyNA(instrumentTypes)) {
    speakmeta <- speakmeta %>% dplyr::filter(instrument_type %in% instrumentTypes)
  } else if ("instrument" %in% meta_cn$name && !anyNA(instruments)) {
    speakmeta <- speakmeta %>% dplyr::filter(instrument %in% instruments)
  }
  if (!anyNA(sources)) {
    if (DBI::dbExistsTable(con, "library_spectra_source")) {
      sourcetbl <- con %>% dplyr::tbl("library_spectra_source")
      speakmeta <- speakmeta %>%
        dplyr::left_join(sourcetbl, by = c("library_spectra_source_id" = "id"), suffix = c("", ".y")) %>%
        dplyr::filter(name.y %in% sources)
    } else if (DBI::dbExistsTable(con, "source")) {
      sourcetbl <- con %>% dplyr::tbl("source")
      speakmeta <- speakmeta %>%
        dplyr::left_join(sourcetbl, by = c("sourceid" = "id"), suffix = c("", ".y")) %>%
        dplyr::filter(name.y %in% sources)
    }
  }
  if ("retention_time" %in% meta_cn$name && !anyNA(rtrange)) {
    speakmeta <- speakmeta %>% dplyr::filter(retention_time >= rtrange[1] & retention_time <= rtrange[2])
  }
  if ("spectrum_type" %in% meta_cn$name && !anyNA(spectraTypes)) {
    if ("av_all" %in% spectraTypes) {
      spectraTypes[spectraTypes == "av_all"] <- "all"
    }
    speakmeta <- speakmeta %>% dplyr::filter(spectrum_type %in% spectraTypes)
  }
  # print("speakmeta")
  return(speakmeta)
}


filterPrecursors <- function(q_pid, q_speakmeta, l_speakmeta, q_ppmPrec, l_ppmPrec) {
  return(l_speakmetaFiltered)
}

F_FixedMatching <- function(top = top, bottom = bottom, mztol = 0.005) {## Fixed matching with Da unit
  for (i in 1:nrow(bottom)) {
    top[, 1][abs(bottom[, 1][i] - top[, 1]) < mztol] <- bottom[, 1][i]
  }
  alignment <- merge(top, bottom, by = 1, all = TRUE) # use the bottom (library) as reference
  alignment[, c(2, 3)][is.na(alignment[, c(2, 3)])] <- 0 # convert NAs to zero (R-Help, Sept. 15, 2004, John Fox)
  names(alignment) <- c("mz", "intensity.top", "intensity.bottom")
  # print(alignment)
  return(alignment)
}
F_FixedMatchingPPM <- function(top = top, bottom = bottom, mztol = 5) {## Fixed matching with PPM unit
  for (i in 1:nrow(bottom)) {
    top[, 1][abs(bottom[, 1][i] - top[, 1]) < mztol * bottom[, 1][i] / 1000000] <- bottom[, 1][i]
  }
  alignment <- merge(top, bottom, by = 1, all = TRUE) # use the bottom (library) as reference
  alignment[, c(2, 3)][is.na(alignment[, c(2, 3)])] <- 0 # convert NAs to zero (R-Help, Sept. 15, 2004, John Fox)
  names(alignment) <- c("mz", "intensity.top", "intensity.bottom")
  # print(alignment)
  return(alignment)
}
getScanPeaksSqlite <- function(con, spectraFilter = TRUE, spectraTypes = NA, raThres = NA, pids = NA) {
  if (DBI::dbExistsTable(con, "s_peaks")) {
    speaks <- con %>% dplyr::tbl("s_peaks")
  } else if (DBI::dbExistsTable(con, "library_spectra")) {
    # old sqlite format
    speaks <- con %>% dplyr::tbl("library_spectra")
  } else {
    stop("No spectra available")
  }
  cn <- getPeakCols(con)
  if ("pid" %in% cn$name && !anyNA(pids)) {
    speaks <- speaks %>% dplyr::filter(pid %in% pids)
  } else if ("library_spectra_meta_id" %in% cn$name && !anyNA(pids)) {
    speaks <- speaks %>% dplyr::filter(library_spectra_meta_id %in% pids)
  }
  if ("pass_flag" %in% cn$name && spectraFilter) {
    speaks <- speaks %>% dplyr::filter(pass_flag == TRUE)
  }
  if ("type" %in% cn$name && !anyNA(spectraTypes)) {
    speaks <- speaks %>% dplyr::filter(type %in% spectraType)
  }
  if ("ra" %in% cn$name && !is.na(raThres)) {
    speaks <- speaks %>% dplyr::filter(ra > raThres)
  }
  return(speaks)
}


getSmeta <- function(con, pids = NA) {
  if (DBI::dbExistsTable(con, "s_peak_meta")) {
    speakmeta <- con %>% dplyr::tbl("s_peak_meta")
    if (!anyNA(pids)) {
      speakmeta <- speakmeta %>% dplyr::filter(pid %in% pids)
    }
  } else if (DBI::dbExistsTable(con, "library_spectra_meta")) {
    # old sqlite format
    speakmeta <- con %>% dplyr::tbl("library_spectra_meta")
    if (!anyNA(pids)) {
      speakmeta <- speakmeta %>% dplyr::filter(id %in% pids)
    }
  } else {
    stop("No meta data for spectra available")
  }
  return(speakmeta)
}

getMetaCols <- function(con) {
  if (DBI::dbExistsTable(con, "s_peak_meta")) {
    meta_cn <- DBI::dbGetQuery(con, "PRAGMA table_info(s_peak_meta)")
  } else {
    meta_cn <- DBI::dbGetQuery(con, "PRAGMA table_info(library_spectra_meta)")
  }
  return(meta_cn)
}
getPeakCols <- function(con) {
  if (DBI::dbExistsTable(con, "s_peak_meta")) {
    cn <- DBI::dbGetQuery(con, "PRAGMA table_info(s_peaks)")
  } else {
    cn <- DBI::dbGetQuery(con, "PRAGMA table_info(library_spectra)")
  }
  return(cn)
}

#' @title MSsim
#' @description
#' Similarity score calculation for a data frame
#' @export
MSsim <- function(alignment) { ## similarity score calculation
  u <- alignment[, 2]
  v <- alignment[, 3]
  return(as.vector((u %*% v) / (sqrt(sum(u^2)) * sqrt(sum(v^2)))))
}

F_MergeTarget <- function(Target){
  Target<-data.table(Target)
  TargetUni<-unique(Target,by='mz')
  Container<-data.table()
  for (i in 1:nrow(TargetUni)){
    mzUni<-TargetUni[i,]$mz
    DT<-data.table('mz'=mzUni,
                   'intensity.top'= sum(Target[mz==mzUni,]$intensity.top), # calculate the sum of the intensity top
                   'intensity.bottom'=TargetUni[i,]$intensity.bottom)
    Container<-rbind(Container,DT)
  }
  return(as.data.frame(Container))
}

F_MergeDecoy <- function(Decoy){
  Decoy<-data.table(Decoy)
  DecoyUni<-unique(Decoy,by='mz')
  Container<-data.table()
  for (i in 1:nrow(DecoyUni)){
    mzUni<-DecoyUni[i,]$mz
    DT<-data.table('mz'=mzUni,
                   'intensity.top'= mean(Decoy[mz==mzUni,]$intensity.top), # calculate the mean of the intensity top
                   'intensity.bottom'=DecoyUni[i,]$intensity.bottom)
    Container<-rbind(Container,DT)
  }
  return(as.data.frame(Container))
}

pullPid <- function(sp, pids) {
  tble <- sp %>% dplyr::collect()
  nms <- colnames(tble)
  if ("pid" %in% nms) {
    pids <- tble$pid
  } else {
    pids <- tble$id
  }
  return(pids)
}
queryVlibrary <- function(q_pid, l_pids, dbDetails) {
  # q_pid <- 5214
  print(paste(q_pid, 'mztol', dbDetails$mztol , 'libsch', dbDetails$libsch, 'decoy' , dbDetails$decoy))
  q_con <- DBI::dbConnect(RSQLite::SQLite(), dbDetails$q$pth)
  l_con <- DBI::dbConnect(RSQLite::SQLite(), dbDetails$l$pth)
  q_speakmetai <- getSmeta(q_con, q_pid) %>% dplyr::collect()
  # print(q_speakmetai)
  ## Here we extract the MS2 peaks according to the q_pid, important for calculating the dot product
  q_speaksi <- getScanPeaksSqlite(q_con, pids = q_pid) %>% dplyr::collect()
  # if no peaks, then skip
  if (nrow(q_speaksi) == 0) {
    return(NULL)
  }
  l_speakmeta <- getSmeta(l_con, l_pids) %>% dplyr::collect()
  l_speaks <- getScanPeaksSqlite(l_con, pids = l_pids) %>% dplyr::collect()
  if (dbDetails$usePrecursors) {
    q_precMZ <- q_speakmetai$precursor_mz
    # Check if precursors are within tolerance
    # We have ppm tolerances for both the library and the query
    q_precMZlo <- q_precMZ - ((q_precMZ * 0.000001) * dbDetails$q_ppmPrec)
    q_precMZhi <- q_precMZ + ((q_precMZ * 0.000001) * dbDetails$q_ppmPrec)
    # Search against the range for the library
    l_fspeakmeta <- l_speakmeta %>%
      dplyr::filter((q_precMZhi >= precursor_mz - ((precursor_mz * 0.000001) * dbDetails$l_ppmPrec)) &
                      (precursor_mz + ((precursor_mz * 0.000001) * dbDetails$l_ppmPrec) >= q_precMZlo)) %>%
      dplyr::collect()
  } else {
    l_fspeakmeta <- l_speakmeta %>% dplyr::collect()
  }
  if (!is.na(dbDetails$rttol)) {
    l_fspeakmeta <- l_fspeakmeta %>% dplyr::filter(abs(retention_time - q_speakmetai$retention_time) < rttol)
  }
  if (nrow(l_fspeakmeta) == 0) {
    return(NULL)
  }
  l_fpids <- l_fspeakmeta$id

  ## Use this loop to match each library MS2 to the selected query MS2 without parallel
  searched <- foreach(l_fpid = l_fpids, .packages = c('data.table', 'dplyr')) %do% {
    l_speaksi <- data.frame(l_speaks %>% dplyr::filter(library_spectra_meta_id == l_fpid) %>% dplyr::collect())
    l_speakmetai <- data.frame(l_speakmeta %>% dplyr::filter(id == l_fpid) %>% dplyr::collect())
    am <- alignAndMatch(q_speaksi,l_speaksi,dbDetails)
    return(c(am,
             "query_rt" = as.numeric(q_speakmetai$retention_time),
             "library_rt" = as.numeric(l_speakmetai$retention_time),
             "query_qpid" = q_pid,
             "library_lpid" = l_speakmetai$id,
             "query_accession" = q_speakmetai$accession,
             "library_accession" = l_speakmetai$accession,
             "query_precursor_mz" = q_speakmetai$precursor_mz,
             "library_precursor_mz" = l_speakmetai$precursor_mz,
             "query_inchikey" = q_speakmetai$inchikey_id,
             "library_inchikey" = l_speakmetai$inchikey_id,
             "query_entry_name" = q_speakmetai$name,
             "library_entry_name" = l_speakmetai$name,
             "library_precursor_type" = l_speakmetai$precursor_type,
             "library_instrument_type" = l_speakmetai$instrument_type,
             "mztol" = dbDetails$mztol,
             "libsch" = dbDetails$libsch
    ))
  }
  on.exit(DBI::dbDisconnect(q_con))
  on.exit(DBI::dbDisconnect(l_con))
  return(bind_rows(searched))
}







