# This file contains functions used for the spectra processing
# The functions are sorted alphabetically

F_CalPMMT <- function(mz, MRP = 17500, RefMZ = 200) {
  A <- 1 / (MRP * (RefMZ^0.5))
  B <- A / 2.35482
  return(B * mz^1.5 + mz / 1000000)
}

F_CountCutoff <- function(DT) {
  setorder(DT, -dpc)
  CountCutoff.rbind <- data.table()
  for (id in seq(0.01, 1, 0.01)) {
    CountCutoff <- data.table("Cutoff" = id, "Count" = nrow(DT[dpc >= id]))
    CountCutoff.rbind <- rbind(CountCutoff.rbind, CountCutoff)
  }
  return(CountCutoff.rbind)
}
F_CountSpectraPerInchikey <- function(Inchikey) {
  print(paste0(Inchikey, nrow(Meta[inchikey_14_precursor_type == Inchikey])))
  return(nrow(Meta[inchikey_14_precursor_type == Inchikey]))
}
F_CountPeaksPerSpectra<-function(z){
  Len<-nrow(library_spectra[library_spectra_meta_id==z])
  print(paste(z,Len))
  return(Len)
  # Test<-TarDecUni[1:2]
  # Test[,SpectraLength:=sapply(Test$library_lpid,
  #                             F_CountSpectraLength)]
}
F_CountIdentical <- function(DT) {
  DT <- unique(DT, by = c("query_inchikey14", "library_inchikey14", "library_precursor_mz"))
  Uni <- max(nrow(unique(DT, by = "query_accession")), nrow(unique(DT, by = "library_accession")))
  print(paste("Identical pairs", nrow(DT) / 2, "Spectra number", Uni / 2))
}
F_CountCandidates<-function(z){
  Len<-nrow(TarDec[query_qpid==z])
  print(paste(z,Len))
  return(Len)
}


F_CSV <- function(Data) {
  write.csv(Data, paste("Test.csv", seq = ""), row.names = F)
}

F_DynamicMatching <- function(top = top, bottom = bottom, RefMZ = 200, MRP = 17500) { ## Dynamic matching
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


#' @export
F_DotProduct <- function(Q, L) {
  return(as.vector((Q %*% L) / (sqrt(sum(Q^2)) * sqrt(sum(L^2)))))
}
F_ExtractGNPSMS2 <- function(x) { # Extract MS2 according to index
  ## Extract all MS2 spectra that were associated with the candidate using the Index
  # print(x)
  Spectra <- GNPSSpectra[[x]]
  PeakList_matrix <- matrix(c(Spectra@mz, Spectra@intensity), ncol = 2, byrow = FALSE)
  DT <- data.table(PeakList_matrix)
  DT[, Index := x]
  # print(PeakList_matrix)
  return(DT)
}
F_getPEPFromScoreLambda <- function(targetScores, decoyScores, FigName) {
  library(reticulate)
  # source_python("inst/extdata/qvality_extracted.py")
  source_python("d:/github/dynamic/inst/extdata/qvality_extracted.py")
  Lambda <- getPEPFromScoreLambda(targetScores, decoyScores, FigName)
  # LambdaDynamic<-F_getPEPFromScoreLambda(TarDec$dpD, TarDec$dpD.decoy,'DynamicCurve')
}
joinPeaks <- function(x, y, type = "outer", tolerance = 0, ppm = 10, ...) {
  map <- MsCoreUtils::join(x[, 1], y[, 1], type = type, tolerance = tolerance,
                           ppm = ppm, ...)
  list(x = x[map$x, , drop = FALSE], y = y[map$y, , drop = FALSE])
}
MSsim <- function(alignment) { ## similarity score calculation
  alignment <- data.frame(alignment)
  # print(alignment)
  # score <- as.numeric(0)
  # if ((sum(alignment[, 2]) > 0) & (sum(alignment[, 3]) > 0)) {
  # alignment <- alignment[alignment[, 1] >= x.threshold, ]
  u <- alignment[, 2]
  v <- alignment[, 3]
  score <- as.vector((u %*% v) / (sqrt(sum(u^2)) * sqrt(sum(v^2))))
  # }
  return(score)
}
maxTic <- function(z) {
  z[[which.max(lapply(intensity(z), sum))]]
}
F_MergeLibrarySearch <- function(Dir) {
  FileNames <- dir(Dir, full.names = TRUE, recursive = TRUE)
  length(FileNames)
  Result.rbind <- data.table::data.table()
  for (id in FileNames) {
    # id <- FileNames[3916]
    FileNames[3911]
    ResultId <- data.table::data.table(read.csv(id))
    query_id <- sub(id, pattern = "_.*", replacement = "", perl = TRUE) %>%
      sub(id, pattern = ".*/", replacement = "", perl = TRUE)
    ResultId[, query_id := query_id]
    mztol <- sub(id, pattern = ".csv", replacement = "", perl = TRUE) %>%
      sub(id, pattern = ".*_", replacement = "", perl = TRUE)
    ResultId[, mztol := mztol]
    Result.rbind <- rbind(Result.rbind, ResultId, fill = TRUE)
  }
  return(Result.rbind)
}
.onAttach <- function(...) {
  packageStartupMessage("\nUse 'msGUI()' to start the GUI program.")
  # msGUI::msGUI()
}
F_peaks_compare_All <- function(x, y, MAPFUN = joinPeaks, FUN = ndotproduct,
                                tolerance = 0, ppm = 10,   ...) {
  ## Without precursor selection
  mat <- matrix(NA_real_, nrow = length(x), ncol = length(y), dimnames = list(names(x), names(y)))
  for (i in seq_along(x)) {
    for (j in seq_along(y)) {
      peak_map <- MAPFUN(x[[i]], y[[j]], tolerance = tolerance, ppm = ppm)
      mat[i, j] <- FUN(peak_map[[1L]], peak_map[[2L]])
    }
  }
  return(mat)
}
F_peaks_compare_Sel <- function(x, y, MAPFUN = joinPeaks, FUN = ndotproduct, UniInchikey,
                                tolerance = 0, ppm = 10, MRP = 70000, RefMZ = 200,  ...) {
  ## Dynamic tol for precursor selection
  mat <- matrix(NA_real_, nrow = length(x), ncol = length(y), dimnames = list(names(x), names(y)))
  A <- 1 / (MRP * (RefMZ^0.5))
  B <- A / 2.35482
  for (i in seq_along(x)) {
    for (j in seq_along(y)) {
      # i <- 1
      # j <- 2
      peak_map <- MAPFUN(x[[i]], y[[j]], tolerance = tolerance, ppm = ppm)
      # UniInchikey[inchikey_14_precursor_mz=="OUSYWCQYMPDAEO_188.0818",]
      InchiX <- names(x)[[i]]
      InchiY <- names(y)[[j]]
      Pmz_x <- UniInchikey[inchikey_14_precursor_mz==InchiX ,]$precursor_mz
      Pmz_y <- UniInchikey[inchikey_14_precursor_mz==InchiY,]$precursor_mz
      # mat[i, j] <- 0
      Tol_x <- B * Pmz_x ^1.5 + Pmz_x / 1000000
      Tol_y <- B * Pmz_y ^1.5 + Pmz_y / 1000000
      print(paste0(names(x)[[i]],'_Tol_x_',Tol_x))
      # print(paste0(names(y)[[j]],Pmz_y,'_Tol_y_',Tol_y))
      if (abs(Pmz_x - Pmz_y) < Tol_x | abs(Pmz_x - Pmz_y) < Tol_y) {
        mat[i, j] <- FUN(peak_map[[1L]], peak_map[[2L]])}
    }
  }
  return(mat)
}
F_PrecursorFilt <- function(DT) {
  print(paste0("before filtration ", nrow(DT)))
  DT[, MS1Tol := F_CalPMMT(library_precursor_mz, 70000, 200)]
  DT <- DT[MZDiff < MS1Tol]
  print(paste0("after filtration ", nrow(DT)))
  return(DT)
}
F_Selection<-function(df,greater,lessorequal=1){
  df[df<=greater]<-NA # remove the cells less than the lower limit
  df<-df[apply(df, 1, function(x) !all(is.na(x))),]
  df<-df[, colSums(is.na(df)) != nrow(df)]
  df[df>lessorequal]<-NA # remove the cells greater than the upper limit
  df<-df[apply(df, 1, function(x) !all(is.na(x))),]
  df<-df[, colSums(is.na(df)) != nrow(df)]
  df[is.na(df)] <- 0
  return(df)
}
F_TopCutoff <- function(DT, TopCutoff){
  TarDynamic <- setorder(DT[mztol=='NA'], -dpc) %>% .[, head(.SD, TopCutoff), by='query_qpid'] %>% .[, Database:='Target']
  Tar0.005 <- setorder(DT[mztol==0.005], -dpc) %>% .[, head(.SD, TopCutoff), by='query_qpid'] %>% .[, Database:='Target']
  Tar0.028 <- setorder(DT[mztol==0.028], -dpc) %>% .[, head(.SD, TopCutoff), by='query_qpid'] %>% .[, Database:='Target']
  Tar0.050 <- setorder(DT[mztol==0.050], -dpc) %>% .[, head(.SD, TopCutoff), by='query_qpid'] %>% .[, Database:='Target']
  
  DecDynamic <- setorder(DT[mztol=='NA'], -dpc.decoy) %>% .[, head(.SD, TopCutoff), by='query_qpid'] %>% .[, Database:='Decoy']
  Dec0.005 <- setorder(DT[mztol==0.005], -dpc.decoy) %>% .[, head(.SD, TopCutoff), by='query_qpid'] %>% .[, Database:='Decoy']
  Dec0.028 <- setorder(DT[mztol==0.028], -dpc.decoy) %>% .[, head(.SD, TopCutoff), by='query_qpid'] %>% .[, Database:='Decoy']
  Dec0.050 <- setorder(DT[mztol==0.050], -dpc.decoy) %>% .[, head(.SD, TopCutoff), by='query_qpid'] %>% .[, Database:='Decoy']
  Top <- rbind(TarDynamic, Tar0.005, Tar0.028, Tar0.050, DecDynamic, Dec0.005, Dec0.028, Dec0.050)
  setorder(Top, -dpc, query_qpid)
  # print(Top[query_id== 9762])
  return(Top)
}
