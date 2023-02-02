# This file contains functions used for the spectra processing
# The functions are sorted alphabetically
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

#' @title F_CountCutoff
#' @description
#' Count how many entries are above the specified cutoffs
#' @export
F_CountCutoff <- function(DT) {
  setorder(DT, -dpc)
  CountCutoff.rbind <- data.table()
  for (id in seq(0.01, 1, 0.01)) {
    CountCutoff <- data.table("Cutoff" = id, "Count" = nrow(DT[dpc >= id]))
    CountCutoff.rbind <- rbind(CountCutoff.rbind, CountCutoff)
  }
  return(CountCutoff.rbind)
}

#' @title F_CountSpectraPerInchikey
#' @description
#' Count the number of spectra for each inchikey (first 14 charactors)
#' @export
F_CountSpectraPerInchikey <- function(Inchikey) {
  print(paste0(Inchikey, nrow(Meta[inchikey_14_precursor_type == Inchikey])))
  return(nrow(Meta[inchikey_14_precursor_type == Inchikey]))
}

#' @title F_CountPeaksPerSpectra
#' @description
#' Count the number of peaks (spectra length) for each Spectra
#' @export
F_CountPeaksPerSpectra<-function(z){
  Len<-nrow(library_spectra[library_spectra_meta_id==z])
  print(paste(z,Len))
  return(Len)
  # Test<-TarDecUni[1:2]
  # Test[,SpectraLength:=sapply(Test$library_lpid, F_CountSpectraLength)]
}

#' @title F_CountIdentical
#' @description
#' Count how many pairs of query and library MS/MS are identical
#' @export
F_CountIdentical <- function(DT) {
  DT <- unique(DT, by = c("query_inchikey14", "library_inchikey14", "library_precursor_mz"))
  Uni <- max(nrow(unique(DT, by = "query_accession")), nrow(unique(DT, by = "library_accession")))
  print(paste("Identical pairs", nrow(DT) / 2, "Spectra number", Uni / 2))
}

#' @title F_CountCandidates
#' @description
#' Count how many candidates per query
#' @export
F_CountCandidates <- function(z) {
  Len <- nrow(TarDec[query_qpid == z])
  print(paste(z, Len))
  return(Len)
}

#' @title F_CSV
#' @description
#' Easy way to export the .csv files, used for checking results quickly.
#' @export
F_CSV <- function(Data) {
  write.csv(Data, paste("Test.csv", seq = ""), row.names = F)
}

#' @title F_DotProduct
#' @description
#' Calculate the dot product cosine using two vector
#' @export
F_DotProduct <- function(Q, L) {
  return(as.vector((Q %*% L) / (sqrt(sum(Q^2)) * sqrt(sum(L^2)))))
}

#' @title F_ExtractGNPSMS2
#' @description
#' Extract all MS2 spectra that were associated with the candidate using the Index
#' @export
F_ExtractGNPSMS2 <- function(x) { 
  # print(x)
  Spectra <- GNPSSpectra[[x]]
  PeakList_matrix <- matrix(c(Spectra@mz, Spectra@intensity), ncol = 2, byrow = FALSE)
  DT <- data.table(PeakList_matrix)
  DT[, Index := x]
  # print(PeakList_matrix)
  return(DT)
}

#' @title F_getPEPFromScoreLambda
#' @import reticulate
#' @description
#' This function is modified based on the qvality.py file in [triqler] (https://github.com/statisticalbiotechnology/triqler/blob/master/triqler/qvality.py)
#' Get PEP score according to the target scores and decoy scores, a regression figure is exported as well
#' @export
F_getPEPFromScoreLambda <- function(targetScores, decoyScores, FigName) {
  # library(reticulate)
  # source_python("inst/extdata/qvality_extracted.py")
  reticulate::source_python("inst/extdata/qvality_extracted.py")
  Lambda <- getPEPFromScoreLambda(targetScores, decoyScores, FigName)
}

#' @title F_MergeLibrarySearch
#' @description
#' Merge the .csv files after library search into a single data table
#' @import data.table
#' @export
F_MergeLibrarySearch <- function(Dir){
  FileNames<-dir(Dir, full.names = TRUE, recursive = TRUE)
  print(paste('File length',length(FileNames)))
  Result.rbind<-data.table::data.table()
  for (id in FileNames) {
    # id <- FileNames[3916]
    ResultId <- data.table::data.table(read.csv(id))
    id.short <- sub(".*/","",id) %>% sub(".csv","",.)
    id.split <-  unlist(strsplit(id.short, '_'))
    ResultId[, query_id := id.split[1]]
    ResultId[, mztol := id.split[2]]
    ResultId[, libsch := id.split[3]]
    Result.rbind<-rbind(Result.rbind,ResultId,fill=TRUE)
  }
  return(Result.rbind)
}

.onAttach <- function(...) {
  packageStartupMessage("\nUse dynamic mass tolerance for peak matching")
}

#' @title F_PrecursorFilt
#' @description
#' Filter the results based on the precursor m/z difference
#' @export
F_PrecursorFilt <- function(DT) {
  print(paste0("before filtration ", nrow(DT)))
  DT[, MS1Tol := F_CalPMMT(library_precursor_mz, 70000, 200)]
  DT <- DT[MZDiff < MS1Tol]
  print(paste0("after filtration ", nrow(DT)))
  return(DT)
}

#' @title F_Selection
#' @description
#' Selection of the heat map cells based on the similarity scores
#' @export
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

