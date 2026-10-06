#### Xiaodong 2026.10 ####
# The computational time to process a given number of spectra is an important aspect
# of a new algorithm. In order to facilitate the usage of potentially interested users,
# we have tested the computational time of 1,500 MS/MS spectra as query
# and 100,000 spectra as library on a computer equuipped with AMD Ryzen 9 9950X 16-core CPU with 96 GB RAM.
# Here we need to load the SpectralMatching2.R function in the package
library(dplyr)
library(data.table)
library(foreach)
library(parallel)
source('d:/github/dynamic/R/SpectralMatching2.R')
Output <- 'D:/Onedrive/github/dynamic/output'
setwd(Output)
Path <- 'd:/OneDrive/github/dynamic/input/Sqlite/' # here the user needs to change into your own local path accordingly
#### Access the library information ####
l_dbPthValue <- paste0(Path,'GNPS_RepairMeta_FilPeaks5_Orb.sqlite') #
con <- DBI::dbConnect(RSQLite::SQLite(), l_dbPthValue)
library_spectra_meta <- con %>% dplyr::tbl("library_spectra_meta") %>% dplyr::collect() %>% as.data.table(.)
metab_compound <- con %>% dplyr::tbl("metab_compound") %>% dplyr::collect() %>% as.data.table(.)
library_spectra <- con %>% dplyr::tbl("library_spectra") %>% dplyr::collect() %>% as.data.table(.)
Meta <- library_spectra_meta
names(Meta) # check the names of the Meta
nrow(Meta) # 323203
set.seed(123)
#### Create the query sqlite database with 1500 Spectra ####
Meta1500 <- Meta[sample(.N, 1500)] # random selection of 1500 as query
con_Query1500 <- DBI::dbConnect(RSQLite::SQLite(), paste0(Path,'GNPS_Query1500.sqlite'))
DBI::dbWriteTable(con_Query1500, name = "library_spectra_meta", value = Meta1500,overwrite=TRUE)
names(library_spectra)
library_spectra_Query1500 <- library_spectra[library_spectra_meta_id %in% Meta1500$id, ]
unique(library_spectra_Query1500, by = "library_spectra_meta_id") # to double check it is 1500 spectra
DBI::dbWriteTable(con_Query1500, name = "library_spectra", value = library_spectra_Query1500,overwrite=TRUE)
DBI::dbWriteTable(con_Query1500, name = "metab_compound", value = metab_compound,overwrite=TRUE)
#### Create the library sqlite database with 100000 Spectra ####
Meta100000 <- Meta[sample(.N, 100000)]
con_Query2 <- DBI::dbConnect(RSQLite::SQLite(), paste0(Path,'GNPS_Library100000.sqlite'))
DBI::dbWriteTable(con_Query2, name = "library_spectra_meta", value = Meta100000,overwrite=TRUE)
names(library_spectra)
library_spectra_Query2 <- library_spectra[library_spectra_meta_id %in% Meta100000$id, ]
unique(library_spectra_Query2, by = "library_spectra_meta_id") # to double check it is 100000 spectra
DBI::dbWriteTable(con_Query2, name = "library_spectra", value = library_spectra_Query2,overwrite=TRUE)
DBI::dbWriteTable(con_Query2, name = "metab_compound", value = metab_compound,overwrite=TRUE)

#### Search query againtst library to test the computation time ####
Pth_Query1500 <- paste0(Path,"GNPS_Query1500.sqlite")
Pth_Library100000 <- paste0(Path,"GNPS_Library100000.sqlite")
## Precursor selection 1 core, Finished in 13.9 s (0.0 s per query spectrum)
Matched <- SpectralMatching2(q_dbPth = Pth_Query1500, l_dbPth = Pth_Library100000, cores = 1, usePrecursors = TRUE)
## Precursor selection 8 core, 162306 hits in 49.8 s (0.033 s per query spectrum)
Matched <- SpectralMatching2(q_dbPth = Pth_Query1500, l_dbPth = Pth_Library100000, cores = 8, usePrecursors = TRUE)
## No precursor selection 1 core, Finished in 2828.9 s (1.886 s per query spectrum)
Matched <- SpectralMatching2(q_dbPth = Pth_Query1500, l_dbPth = Pth_Library100000, cores = 1, usePrecursors = FALSE)
## No precursor selection 8 core, Finished in 497.5 s (0.332 s per query spectrum)
Matched <- SpectralMatching2(q_dbPth = Pth_Query1500, l_dbPth = Pth_Library100000, cores = 8, usePrecursors = FALSE)

