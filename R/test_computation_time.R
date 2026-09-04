#### Xiaodong 2023.02 ####
# The computational time to process a given number of spectra is an important aspect
# of a new algorithm. In order to facilitate the usage of potentially interested users,
# we have tested the computational time of 1,500 MS/MS spectra as query
# and 100,000 spectra as library on a laptop equipped with an Intel Core i7-8550U CPU,
# 1 TB HDD, and 16GB RAM.
# Here we need to load the SpectralMatching2.R function in the package
library(dplyr)
library(data.table)
library(foreach)
library(parallel)
#### Access the library information ####
Path <- 'd:/OneDrive/github/dynamic/input/Sqlite/' # here the user needs to change into your own local path accordingly
l_dbPthValue <- paste0(Path,'MoNA-export-All_LC-MS-MS_Orbitrap_MetaFil_0410.sqlite') #
con <- DBI::dbConnect(RSQLite::SQLite(), l_dbPthValue)
library_spectra_meta <- con %>%
  dplyr::tbl("library_spectra_meta") %>%
  dplyr::collect() %>%
  as.data.table(.)
metab_compound <- con %>%
  dplyr::tbl("metab_compound") %>%
  dplyr::collect() %>%
  as.data.table(.)
library_spectra <- con %>%
  dplyr::tbl("library_spectra") %>%
  dplyr::collect() %>%
  as.data.table(.)
Meta <- library_spectra_meta
names(Meta) # check the names of the Meta
nrow(Meta) # 57865
set.seed(123)
#### Create the query sqlite database with 1500 Spectra ####
Meta1500 <- Meta[sample(.N, 1500)] # random selection of 1500 as query
con_Query1500 <- DBI::dbConnect(RSQLite::SQLite(), paste0(Path,'Query1500.sqlite'))
DBI::dbWriteTable(con_Query1500, name = "library_spectra_meta", value = Meta1500)
names(library_spectra)
library_spectra_Query1500 <- library_spectra[library_spectra_meta_id %in% Meta1500$id, ]
unique(library_spectra_Query1500, by = "library_spectra_meta_id") # to double check it is 1500 spectra
DBI::dbWriteTable(con_Query1500, name = "library_spectra", value = library_spectra_Query1500)
DBI::dbWriteTable(con_Query1500, name = "metab_compound", value = metab_compound)
#### Create the query sqlite database with 2 Spectra ####
Meta2 <- Meta1500[1:2]
con_Query2 <- DBI::dbConnect(RSQLite::SQLite(), paste0(Path,'Query2.sqlite'))
DBI::dbWriteTable(con_Query2, name = "library_spectra_meta", value = Meta2,overwrite=TRUE)
names(library_spectra)
library_spectra_Query2 <- library_spectra[library_spectra_meta_id %in% Meta2$id, ]
unique(library_spectra_Query2, by = "library_spectra_meta_id") # to double check it is 1500 spectra
DBI::dbWriteTable(con_Query2, name = "library_spectra", value = library_spectra_Query2,overwrite=TRUE)
DBI::dbWriteTable(con_Query2, name = "metab_compound", value = metab_compound,overwrite=TRUE)

#### Search query againtst library to test the computation time ####
Pth_Query2 <- paste0(Path,"Query2.sqlite")
Pth_Query1500 <- paste0(Path,"Query1500.sqlite")
Pth_Library <- l_dbPthValue
## Searching 57865 library with 2 query, 1 core, 128 hits in 6.3 s (3.150 s per query spectrum)
Matched <- SpectralMatching2(q_dbPth = Pth_Query2, l_dbPth = Pth_Library, cores = 1)
## Searching 57865 library with 2 query, 4 core, 128 hits in 14.8 s (7.375 s per query spectrum)
Matched <- SpectralMatching2(q_dbPth = Pth_Query2, l_dbPth = Pth_Library, cores = 4)
## Searching 57865 library with 1500 query, 1 core, 107021 hits in 43.7 s (0.029 s per query spectrum)
Matched <- SpectralMatching2(q_dbPth = Pth_Query1500, l_dbPth = Pth_Library, cores = 1)
## Searching 57865 library with 1500 query, 4 core, 107021 hits in 27.0 s (0.018 s per query spectrum)
Matched <- SpectralMatching2(q_dbPth = Pth_Query1500, l_dbPth = Pth_Library, cores = 4)


