#### Xiaodong 2023.02 ####
# The computational time to process a given number of spectra is an important aspect 
# of a new algorithm. In order to facilitate the usage of potentially interested users, 
# we have tested the computational time of 1,500 MS/MS spectra as query 
# and 100,000 spectra as library on a laptop equipped with an Intel Core i7-8550U CPU, 
# 1 TB HDD, and 16GB RAM. The reported time using the tictoc R package is 5.6 hours 
# for annotation of 1,500 spectra, resulting in around 14 seconds per spectra on average.
# The library used in this script can be found in https://doi.org/10.5281/zenodo.7549795
# The 20210423_mona.sqlite file is originally downloaded via https://github.com/computational-metabolomics/msp2db/releases/tag/v0.0.14-mona-23042021
library(dplyr)
library(data.table)
#### Access the library information ####
Path <- 'd:/github/dynamic/input/sqlite/' # here the user needs to change into your own local path accordingly
l_dbPthValue <- paste0(Path,'20210423_mona.sqlite') # 
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
nrow(Meta) # 661421
set.seed(123)
#### Create the query sqlite database with 2 Spectra ####
Meta2 <- Meta1500[1:2]
con_Query2 <- DBI::dbConnect(RSQLite::SQLite(), paste0(Path,'Query2.sqlite'))
DBI::dbWriteTable(con_Query2, name = "library_spectra_meta", value = Meta2,overwrite=TRUE)
names(library_spectra)
library_spectra_Query2 <- library_spectra[library_spectra_meta_id %in% Meta2$id, ]
unique(library_spectra_Query2, by = "library_spectra_meta_id") # to double check it is 1500 spectra
DBI::dbWriteTable(con_Query2, name = "library_spectra", value = library_spectra_Query2,overwrite=TRUE)
DBI::dbWriteTable(con_Query2, name = "metab_compound", value = metab_compound,overwrite=TRUE)
#### Create the query sqlite database with 1500 Spectra ####
Meta1500 <- Meta[sample(.N, 1500)] # random selection of 1500 as query
con_Query1500 <- DBI::dbConnect(RSQLite::SQLite(), paste0(Path,'Query1500.sqlite'))
DBI::dbWriteTable(con_Query1500, name = "library_spectra_meta", value = Meta1500)
names(library_spectra)
library_spectra_Query1500 <- library_spectra[library_spectra_meta_id %in% Meta1500$id, ]
unique(library_spectra_Query1500, by = "library_spectra_meta_id") # to double check it is 1500 spectra
DBI::dbWriteTable(con_Query1500, name = "library_spectra", value = library_spectra_Query1500)
DBI::dbWriteTable(con_Query1500, name = "metab_compound", value = metab_compound)
#### Create the library sqlite database with 100000 Spectra ####
Meta10e5 <- Meta[sample(.N, 100000)] # random selection
con_Library10e5 <- DBI::dbConnect(RSQLite::SQLite(), paste0(Path, "Library10e5.sqlite"))
DBI::dbWriteTable(con_Library10e5, name = "library_spectra_meta", value = Meta10e5)
library_spectra_Library10e5 <- library_spectra[library_spectra_meta_id %in% Meta10e5$id, ]
unique(library_spectra_Library10e5, by = "library_spectra_meta_id") # to double check it is 10e5 spectra
DBI::dbWriteTable(con_Library10e5, name = "library_spectra", value = library_spectra_Library10e5)
DBI::dbWriteTable(con_Library10e5, name = "metab_compound", value = metab_compound)
#### Create the query sqlite database with 3 Spectra using priviously created query####
con <- DBI::dbConnect(RSQLite::SQLite(), 'd:/github/dynamic/input/msp/Query.sqlite')
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
nrow(Meta) # 661421
setorder(Meta,id)
Meta3 <- Meta[1:3]
con_Query3 <- DBI::dbConnect(RSQLite::SQLite(), paste0(Path,'Query3.sqlite'))
DBI::dbWriteTable(con_Query3, name = "library_spectra_meta", value = Meta3,overwrite=TRUE)
library_spectra_Query3 <- library_spectra[library_spectra_meta_id %in% Meta3$id, ]
unique(library_spectra_Query3, by = "library_spectra_meta_id") # to double check it is 1500 spectra
DBI::dbWriteTable(con_Query3, name = "library_spectra", value = library_spectra_Query3,overwrite=TRUE)
DBI::dbWriteTable(con_Query3, name = "metab_compound", value = metab_compound,overwrite=TRUE)

#### Search query againtst library to test the computation time ####
MRPValue <<- 17500
RefMZValue <<- 200
Pth_Query2 <- paste0(Path,"Query2.sqlite")
Pth_Query3 <- paste0(Path,"Query3.sqlite")
Pth_Query10 <- paste0(Path,"Query10.sqlite")
Pth_Query1500 <- paste0(Path,"Query1500.sqlite")
Pth_Library10e5 <- paste0(Path,"Library10e5.sqlite")
setwd("output")
# for test begin
library(dplyr)
library(data.table)
library(parallel)
library(foreach)
q_pol = l_pol = q_instrumentTypes = l_instrumentTypes = q_instruments = l_instruments = q_sources = l_sources = q_pids = l_pids =q_accessions = l_accessions =rttol= mztol=NA;
q_rtrange = l_rtrange = c(NA, NA)
q_dbType = l_dbType = "sqlite"
usePrecursors = TRUE
updateDb = copyDb = decoy = FALSE 
outPth = "sm_result.sqlite"
q_ppmPrec = l_ppmPrec = 5;
raW = 1; mzW = 0; MRPValue = 17500; RefMZValue = 200; HalfLenth= 75
libsch = 'reverse'

Matched <- SpectralMatching(q_dbPth = Pth_Query3, l_dbPth = Pth_Sensus,cores = 3, l_ppmPrec = 5,decoy = FALSE)
head(Matched)
# for test end
# Here we used the for loop to perform the library search per query, some of the query has the related
# candidates in the database, some of the query do not have, thus returned no results
# The total time is 20272.5 seconds, thus 14 seconds per spectra on average.


