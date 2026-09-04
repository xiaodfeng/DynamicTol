#### Create the Xcorr decoy database based on Query of Mona ####
# The scripts are used for creating the decoy database based on Mona. To do this,
# Transfer the .sqlite database into .ms database suitable for sirus https://boecker-lab.github.io/docs.sirius.github.io/io/#input
# Analyze the .ms database
# Collect the decoys from the sirius and write back into the .sqlite database.
generate_decoy <- function(id, top_tmp, PMMTmax, PMMTmin) {
  # Define a function to calculate decoy mz values for a given id
  decoy <- data.table(mz = top_tmp$mz + 1.5 * (id / abs(id)) * PMMTmax + id * PMMTmin,
                      i = top_tmp$i,
                      library_spectra_meta_id = top_tmp$library_spectra_meta_id)
  return(decoy)
}
F_MakeXcorrDecoyDB <- function(inPth, outPth, noise_cutoff=0){
  ## Access the .sqlite database
  con_Query <- DBI::dbConnect(RSQLite::SQLite(),inPth)
  library_spectra_meta_Query <- con_Query %>% dplyr::tbl("library_spectra_meta") %>% dplyr::collect() %>% as.data.table(.)
  metab_compound_Query <- con_Query %>% dplyr::tbl("metab_compound") %>% dplyr::collect() %>% as.data.table(.)
  library_spectra_source_Query <- con_Query %>% dplyr::tbl("library_spectra_source") %>% dplyr::collect() %>% as.data.table(.)
  library_spectra_Query <- con_Query %>% dplyr::tbl("library_spectra") %>% dplyr::collect() %>% as.data.table(.)
  DBI::dbDisconnect(con_Query)
  ## Create the decoy
  MetaIds <- unique(library_spectra_Query$library_spectra_meta_id)
  if (noise_cutoff>0) {
    xcorr.rbind <- foreach(keyword = MetaIds) %do% {
    # keyword <- 47962
    Spec <- library_spectra_Query[library_spectra_meta_id == keyword]
    Spec$i <- 100 * Spec$i / max(Spec$i)
    Spec <- Spec[i > noise_cutoff]
    PMMTmin<-F_CalPMMT(min(Spec$mz)) #
    PMMTmax<-F_CalPMMT(max(Spec$mz)) #
    ids <- c(seq(-75, -1), seq(1, 75))
    decoy_list <- lapply(ids, generate_decoy, top_tmp = Spec, PMMTmax = PMMTmax, PMMTmin = PMMTmin)
    decoy.rbind <- do.call(rbind, decoy_list)
    decoy.rbind <- decoy.rbind [mz>0] # %>% setorder(.,mz) # exclude the entries with mz less than zero
    decoy.rbind <- unique(decoy.rbind,by=c('mz','i'))  # unique
    return(decoy.rbind)
  }
    }else{
    ## Without noise filtration
    xcorr.rbind <- foreach(keyword = MetaIds) %do% {
        # keyword <- 57628
        Spec <- library_spectra_Query[library_spectra_meta_id == keyword]
        PMMTmin<-F_CalPMMT(min(Spec$mz)) #
        PMMTmax<-F_CalPMMT(max(Spec$mz)) #
        ids <- c(seq(-75, -1), seq(1, 75))
        decoy_list <- lapply(ids, generate_decoy, top_tmp = Spec, PMMTmax = PMMTmax, PMMTmin = PMMTmin)
        decoy.rbind <- do.call(rbind, decoy_list)
        decoy.rbind <- decoy.rbind [mz>0] # %>% setorder(.,mz) # exclude the entries with mz less than zero
        decoy.rbind <- unique(decoy.rbind,by=c('mz','i'))  # unique
        return(decoy.rbind)
      }
  }
  xcorr.rbind <- bind_rows(xcorr.rbind) # Combine the results from lists
  ## Write out the database
  Pth_QueryDecoy <- outPth
  con_QueryDecoy <- DBI::dbConnect(RSQLite::SQLite(),Pth_QueryDecoy)
  DBI::dbWriteTable(con_QueryDecoy, name='library_spectra_meta', value=library_spectra_meta_Query,overwrite=TRUE)
  DBI::dbWriteTable(con_QueryDecoy, name='metab_compound',value=metab_compound_Query,overwrite=TRUE)
  DBI::dbWriteTable(con_QueryDecoy, name='library_spectra_source', value=library_spectra_source_Query,overwrite=TRUE)
  DBI::dbWriteTable(con_QueryDecoy, name='library_spectra', value=xcorr.rbind,overwrite=TRUE)
  DBI::dbDisconnect(con_QueryDecoy)
}

F_MakeXcorrDecoyDB(inPth='d:/OneDrive/github/dynamic/input/Sqlite/Query500Pos500Neg_MonaOrb_1117.sqlite',
                   outPth='d:/OneDrive/github/dynamic/input/Sqlite/Query500Pos500Neg_MonaOrb_1117_XcorrCutoff.sqlite',
                   , noise_cutoff=1)
F_MakeXcorrDecoyDB(inPth='d:/OneDrive/github/dynamic/input/Sqlite/Query500Pos500Neg_MonaOrb_1117.sqlite',
                   outPth='d:/OneDrive/github/dynamic/input/Sqlite/Query500Pos500Neg_MonaOrb_1117_Xcorr.sqlite',
                   , noise_cutoff=0)
## Perform the library search using the created  decoy
## For query decoy
Pth_Query <- "d:/OneDrive/github/dynamic/input/Sqlite/Query500Pos500Neg_MonaOrb_1117_XcorrCutoff.sqlite"
Pth_Library <- "d:/OneDrive/github/dynamic/input/Sqlite/LibraryNoNeg_MonaOrb_1117_Sensus.sqlite"
Dir <- "d:/OneDrive/github/dynamic/output/LibrarySearch/MonaOrb/20260826/QueryXcorrCutoff.VS.Sensus"
dir.create(Dir, recursive = TRUE)
setwd(Dir)
Test <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, mztol = NA, cores = 2,q_pids = 115,decoy = FALSE)
Test
XcorrNA <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, mztol = NA,cores = 4,decoy = FALSE)
#### creating the decoy database based on Mona using sirius ####
# The scripts are used for creating the decoy database based on Mona using sirius. To do this,
# Transfer the .sqlite database into .ms database suitable for sirus https://boecker-lab.github.io/docs.sirius.github.io/io/#input
# Analyze the .ms database using sirius
# Collect the decoys from the sirius and write back into the .sqlite database.

#### Transfer the .sqlite database into .ms database suitable for sirus
## Access the .sqlite database
con_Sensus <- DBI::dbConnect(RSQLite::SQLite(),'d:/OneDrive/github/dynamic/input/Sqlite/LibraryNoNeg_MonaOrb_1117_Sensus.sqlite')
library_spectra_meta_Sensus <- con_Sensus %>% dplyr::tbl("library_spectra_meta") %>% dplyr::collect() %>% as.data.table(.)
metab_compound_Sensus <- con_Sensus %>% dplyr::tbl("metab_compound") %>% dplyr::collect() %>% as.data.table(.)
library_spectra_source_Sensus <- con_Sensus %>% dplyr::tbl("library_spectra_source") %>% dplyr::collect() %>% as.data.table(.)
library_spectra_Sensus <- con_Sensus %>% dplyr::tbl("library_spectra") %>% dplyr::collect() %>% as.data.table(.)
DBI::dbDisconnect(con_Sensus)

con_Query <- DBI::dbConnect(RSQLite::SQLite(),'d:/OneDrive/github/dynamic/input/Sqlite/Query500Pos500Neg_MonaOrb_1117.sqlite')
library_spectra_meta_Query <- con_Query %>% dplyr::tbl("library_spectra_meta") %>% dplyr::collect() %>% as.data.table(.)
metab_compound_Query <- con_Query %>% dplyr::tbl("metab_compound") %>% dplyr::collect() %>% as.data.table(.)
library_spectra_source_Query <- con_Query %>% dplyr::tbl("library_spectra_source") %>% dplyr::collect() %>% as.data.table(.)
library_spectra_Query <- con_Query %>% dplyr::tbl("library_spectra") %>% dplyr::collect() %>% as.data.table(.)
DBI::dbDisconnect(con_Query)




## Create the decoy sirius based on library_spectra_Sensus
InputMs <- "D:/github/sirius/input/LibraryNoNeg_MonaOrb_1117_Sensus.ms" # Open the MS file for writing
DT <- library_spectra_meta_Sensus
names(DT)
metab_compound_Sensus # 4747
metab_compound_Sensus_Uni <- unique(metab_compound_Sensus,by="inchikey_id") # 4247
DT <- left_join(DT,metab_compound_Sensus_Uni[,c("inchikey_id","molecular_formula")],by="inchikey_id")
msFile <- file(InputMs, "w")
for (Ind in DT$id) {
  # Ind <- 8789
  print(DT[id==Ind,]$name)
  spectrum <- list(
    compound = DT[id==Ind,]$name,
    formula = DT[id==Ind,]$molecular_formula,
    parentmass = DT[id==Ind,]$precursor_mz,
    ionization = DT[id==Ind,]$precursor_type,
    mz2Values = library_spectra_Sensus[library_spectra_meta_id==Ind,]$mz,
    intensity2Values = library_spectra_Sensus[library_spectra_meta_id==Ind,]$i
  )
  writeMSSpectrum(msFile,spectrum)
}
close(msFile)

## Create the decoy sirius based on library_spectra_Query
InputMs <- "D:/github/sirius/input/Query500Pos500Neg_MonaOrb_1117.ms" # Open the MS file for writing
DT <- library_spectra_meta_Query
names(DT)
metab_compound_Query # 1000
metab_compound_Query_Uni <- unique(metab_compound_Query,by="inchikey_id") # 1000
DT <- left_join(DT,metab_compound_Query_Uni[,c("inchikey_id","molecular_formula")],by="inchikey_id")
msFile <- file(InputMs, "w")
for (Ind in DT$id) {
  # Ind <- 57540
  # names(DT)
  print(DT[id==Ind,]$name)
  spectrum <- list(
    compound = DT[id==Ind,]$name,
    formula = DT[id==Ind,]$molecular_formula,
    parentmass = DT[id==Ind,]$precursor_mz,
    ionization = DT[id==Ind,]$precursor_type,
    mz2Values = library_spectra_Query[library_spectra_meta_id==Ind,]$mz,
    intensity2Values = library_spectra_Query[library_spectra_meta_id==Ind,]$i
  )
  writeMSSpectrum(msFile,spectrum)
}
close(msFile)

#### Analyze the .ms database using sirius
## For library
OutDir <- 'D:/github/sirius/output/LibraryNoNeg_MonaOrb_1117_Sensus'
dir.create(OutDir) # need to create the folder first, otherwise will have no folders
sirius_cmd <- paste('"D:/github/sirius/sirius.exe"',
                    '-i',"D:/github/sirius/input/LibraryNoNeg_MonaOrb_1117_Sensus.ms",
                    '-o', OutDir,
                    # '--processors 5',
                    'formula -p orbitrap passatutto write-summaries')
system(sirius_cmd)
## For query
OutDir <- 'D:/github/sirius/output/Query500Pos500Neg_MonaOrb_1117'
dir.create(OutDir) # need to create the folder first, otherwise will have no folders
sirius_cmd <- paste('"D:/github/sirius/sirius.exe"',
                    '-i','D:/github/sirius/input/Query500Pos500Neg_MonaOrb_1117.ms',
                    '-o', OutDir,
                    # '--processors 5',
                    'formula -p orbitrap passatutto write-summaries')
system(sirius_cmd)
#### Collect the decoys from the sirius to replace the library_spectra_Sensus
## For library
formula_identifications <- fread('D:/github/sirius/output/LibraryNoNeg_MonaOrb_1117_Sensus/formula_identifications.tsv')
decoys.rbind = data.table()
for (ind in formula_identifications$id) {
  # ind <- "5490_LibraryNoNeg_MonaOrb_1117_Sensus_Khayanthone"
  print(ind)
  NRow <- sub('_LibraryNoNeg.*','',ind) %>% as.numeric(.)
  DecoyDir <- dir(paste0('D:/github/sirius/output/LibraryNoNeg_MonaOrb_1117_Sensus/', ind,'/decoys/'), full.names = TRUE, recursive = TRUE)
  if (length(DecoyDir)>0) {
    decoy <- fread(DecoyDir) %>% setnames(.,'rel.intensity','i') %>% .[,c('mz','i')] %>% # extract the MS2
      .[,library_spectra_meta_id:= DT[NRow,]$id] #%>%   Add the meta_id index
    # .[,inchikey_14_precursor_mz:= DT[NRow,]$inchikey_14_precursor_mz] #  Add the meta_id index and inchikey_14_precursor_mz
    decoys.rbind <- rbind(decoys.rbind, decoy)
  }
}
## For query
formula_identifications <- fread('D:/github/sirius/output/Query500Pos500Neg_MonaOrb_1117/formula_identifications.tsv')
decoys.rbind = data.table()
for (ind in formula_identifications$id) {
  # ind <- '6_Query500Pos500Neg_MonaOrb_1117_methyl5Z-5-ethylidene-4-2-2R3S4S5R6R'
  print(ind)
  NRow <- sub('_Query.*','',ind) %>% as.numeric(.)
  DecoyDir <- dir(paste0('D:/github/sirius/output/Query500Pos500Neg_MonaOrb_1117/', ind,'/decoys/'), full.names = TRUE, recursive = TRUE)
  if (length(DecoyDir)>0) {
    decoy <- fread(DecoyDir) %>% setnames(.,'rel.intensity','i') %>% .[,c('mz','i')] %>% # extract the MS2
      .[,library_spectra_meta_id:= DT[NRow,]$id] # %>% #  Add the meta_id index
    # .[,inchikey_14_precursor_mz:= DT[NRow,]$inchikey_14_precursor_mz] #  Add the meta_id index and inchikey_14_precursor_mz
    decoys.rbind <- rbind(decoys.rbind, decoy)
  }
}

#### write back into the .sqlite database
## For library
Pth_SensusDecoy <- 'd:/OneDrive/github/dynamic/input/Sqlite/LibraryNoNeg_MonaOrb_1117_Sensus_Decoy.sqlite'
con_SensusDecoy <- DBI::dbConnect(RSQLite::SQLite(),Pth_SensusDecoy)
DBI::dbWriteTable(con_SensusDecoy, name='library_spectra_meta', value=library_spectra_meta_Sensus)
DBI::dbWriteTable(con_SensusDecoy, name='metab_compound',value=metab_compound_Sensus)
DBI::dbWriteTable(con_SensusDecoy, name='library_spectra_source', value=library_spectra_source_Sensus)
DBI::dbWriteTable(con_SensusDecoy, name='library_spectra', value=decoys.rbind)
DBI::dbDisconnect(con_SensusDecoy)
## For query
Pth_QueryDecoy <- 'd:/OneDrive/github/dynamic/input/Sqlite/Query500Pos500Neg_MonaOrb_1117_Decoy.sqlite'
con_QueryDecoy <- DBI::dbConnect(RSQLite::SQLite(),Pth_QueryDecoy)
DBI::dbWriteTable(con_QueryDecoy, name='library_spectra_meta', value=library_spectra_meta_Query)
DBI::dbWriteTable(con_QueryDecoy, name='metab_compound',value=metab_compound_Query)
DBI::dbWriteTable(con_QueryDecoy, name='library_spectra_source', value=library_spectra_source_Query)
DBI::dbWriteTable(con_QueryDecoy, name='library_spectra', value=decoys.rbind)
DBI::dbDisconnect(con_QueryDecoy)

#### Perform the library search using the created sensus decoy
## For library decoy
Pth_Query <- "d:/OneDrive/github/dynamic/input/Sqlite/Query500Pos500Neg_MonaOrb_1117.sqlite"
Pth_Library <- "d:/OneDrive/github/dynamic/input/Sqlite/LibraryNoNeg_MonaOrb_1117_Sensus_Decoy.sqlite"
Dir <- "d:/OneDrive/github/dynamic/output/LibrarySearch/MonaOrb/20260826/Query.VS.SiriusDecoy"
dir.create(Dir, recursive = TRUE)
setwd(Dir)
Test <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library,
                         mztol = NA, cores = 1,q_pids = 115,decoy = FALSE)
setorder(Test,-dpc)
Test
SiriusNA <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library,cores = 5)
## For query decoy
Pth_Query <- "d:/OneDrive/github/dynamic/input/Sqlite/Query500Pos500Neg_MonaOrb_1117_Decoy.sqlite"
Pth_Library <- "d:/OneDrive/github/dynamic/input/Sqlite/LibraryNoNeg_MonaOrb_1117_Sensus.sqlite"
Dir <- "d:/OneDrive/github/dynamic/output/LibrarySearch/MonaOrb/20260826/QueryDecoy.VS.Sensus"
dir.create(Dir, recursive = TRUE)
setwd(Dir)
SiriusNA <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, mztol = NA,cores = 4,decoy = FALSE)
#### creating the query and decoy database based on GNPS using sirius ####
# The scripts are used for creating the query and decoy database based on GNPS using sirius. To do this,
# In pycharm, transfer the OrbRemovePos.msp into .sqlite format using msp2db python package.
# In Rstudio, transfer the .sqlite database into .ms database suitable for sirus https://boecker-lab.github.io/docs.sirius.github.io/io/#input
# Analyze the .ms database using sirius
# Collect the decoys from the sirius and write back into the .sqlite database.
#### Prepare the query .sqlite database
# Dowload the .mgf file from https://massive.ucsd.edu/ProteoSAFe/dataset_files.jsp?task=b753bf1e39cb4875bdf3b786e747bc15#%7B%22table_sort_history%22%3A%22main.collection_asc%22%2C%22main.attachment_input%22%3A%22updates%2F2022-09-14_pmallard_c2fb482f%22%7D
Output <- 'e:/dynamic/DynamicTol/output'
MRPValue <- 17500
RefMZValue <- 200
#### QueryAnnoOrb0.002.VS.OrbRemovePosPlusAlignedInstrKnown0.2_86662
Pth_Query <- 'e:/dynamic/DynamicTol/input/Sqlite/QueryAnnoOrb0.002.sqlite' # as query
Pth_DecoySirius <- 'e:/dynamic/DynamicTol/input/Sqlite/QueryAnnoOrb0.002Sirius.sqlite' # as query
Pth_Library <- 'e:/dynamic/DynamicTol/input/Sqlite/OrbRemovePosPlusAlignedInstrKnown0.2.sqlite' #
setwd('e:/dynamic/DynamicTol/output/LibrarySearch/Plant/QueryAnnoOrb0.002.VS.OrbRemovePosPlusAlignedInstrKnown0.2_86662/Target')
# Test <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 2,mztol=NA,
#                          q_pids=100525, l_pids = 187578,decoy = TRUE) # for test
# Test
PlantQuery5 <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 5,decoy = TRUE)
PlantQuery10 <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 10,decoy = TRUE)
PlantQuery0.005 <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 0.005,decoy = TRUE)
PlantQuery0.028 <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 0.028,decoy = TRUE)
PlantQuery0.050 <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 0.050,decoy = TRUE)
PlantQueryNA <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7,decoy = TRUE)
setwd('e:/dynamic/DynamicTol/output/LibrarySearch/Plant/QueryAnnoOrb0.002.VS.OrbRemovePosPlusAlignedInstrKnown0.2_86662/Decoy')
PlantDecoySirius <- SpectralMatching(q_dbPth=Pth_DecoySirius, l_dbPth=Pth_Library, cores = 7,decoy = FALSE)
#### QueryAnnoOrb0.002.VS.Sensus0.2_13101
Pth_Query <- 'e:/dynamic/DynamicTol/input/Sqlite/QueryAnnoOrb0.002.sqlite' # as query
Pth_DecoySirius <- 'e:/dynamic/DynamicTol/input/Sqlite/QueryAnnoOrb0.002Sirius.sqlite' # as query
Pth_Library <- 'e:/dynamic/DynamicTol/input/Sqlite/GnpsSensus0.2.sqlite' #
setwd('e:/dynamic/DynamicTol/output/LibrarySearch/Plant/QueryAnnoOrb0.002.VS.Sensus0.2_13101/Target')
PlantQueryNA <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 1, q_pids=100525) # for test

setwd('e:/dynamic/DynamicTol/output/LibrarySearch/Plant/QueryAnnoOrb0.002.VS.Sensus0.2_13101/Decoy')
PlantDecoySirius <- SpectralMatching(q_dbPth=Pth_DecoySirius, l_dbPth=Pth_Library, cores = 7,decoy = FALSE)


#### QueryAnnoOrb2314.VS.OrbRemovePosPlusAlignedInstr185865
Pth_Query <- 'e:/dynamic/DynamicTol/input/Sqlite/QueryAnnoOrb.sqlite' # as query
Pth_DecoySirius <- 'e:/dynamic/DynamicTol/input/Sqlite/QueryAnnoOrbSirius.sqlite' # as query
Pth_Library <- 'e:/dynamic/DynamicTol/input/Sqlite/OrbRemovePosPlusAlignedInstr.sqlite' #
setwd('e:/dynamic/DynamicTol/output/LibrarySearch/Plant/QueryAnnoOrb2314.VS.OrbRemovePosPlusAlignedInstr185865/Target')
# PlantQueryNA <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 1, q_pids=1) # for test
PlantQuery5 <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 5,decoy = TRUE)
PlantQuery10 <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 10,decoy = TRUE)
PlantQuery0.005 <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 0.005,decoy = TRUE)
PlantQuery0.028 <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 0.028,decoy = TRUE)
PlantQuery0.050 <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 0.050,decoy = TRUE)
PlantQueryNA <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7,decoy = TRUE)
setwd('e:/dynamic/DynamicTol/output/LibrarySearch/Plant/QueryAnnoOrb2314.VS.OrbRemovePosPlusAlignedInstr185865/Decoy')
PlantDecoySirius <- SpectralMatching(q_dbPth=Pth_DecoySirius, l_dbPth=Pth_Library, cores = 7,decoy = FALSE)

#### QueryAnnoOrb2314.VS.Sensus22931
Pth_Query <- 'e:/dynamic/DynamicTol/input/Sqlite/QueryAnnoOrb.sqlite' # as query
Pth_DecoySirius <- 'e:/dynamic/DynamicTol/input/Sqlite/QueryAnnoOrbSirius.sqlite' # as query
Pth_Library <- 'e:/dynamic/DynamicTol/input/Sqlite/GnpsSensus.sqlite' #
setwd('e:/dynamic/DynamicTol/output/LibrarySearch/Plant/QueryAnnoOrb2314.VS.Sensus22931/Target')
# PlantQueryNA <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 1, q_pids=1) # for test
PlantQuery5 <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 5,decoy = TRUE)
PlantQuery10 <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 10,decoy = TRUE)
PlantQuery0.005 <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 0.005,decoy = TRUE)
PlantQuery0.028 <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 0.028,decoy = TRUE)
PlantQuery0.050 <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 0.050,decoy = TRUE)
PlantQueryNA <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7,decoy = TRUE)
setwd('e:/dynamic/DynamicTol/output/LibrarySearch/Plant/QueryAnnoOrb2314.VS.Sensus22931/Decoy')
PlantDecoySirius <- SpectralMatching(q_dbPth=Pth_DecoySirius, l_dbPth=Pth_Library, cores = 7,decoy = FALSE)




#### QueryAnnoOrb2314.VS.OrbRemovePosPlusAlignedInstrKnown0.2_86662
Pth_Query <- 'e:/dynamic/DynamicTol/input/Sqlite/QueryAnnoOrb.sqlite' # as query
Pth_DecoySirius <- 'e:/dynamic/DynamicTol/input/Sqlite/QueryAnnoOrbSirius.sqlite' # as query
Pth_Library <- 'e:/dynamic/DynamicTol/input/Sqlite/OrbRemovePosPlusAlignedInstrKnown0.2.sqlite' #
setwd('e:/dynamic/DynamicTol/output/LibrarySearch/Plant/QueryAnnoOrb2314.VS.OrbRemovePosPlusAlignedInstrKnown0.2_86662/Target')
# PlantQueryNA <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 1, q_pids=1) # for test
PlantQuery5 <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 5,decoy = TRUE)
PlantQuery10 <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 10,decoy = TRUE)
PlantQuery0.005 <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 0.005,decoy = TRUE)
PlantQuery0.028 <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 0.028,decoy = TRUE)
PlantQuery0.050 <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 0.050,decoy = TRUE)
PlantQueryNA <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7,decoy = TRUE)
setwd('e:/dynamic/DynamicTol/output/LibrarySearch/Plant/QueryAnnoOrb2314.VS.OrbRemovePosPlusAlignedInstrKnown0.2_86662/Decoy')
PlantDecoySirius <- SpectralMatching(q_dbPth=Pth_DecoySirius, l_dbPth=Pth_Library, cores = 7,decoy = FALSE)

#### QueryAnnoOrb2314.VS.Sensus0.2_13101
Pth_Query <- 'e:/dynamic/DynamicTol/input/Sqlite/QueryAnnoOrb.sqlite' # as query
Pth_DecoySirius <- 'e:/dynamic/DynamicTol/input/Sqlite/QueryAnnoOrbSirius.sqlite' # as query
Pth_Library <- 'e:/dynamic/DynamicTol/input/Sqlite/GnpsSensus0.2.sqlite' #
setwd('e:/dynamic/DynamicTol/output/LibrarySearch/Plant/QueryAnnoOrb2314.VS.Sensus0.2_13101/Target')
# PlantQueryNA <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 1, q_pids=1) # for test
PlantQuery5 <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 5,decoy = TRUE)
PlantQuery10 <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 10,decoy = TRUE)
PlantQuery0.005 <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 0.005,decoy = TRUE)
PlantQuery0.028 <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 0.028,decoy = TRUE)
PlantQuery0.050 <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 0.050,decoy = TRUE)
PlantQueryNA <- SpectralMatching(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7,decoy = TRUE)
setwd('e:/dynamic/DynamicTol/output/LibrarySearch/Plant/QueryAnnoOrb2314.VS.Sensus0.2_13101/Decoy')
PlantDecoySirius <- SpectralMatching(q_dbPth=Pth_DecoySirius, l_dbPth=Pth_Library, cores = 7,decoy = FALSE)



#### Prepare the consensus database based on OrbRemovePosPlusAlignedInstr
## Start creating conSensus library_spectra
MetaL <- library_spectra_meta_OrbRemovePosPlusAlignedInstr
MetaL$inchikey_14 <- sub(MetaL$inchikey_id, pattern = "-.*",replacement = "",perl = TRUE)
MetaL[,inchikey_14_precursor_mz:=paste0(inchikey_14,'_',precursor_mz)]
library_spectra_Sensus <- data.table()
for (inchi in unique(MetaL$inchikey_14_precursor_mz)) {
  ## Extract the meta id related to each unique inchikey
  # inchi <- 'GWCQNKRMTGVYIZ_357.197'
  print(inchi)
  Selected <- MetaL[inchikey_14_precursor_mz==inchi]
  ## Extract MS2 spectra related to each meta id
  ListMS2 <- list()
  for (j in 1:nrow(Selected)){
    # j <- 1
    SelectedMS2 <- library_spectra_OrbRemovePosPlusAlignedInstr[library_spectra_meta_id==Selected[j,]$id]
    ListMS2[[j]]<- SelectedMS2[,c('mz','i')] %>% setnames(.,'i','intensity') %>% as.matrix(.)
  }
  ## Combine peaks
  Combined <- combinePeaks(ListMS2,ppm = 10, peaks = 'intersect',minProp = 0.5)
  ## Add meta information
  CombinedMeta <- as.data.table(Combined) %>% setnames(.,'intensity','i') %>%
    .[,library_spectra_meta_id:=min(Selected$id)] %>% .[,inchikey_14_precursor_mz:=inchi]
  library_spectra_Sensus<- rbind(library_spectra_Sensus, CombinedMeta)
}
unique(library_spectra_Sensus,by='inchikey_14_precursor_mz')
MetaL$inchikey_14_precursor_mz

## Creat Sensus.sqlite
con_Sensus <- DBI::dbConnect(RSQLite::SQLite(),'e:/dynamic/DynamicTol/input/Sqlite/GnpsSensus.sqlite')
library_spectra_meta_Sensus <- MetaL
setorder(library_spectra_meta_Sensus,id) # small to big
library_spectra_meta_Sensus <- unique(library_spectra_meta_Sensus, by= 'inchikey_14_precursor_mz')
DBI::dbWriteTable(con_Sensus, name = "library_spectra_meta", value = library_spectra_meta_Sensus, overwrite = T)
library_spectra_source_Sensus <- data.frame(id=1,
                                            name=paste('ConSensus Database',  format(Sys.time(), "%Y-%m-%d-%I%M%S"), sep='-'),
                                            parsing_software=paste('DBI::dbWriteTable'))
DBI::dbWriteTable(con_Sensus, name='library_spectra_source', value=library_spectra_source_Sensus, overwrite = T)
metab_compound_Sensus <- library_spectra_meta_OrbRemovePosPlusAlignedInstr[inchikey_id %in% MetaL$inchikey_id, ]
DBI::dbWriteTable(con_Sensus, name='metab_compound',value=metab_compound_Sensus, overwrite = T)
DBI::dbWriteTable(con_Sensus, name='library_spectra', value=library_spectra_Sensus, overwrite=T)
DBI::dbDisconnect(con_Sensus)



#### Prepare the consensus database based on OrbRemovePosPlusAlignedInstrKnown0.2
con_OrbAlignedInstrKnown0.2 <- DBI::dbConnect(RSQLite::SQLite(),'e:/dynamic/DynamicTol/input/Sqlite/OrbRemovePosPlusAlignedInstrKnown0.2.sqlite')
library_spectra_meta_OrbAlignedInstrKnown0.2 <- con_OrbAlignedInstrKnown0.2 %>% dplyr::tbl("library_spectra_meta") %>% dplyr::collect() %>% as.data.table(.)
library_spectra_OrbAlignedInstrKnown0.2 <- con_OrbAlignedInstrKnown0.2 %>% dplyr::tbl("library_spectra") %>% dplyr::collect() %>% as.data.table(.)

DBI::dbDisconnect(con_OrbAlignedInstrKnown0.2)

## Start creating conSensus library_spectra
MetaL <- library_spectra_meta_OrbAlignedInstrKnown0.2
MetaL$inchikey_14 <- sub(MetaL$inchikey_id, pattern = "-.*",replacement = "",perl = TRUE)
MetaL[,inchikey_14_precursor_mz:=paste0(inchikey_14,'_',precursor_mz)]
library_spectra_Sensus <- data.table()
for (inchi in unique(MetaL$inchikey_14_precursor_mz)) {
  ## Extract the meta id related to each unique inchikey
  # inchi <- 'GWCQNKRMTGVYIZ_357.197'
  print(inchi)
  Selected <- MetaL[inchikey_14_precursor_mz==inchi]
  ## Extract MS2 spectra related to each meta id
  ListMS2 <- list()
  for (j in 1:nrow(Selected)){
    # j <- 1
    SelectedMS2 <- library_spectra_OrbAlignedInstrKnown0.2[library_spectra_meta_id==Selected[j,]$id]
    ListMS2[[j]]<- SelectedMS2[,c('mz','i')] %>% setnames(.,'i','intensity') %>% as.matrix(.)
  }
  ## Combine peaks
  Combined <- combinePeaks(ListMS2,ppm = 10, peaks = 'intersect',minProp = 0.5)
  ## Add meta information
  CombinedMeta <- as.data.table(Combined) %>% setnames(.,'intensity','i') %>%
    .[,library_spectra_meta_id:=min(Selected$id)] %>% .[,inchikey_14_precursor_mz:=inchi]
  library_spectra_Sensus<- rbind(library_spectra_Sensus, CombinedMeta)
}
unique(library_spectra_Sensus,by='inchikey_14_precursor_mz')
MetaL$inchikey_14_precursor_mz

## Create Sensus.sqlite
con_Sensus0.2 <- DBI::dbConnect(RSQLite::SQLite(),'e:/dynamic/DynamicTol/input/Sqlite/Sensus0.2.sqlite')
library_spectra_meta_Sensus <- MetaL
setorder(library_spectra_meta_Sensus,id) # small to big
library_spectra_meta_Sensus <- unique(library_spectra_meta_Sensus, by= 'inchikey_14_precursor_mz')
DBI::dbWriteTable(con_Sensus0.2, name = "library_spectra_meta", value = library_spectra_meta_Sensus, overwrite = T)
library_spectra_source_Sensus <- data.frame(id=1,
                                            name=paste('ConSensus Database',  format(Sys.time(), "%Y-%m-%d-%I%M%S"), sep='-'),
                                            parsing_software=paste('DBI::dbWriteTable'))
DBI::dbWriteTable(con_Sensus0.2, name='library_spectra_source', value=library_spectra_source_Sensus, overwrite = T)
metab_compound_Sensus <- library_spectra_meta_OrbRemovePosPlusAlignedInstr[inchikey_id %in% MetaL$inchikey_id, ]
DBI::dbWriteTable(con_Sensus0.2, name='metab_compound',value=metab_compound_Sensus, overwrite = T)
DBI::dbWriteTable(con_Sensus0.2, name='library_spectra', value=library_spectra_Sensus, overwrite=T)
DBI::dbDisconnect(con_Sensus0.2)



#### Create the decoy database by sirius
## Extract the information from the library
# Pth_QueryAnnoOrb <- 'e:/dynamic/DynamicTol/input/Sqlite/QueryAnnoOrb.sqlite' # original
Pth_QueryAnnoOrb <- 'e:/dynamic/DynamicTol/input/Sqlite/QueryAnnoOrb0.002.sqlite' # with more close m/z and rt of 0.002
con_QueryAnnoOrb <- DBI::dbConnect(RSQLite::SQLite(),Pth_QueryAnnoOrb)
library_spectra_meta_QueryAnnoOrb <- con_QueryAnnoOrb %>% dplyr::tbl("library_spectra_meta") %>% dplyr::collect() %>% as.data.table(.)
library_spectra_QueryAnnoOrb <- con_QueryAnnoOrb %>% dplyr::tbl("library_spectra") %>% dplyr::collect() %>% as.data.table(.)
metab_compound_QueryAnnoOrb <- con_QueryAnnoOrb %>% dplyr::tbl("metab_compound") %>% dplyr::collect() %>% as.data.table(.)
library_spectra_source_QueryAnnoOrb <- con_QueryAnnoOrb %>% dplyr::tbl("library_spectra_source") %>% dplyr::collect() %>% as.data.table(.)
DBI::dbDisconnect(con_QueryAnnoOrb)

Pth_OrbRemovePosPlusAlignedInstr <- 'e:/dynamic/DynamicTol/input/Sqlite/OrbRemovePosPlusAlignedInstr.sqlite'
con_OrbRemovePosPlusAlignedInstr <- DBI::dbConnect(RSQLite::SQLite(),Pth_OrbRemovePosPlusAlignedInstr)
library_spectra_meta_OrbRemovePosPlusAlignedInstr <- con_OrbRemovePosPlusAlignedInstr %>% dplyr::tbl("library_spectra_meta") %>% dplyr::collect() %>% as.data.table(.)
metab_compound_OrbRemovePosPlusAlignedInstr <- con_OrbRemovePosPlusAlignedInstr %>% dplyr::tbl("metab_compound") %>% dplyr::collect() %>% as.data.table(.)
DBI::dbDisconnect(con_OrbRemovePosPlusAlignedInstr)

DT <- library_spectra_meta_QueryAnnoOrb
DT[,inchikey_id_14:=sub(inchikey_id,pattern = "-.*",replacement = "",perl = TRUE)]
metab_compound_OrbRemovePosPlusAlignedInstr_Uni <- metab_compound_OrbRemovePosPlusAlignedInstr[, c('inchikey_id','molecular_formula')] %>%
  .[,inchikey_id_14:=sub(inchikey_id,pattern = "-.*",replacement = "",perl = TRUE)] %>% unique(.,by='inchikey_id_14')
DT <- left_join(DT, metab_compound_OrbRemovePosPlusAlignedInstr_Uni, by='inchikey_id_14')

## Open the MS file for writing
InputMs <- "e:/dynamic/DynamicTol/input/ms/QueryAnnoOrb0.002.ms"
msFile <- file(InputMs, "w")
for (Ind in DT$id) {
  # Ind <- 1
  # names(DT)
  print(DT[id==Ind,]$name)
  spectrum <- list(
    compound = paste0(DT[id==Ind,]$inchikey_id_14, "_",Ind),
    formula = DT[id==Ind,]$molecular_formula,
    parentmass = DT[id==Ind,]$precursor_mz,
    ionization = DT[id==Ind,]$precursor_type,
    mz2Values = library_spectra_QueryAnnoOrb[library_spectra_meta_id==Ind,]$mz,
    intensity2Values = library_spectra_QueryAnnoOrb[library_spectra_meta_id==Ind,]$i
  )
  writeMSSpectrum(msFile,spectrum)
}
close(msFile)

## Analyze the .ms database using sirius
OutDir <- 'e:/dynamic/DynamicTol/output/Sirius/QueryAnnoOrb0.002/'
dir.create(OutDir) # need to create the folder first, otherwise will have no folders
sirius_cmd <- paste('"e:/Program Files/sirius/sirius.exe"',
                    '-i',InputMs,
                    '-o', OutDir,
                    # '--processors 5',
                    'formula -p orbitrap passatutto write-summaries')
system(sirius_cmd)

## Collect the decoys from the sirius to replace the library_spectra
formula_identifications <- fread('e:/dynamic/DynamicTol/output/Sirius/QueryAnnoOrb0.002/formula_identifications.tsv')
decoys.rbind = data.table()
for (ind in formula_identifications$id) {
  # ind <- '501_QueryAnnoOrb0.002_VOCJFEUCJDTXSA_15437'
  print(ind)
  MetaId <- sub('.*\\_','',ind) %>% as.numeric(.)
  Inchikey14 <- sub('.*QueryAnnoOrb0.002\\_','',ind) %>% sub('\\_.','',.)
  # NRow <- sub('_OrbRemovePos.*','',ind) %>% as.numeric(.)
  DecoyDir <- dir(paste0('e:/dynamic/DynamicTol/output/Sirius/QueryAnnoOrb0.002/', ind,'/decoys/'), full.names = TRUE, recursive = TRUE)
  if (length(DecoyDir)>0) {
    decoy <- fread(DecoyDir) %>% setnames(.,'rel.intensity','i') %>% .[,c('mz','i')] %>% # extract the MS2
      .[,library_spectra_meta_id:=   MetaId] %>% #  Add the meta_id index
      .[,inchikey_14:= Inchikey14] #  Add the meta_id index and inchikey_14_precursor_mz
    decoys.rbind <- rbind(decoys.rbind, decoy)
  }
}

## write back into the .sqlite database
Pth_QueryAnnoOrbSirius <- 'e:/dynamic/DynamicTol/input/Sqlite/QueryAnnoOrb0.002Sirius.sqlite'
con_QueryAnnoOrbSirius <- DBI::dbConnect(RSQLite::SQLite(),Pth_QueryAnnoOrbSirius)
DBI::dbWriteTable(con_QueryAnnoOrbSirius, name='library_spectra_meta', value=library_spectra_meta_QueryAnnoOrb)
DBI::dbWriteTable(con_QueryAnnoOrbSirius, name='metab_compound',value=metab_compound_QueryAnnoOrb)
DBI::dbWriteTable(con_QueryAnnoOrbSirius, name='library_spectra_source', value=library_spectra_source_QueryAnnoOrb)
DBI::dbWriteTable(con_QueryAnnoOrbSirius, name='library_spectra', value=decoys.rbind)
DBI::dbDisconnect(con_QueryAnnoOrbSirius)

