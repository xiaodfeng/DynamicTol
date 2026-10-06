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
## NOTE: F_MakeXcorrDecoyDB() pools all 150 shifted copies into ONE spectrum per query.
## Searching that pooled spectrum gives the cosine of the pooled spectrum, which is NOT
## the per-shift average of Eq. (6)-(7) that the manuscript uses for Xcorr and
## XcorrCutoff. It is kept for reference only; the XcorrCutoff decoys used in Fig. 7
## are produced by the SpectralMatching2(..., DecoyCutoff = 1) run below.
F_MakeXcorrDecoyDB(inPth='d:/OneDrive/github/dynamic/input/Sqlite/Query500Pos500Neg_MonaOrb_1117.sqlite',
                   outPth='d:/OneDrive/github/dynamic/input/Sqlite/Query500Pos500Neg_MonaOrb_1117_XcorrCutoff.sqlite',
                   noise_cutoff=1)
F_MakeXcorrDecoyDB(inPth='d:/OneDrive/github/dynamic/input/Sqlite/Query500Pos500Neg_MonaOrb_1117.sqlite',
                   outPth='d:/OneDrive/github/dynamic/input/Sqlite/Query500Pos500Neg_MonaOrb_1117_Xcorr.sqlite',
                   noise_cutoff=0)
## Perform the library search for the XcorrCutoff decoy (manuscript definition):
## ORIGINAL query vs consensus library; peaks < 1% of the base peak are removed only
## before the 150 shifted decoy spectra are generated, so the target score is identical
## to the Xcorr variant and only the decoy changes. MonaOrb.R section 12 reads
## decoy.mean / decoy.mean.dpc from this folder.
Pth_Query <- "d:/OneDrive/github/dynamic/input/Sqlite/Query500Pos500Neg_MonaOrb_1117.sqlite"
Pth_Library <- "d:/OneDrive/github/dynamic/input/Sqlite/LibraryNoNeg_MonaOrb_1117_Sensus.sqlite"
Dir <- "d:/OneDrive/github/dynamic/output/LibrarySearch/MonaOrb/20260826/QueryXcorrCutoff.VS.Sensus"
dir.create(Dir, recursive = TRUE)
setwd(Dir)
Test <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, mztol = NA, cores = 2,q_pids = 115,decoy = TRUE,DecoyCutoff = 1, write = FALSE)
Test
XcorrNA <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, mztol = NA,cores = 4,decoy = TRUE, DecoyCutoff = 1)
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
#  separate DT_Sensus / DT_Query objects. 'DT' was reassigned to the QUERY meta
# table below, but the library decoy collection later indexed DT[NRow,] and therefore
# attached query meta ids to the library decoys.
DT <- DT_Sensus <- library_spectra_meta_Sensus
names(DT)
metab_compound_Sensus # 4747
metab_compound_Sensus_Uni <- unique(metab_compound_Sensus,by="inchikey_id") # 4247
DT <- DT_Sensus <- left_join(DT,metab_compound_Sensus_Uni[,c("inchikey_id","molecular_formula")],by="inchikey_id")
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
DT <- DT_Query <- library_spectra_meta_Query
names(DT)
metab_compound_Query # 1000
metab_compound_Query_Uni <- unique(metab_compound_Query,by="inchikey_id") # 1000
DT <- DT_Query <- left_join(DT,metab_compound_Query_Uni[,c("inchikey_id","molecular_formula")],by="inchikey_id")
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
decoys_Sensus = data.table()
for (ind in formula_identifications$id) {
  # ind <- "5490_LibraryNoNeg_MonaOrb_1117_Sensus_Khayanthone"
  print(ind)
  NRow <- sub('_LibraryNoNeg.*','',ind) %>% as.numeric(.)
  DecoyDir <- dir(paste0('D:/github/sirius/output/LibraryNoNeg_MonaOrb_1117_Sensus/', ind,'/decoys/'), full.names = TRUE, recursive = TRUE)
  if (length(DecoyDir)>0) {
    decoy <- fread(DecoyDir) %>% setnames(.,'rel.intensity','i') %>% .[,c('mz','i')] %>% # extract the MS2
      .[,library_spectra_meta_id:= DT_Sensus[NRow,]$id] #%>%   Add the meta_id index
    # .[,inchikey_14_precursor_mz:= DT[NRow,]$inchikey_14_precursor_mz] #  Add the meta_id index and inchikey_14_precursor_mz
    decoys_Sensus <- rbind(decoys_Sensus, decoy)
  }
}
## For query
formula_identifications <- fread('D:/github/sirius/output/Query500Pos500Neg_MonaOrb_1117/formula_identifications.tsv')
decoys_Query = data.table()
for (ind in formula_identifications$id) {
  # ind <- '6_Query500Pos500Neg_MonaOrb_1117_methyl5Z-5-ethylidene-4-2-2R3S4S5R6R'
  print(ind)
  NRow <- sub('_Query.*','',ind) %>% as.numeric(.)
  DecoyDir <- dir(paste0('D:/github/sirius/output/Query500Pos500Neg_MonaOrb_1117/', ind,'/decoys/'), full.names = TRUE, recursive = TRUE)
  if (length(DecoyDir)>0) {
    decoy <- fread(DecoyDir) %>% setnames(.,'rel.intensity','i') %>% .[,c('mz','i')] %>% # extract the MS2
      .[,library_spectra_meta_id:= DT_Query[NRow,]$id] # %>% #  Add the meta_id index
    # .[,inchikey_14_precursor_mz:= DT[NRow,]$inchikey_14_precursor_mz] #  Add the meta_id index and inchikey_14_precursor_mz
    decoys_Query <- rbind(decoys_Query, decoy)
  }
}

#### write back into the .sqlite database
## For library
Pth_SensusDecoy <- 'd:/OneDrive/github/dynamic/input/Sqlite/LibraryNoNeg_MonaOrb_1117_Sensus_Decoy.sqlite'
con_SensusDecoy <- DBI::dbConnect(RSQLite::SQLite(),Pth_SensusDecoy)
DBI::dbWriteTable(con_SensusDecoy, name='library_spectra_meta', value=library_spectra_meta_Sensus)
DBI::dbWriteTable(con_SensusDecoy, name='metab_compound',value=metab_compound_Sensus)
DBI::dbWriteTable(con_SensusDecoy, name='library_spectra_source', value=library_spectra_source_Sensus)
DBI::dbWriteTable(con_SensusDecoy, name='library_spectra', value=decoys.Sensus)
DBI::dbDisconnect(con_SensusDecoy)
## For query
Pth_QueryDecoy <- 'd:/OneDrive/github/dynamic/input/Sqlite/Query500Pos500Neg_MonaOrb_1117_Decoy.sqlite'
con_QueryDecoy <- DBI::dbConnect(RSQLite::SQLite(),Pth_QueryDecoy)
DBI::dbWriteTable(con_QueryDecoy, name='library_spectra_meta', value=library_spectra_meta_Query)
DBI::dbWriteTable(con_QueryDecoy, name='metab_compound',value=metab_compound_Query)
DBI::dbWriteTable(con_QueryDecoy, name='library_spectra_source', value=library_spectra_source_Query)
DBI::dbWriteTable(con_QueryDecoy, name='library_spectra', value=decoys.Query)
DBI::dbDisconnect(con_QueryDecoy)

#### Perform the library search using the created sensus decoy
## For library decoy
Pth_Query <- "d:/OneDrive/github/dynamic/input/Sqlite/Query500Pos500Neg_MonaOrb_1117.sqlite"
Pth_Library <- "d:/OneDrive/github/dynamic/input/Sqlite/LibraryNoNeg_MonaOrb_1117_Sensus_Decoy.sqlite"
Dir <- "d:/OneDrive/github/dynamic/output/LibrarySearch/MonaOrb/20260826/Query.VS.SiriusDecoy"
dir.create(Dir, recursive = TRUE)
setwd(Dir)
Test <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library,
                          mztol = NA, cores = 1,q_pids = 115,decoy = FALSE,write=FALSE)
setorder(Test,-dpc)
Test
SiriusNA <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library,cores = 5)
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
# Test <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 2,mztol=NA,
#                          q_pids=100525, l_pids = 187578,decoy = TRUE,write=FALSE) # for test
# Test
PlantQuery5 <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 5,decoy = TRUE)
PlantQuery10 <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 10,decoy = TRUE)
PlantQuery0.005 <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 0.005,decoy = TRUE)
PlantQuery0.028 <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 0.028,decoy = TRUE)
PlantQuery0.050 <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 0.050,decoy = TRUE)
PlantQueryNA <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7,decoy = TRUE)
setwd('e:/dynamic/DynamicTol/output/LibrarySearch/Plant/QueryAnnoOrb0.002.VS.OrbRemovePosPlusAlignedInstrKnown0.2_86662/Decoy')
PlantDecoySirius <- SpectralMatching2(q_dbPth=Pth_DecoySirius, l_dbPth=Pth_Library, cores = 7,decoy = FALSE)
#### QueryAnnoOrb0.002.VS.Sensus0.2_13101
Pth_Query <- 'e:/dynamic/DynamicTol/input/Sqlite/QueryAnnoOrb0.002.sqlite' # as query
Pth_DecoySirius <- 'e:/dynamic/DynamicTol/input/Sqlite/QueryAnnoOrb0.002Sirius.sqlite' # as query
Pth_Library <- 'e:/dynamic/DynamicTol/input/Sqlite/GnpsSensus0.2.sqlite' #
setwd('e:/dynamic/DynamicTol/output/LibrarySearch/Plant/QueryAnnoOrb0.002.VS.Sensus0.2_13101/Target')
PlantQueryNA <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 1, q_pids=100525) # for test

setwd('e:/dynamic/DynamicTol/output/LibrarySearch/Plant/QueryAnnoOrb0.002.VS.Sensus0.2_13101/Decoy')
PlantDecoySirius <- SpectralMatching2(q_dbPth=Pth_DecoySirius, l_dbPth=Pth_Library, cores = 7,decoy = FALSE)


#### QueryAnnoOrb2314.VS.OrbRemovePosPlusAlignedInstr185865
Pth_Query <- 'e:/dynamic/DynamicTol/input/Sqlite/QueryAnnoOrb.sqlite' # as query
Pth_DecoySirius <- 'e:/dynamic/DynamicTol/input/Sqlite/QueryAnnoOrbSirius.sqlite' # as query
Pth_Library <- 'e:/dynamic/DynamicTol/input/Sqlite/OrbRemovePosPlusAlignedInstr.sqlite' #
setwd('e:/dynamic/DynamicTol/output/LibrarySearch/Plant/QueryAnnoOrb2314.VS.OrbRemovePosPlusAlignedInstr185865/Target')
# PlantQueryNA <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 1, q_pids=1) # for test
PlantQuery5 <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 5,decoy = TRUE)
PlantQuery10 <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 10,decoy = TRUE)
PlantQuery0.005 <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 0.005,decoy = TRUE)
PlantQuery0.028 <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 0.028,decoy = TRUE)
PlantQuery0.050 <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 0.050,decoy = TRUE)
PlantQueryNA <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7,decoy = TRUE)
setwd('e:/dynamic/DynamicTol/output/LibrarySearch/Plant/QueryAnnoOrb2314.VS.OrbRemovePosPlusAlignedInstr185865/Decoy')
PlantDecoySirius <- SpectralMatching2(q_dbPth=Pth_DecoySirius, l_dbPth=Pth_Library, cores = 7,decoy = FALSE)

#### QueryAnnoOrb2314.VS.Sensus22931
Pth_Query <- 'e:/dynamic/DynamicTol/input/Sqlite/QueryAnnoOrb.sqlite' # as query
Pth_DecoySirius <- 'e:/dynamic/DynamicTol/input/Sqlite/QueryAnnoOrbSirius.sqlite' # as query
Pth_Library <- 'e:/dynamic/DynamicTol/input/Sqlite/GnpsSensus.sqlite' #
setwd('e:/dynamic/DynamicTol/output/LibrarySearch/Plant/QueryAnnoOrb2314.VS.Sensus22931/Target')
# PlantQueryNA <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 1, q_pids=1) # for test
PlantQuery5 <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 5,decoy = TRUE)
PlantQuery10 <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 10,decoy = TRUE)
PlantQuery0.005 <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 0.005,decoy = TRUE)
PlantQuery0.028 <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 0.028,decoy = TRUE)
PlantQuery0.050 <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 0.050,decoy = TRUE)
PlantQueryNA <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7,decoy = TRUE)
setwd('e:/dynamic/DynamicTol/output/LibrarySearch/Plant/QueryAnnoOrb2314.VS.Sensus22931/Decoy')
PlantDecoySirius <- SpectralMatching2(q_dbPth=Pth_DecoySirius, l_dbPth=Pth_Library, cores = 7,decoy = FALSE)

#### QueryAnnoOrb2314.VS.OrbRemovePosPlusAlignedInstrKnown0.2_86662
Pth_Query <- 'e:/dynamic/DynamicTol/input/Sqlite/QueryAnnoOrb.sqlite' # as query
Pth_DecoySirius <- 'e:/dynamic/DynamicTol/input/Sqlite/QueryAnnoOrbSirius.sqlite' # as query
Pth_Library <- 'e:/dynamic/DynamicTol/input/Sqlite/OrbRemovePosPlusAlignedInstrKnown0.2.sqlite' #
setwd('e:/dynamic/DynamicTol/output/LibrarySearch/Plant/QueryAnnoOrb2314.VS.OrbRemovePosPlusAlignedInstrKnown0.2_86662/Target')
# PlantQueryNA <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 1, q_pids=1) # for test
PlantQuery5 <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 5,decoy = TRUE)
PlantQuery10 <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 10,decoy = TRUE)
PlantQuery0.005 <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 0.005,decoy = TRUE)
PlantQuery0.028 <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 0.028,decoy = TRUE)
PlantQuery0.050 <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 0.050,decoy = TRUE)
PlantQueryNA <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7,decoy = TRUE)
setwd('e:/dynamic/DynamicTol/output/LibrarySearch/Plant/QueryAnnoOrb2314.VS.OrbRemovePosPlusAlignedInstrKnown0.2_86662/Decoy')
PlantDecoySirius <- SpectralMatching2(q_dbPth=Pth_DecoySirius, l_dbPth=Pth_Library, cores = 7,decoy = FALSE)

#### QueryAnnoOrb2314.VS.Sensus0.2_13101
Pth_Query <- 'e:/dynamic/DynamicTol/input/Sqlite/QueryAnnoOrb.sqlite' # as query
Pth_DecoySirius <- 'e:/dynamic/DynamicTol/input/Sqlite/QueryAnnoOrbSirius.sqlite' # as query
Pth_Library <- 'e:/dynamic/DynamicTol/input/Sqlite/GnpsSensus0.2.sqlite' #
setwd('e:/dynamic/DynamicTol/output/LibrarySearch/Plant/QueryAnnoOrb2314.VS.Sensus0.2_13101/Target')
# PlantQueryNA <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 1, q_pids=1) # for test
PlantQuery5 <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 5,decoy = TRUE)
PlantQuery10 <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 10,decoy = TRUE)
PlantQuery0.005 <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 0.005,decoy = TRUE)
PlantQuery0.028 <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 0.028,decoy = TRUE)
PlantQuery0.050 <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7, mztol = 0.050,decoy = TRUE)
PlantQueryNA <- SpectralMatching2(q_dbPth=Pth_Query, l_dbPth=Pth_Library, cores = 7,decoy = TRUE)
setwd('e:/dynamic/DynamicTol/output/LibrarySearch/Plant/QueryAnnoOrb2314.VS.Sensus0.2_13101/Decoy')
PlantDecoySirius <- SpectralMatching2(q_dbPth=Pth_DecoySirius, l_dbPth=Pth_Library, cores = 7,decoy = FALSE)

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

