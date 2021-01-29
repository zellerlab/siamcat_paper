# ##############################################################################
#
## Script to prepare all analyses by downloading the data from
##  the SIAMCAT paper Zenodo repository
##  the Pasolli et al. meta-analysis paper
##  the mOTUs database Zenodo repository
##  the Ocean microbiome project
#
# ##############################################################################

library("tidyverse")
library("here")
library("filesstrings")
library("tools")
library("yaml")

if (!dir.exists(here('temp'))){
  dir.create(here('temp'))
}

# ##############################################################################
# Download internal data from Zenodo

# download zipped metadata and unzip
meta.file <- 'https://zenodo.org/api/files/6583a5bd-ff64-4482-a2a0-f86f097eef2e/metadata.zip'
download.file(meta.file, destfile = here('temp', 'metadata.zip'))
ref.md5sum <- '8a962531b6697ad19f84cb3772da5061'
check.md5sum <- md5sum(here('temp', 'metadata.zip'))
if (!ref.md5sum==check.md5sum){
  stop('Something went wrong when downloading the metadata files!')
}
unzip(here('temp', 'metadata.zip'), exdir = here('temp'))
metadata.files <- list.files(here('temp'), pattern = 'meta_.*.tsv')
map(metadata.files, 
    .f = function(x){file.move(here("temp", x), here("data", 'meta'))})

# download zipped mOTU profiles and unzip
motus.file <- 'https://zenodo.org/api/files/6583a5bd-ff64-4482-a2a0-f86f097eef2e/motus.zip'
download.file(motus.file, destfile = here('temp', 'motus.zip'))
ref.md5sum <- 'd89154e4248ad283eaeadd418606b2e2'
check.md5sum <- md5sum(here('temp', 'motus.zip'))
if (!ref.md5sum==check.md5sum){
  stop('Something went wrong when downloading the mOTUs_v2 files!')
}
unzip(here('temp', 'motus.zip'), exdir = here('temp'))
motus.files <- list.files(here('temp'), pattern = 'motus.tsv')
map(motus.files, 
    .f = function(x){file.move(here("temp", x), 
                               here("data", 'features', 'motus'))})

# download zipped eggNOG profiles and unzip
eggNOG.file <- 'https://zenodo.org/api/files/6583a5bd-ff64-4482-a2a0-f86f097eef2e/eggNOG.zip'
download.file(eggNOG.file, destfile = here('temp', 'eggNOG.zip'))
ref.md5sum <- 'fad5f59ffdc6a20316eee17aa010128e'
check.md5sum <- md5sum(here('temp', 'eggNOG.zip'))
if (!ref.md5sum==check.md5sum){
  stop('Something went wrong when downloading the eggNOG files!')
}
unzip(here('temp', 'eggNOG.zip'), exdir = here('temp'))
eggNOG.files <- list.files(here('temp'), pattern = 'meta_.*.tsv')
map(eggNOG.files, 
    .f = function(x){file.move(here("temp", x), 
                               here("data", 'features', 'eggNOG'))})

# ##############################################################################
# mOTUs database (taxonomy file)
if (!file.exists(here('data', 'motus_taxonomy.tsv'))){
  db.file <- "https://zenodo.org/api/files/793e834a-a3c9-4d57-8ebe-f9a462b0dfbb/db_mOTU_v2.5.1.tar.gz"
  download.file(db.file, destfile = here('temp', 'motus_db.tar.gz'))
  ref.md5sum <- '6117b479706ff805fb1bdcd8642b70e4'
  check.md5sum <- md5sum(here('temp', 'motus_db.tar.gz'))
  if (!ref.md5sum==check.md5sum){
    stop('Something went wrong when downloading the mOTUs database!')
  }
  untar(here('temp', 'motus_db.tar.gz'), exdir = here('temp'))
  tax.ref <- read_tsv(here('temp', 'db_mOTU', 'db_mOTU_taxonomy_ref-mOTUs.tsv'))
  tax.meta <- read_tsv(here('temp', 'db_mOTU',
                            'db_mOTU_taxonomy_meta-mOTUs.tsv'))
  colnames(tax.ref)[2] <- 'mOTU_ID'
  colnames(tax.meta)[1] <- 'mOTU_ID'
  tax.full <- bind_rows(tax.ref, tax.meta)
  write_tsv(tax.full, here('data', 'motus_taxonomy.tsv'))
}

# ##############################################################################
# Download data from Pasolli et al.
if (!file.exists(here('data', 'pasolli', 'pasolli.txt.bz2'))){
  feat.file <- "https://bitbucket.org/CibioCM/metaml/raw/0bbf049bf5a9e62a2533a25c48d9a413e6420a8b/data/abundance.txt.bz2"
  download.file(feat.file, destfile = here('temp', 'pasolli.txt.bz2'))
  file.move(here('temp', 'pasolli.txt.bz2'),
            here('data', 'pasolli'))
}

# ##############################################################################
# Download Ocean data

# features
if (!file.exists(here('data', 'ocean', 'mitags_tab_genus.tsv.gz'))){
  feat.file <- 'https://www.ebi.ac.uk/biostudies/files/S-BSST297/u/OM-RGC_v2_taxonomic_profiles.tar.gz'
  download.file(feat.file, destfile = here('temp', 'ocean.tar.gz'))
  untar(here('temp', 'ocean.tar.gz'), exdir = here('temp'))
  file.move(here('temp', 'OM-RGC_v2_taxonomic_profiles',
                 'mitags_tab_genus.tsv.gz'),
            here('data', 'ocean'))
}
# metadata
if (!file.exists(here('data', 'ocean', 'ocean.xlsx'))){
  meta.file <- 'https://zenodo.org/record/3539258/files/Salazar_et_al_2019_Suppl_Info.xlsx?download=1'
  download.file(meta.file, destfile = here('temp', 'ocean.xlsx'))
  file.move(here('temp', 'ocean.xlsx'),
            here('data', 'ocean'))
}

# ##############################################################################
# clean the temp directory again
unlink(here('temp'), recursive = TRUE)
