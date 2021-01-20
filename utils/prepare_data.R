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
# Download SIAMCAT data from Zenodo
# TODO download zipped metadata and unzip
# download zipped mOTU profiles and unzip
# download zipped eggNOG profiles and unzip

datasets <- c('CN-T2D', 'SE-T2D', 'DE-PD', 'metaHIT-IBD', 'CN-RA', 'FR-CRC',
              'AT-CRC', 'CN-CRC', 'US-CRC', 'DE-CRC', 'SP-NAFLD', 'CN-AS',
              'CN-ACVD', 'Lewis-IBD', 'He-IBD', 'Franzosa-IBD', 'JP-CRC',
              'US-NAFLD', 'KZ-MS', 'HMP2-IBD')
for (d in datasets){
  # metadata
  if (!file.exists(here('data', 'meta', paste0(d, '.tsv')))){
    file.move(here('temp', paste0(d, '.tsv')),
              here('data', 'meta'))
  }
  # mOTUs
  if (!file.exists(here('data', 'features', paste0(d, '_motus.tsv')))){
    file.move(here('temp', paste0(d, '_motus.tsv')),
              here('data', 'features'))
  }
  # eggNOG
  if (!file.exists(here('data', 'features', paste0(d, '_eggNOG.tsv')))){
    file.move(here('temp', paste0(d, '_eggNOG.tsv')),
              here('data', 'features'))
  }
}

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
