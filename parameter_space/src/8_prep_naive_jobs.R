# ##############################################################################
#
##  Script to produce scripts to train the naive models
#
# ##############################################################################

.libPaths(c('/g/scb2/zeller/SHARED/software/R/3.5', .libPaths()))

# packages
library("tidyverse")
library("here")

# ##############################################################################
# functions to expand parameters
# function to prepare submission scripts
prep.submission.script <- function(dataset, type, case, ml.method){
  
  script.file <- here("parameter_space", "submission_temp",
                      paste0('script_', dataset, '_', type, '_', case, 
                             '_', ml.method, '.sh'))
  
  lines <- c(
    "#!/bin/bash", "#SBATCH -A zeller",
    "#SBATCH -t 10:00:00", "#SBATCH --mem 3G",
    paste0("#SBATCH -o /scratch/jawirbel/siamcat/job.naive.", dataset, 
           '_', type, '_', case, '_',
           ml.method,  ".out"),
    paste0("#SBATCH -e /scratch/jawirbel/siamcat/job.naive.", dataset, 
           '_', type, '_', case, '_',
           ml.method,  ".err"),
    "", "module load R/3.5.1-foss-2017b-X11-20171023", "",
    paste0("Rscript ../src/9_train_best_models.R --dataset ", dataset, 
           " --type ", type, " --case ", case, ' --ml_method ', ml.method))
  
  tmp <- file(script.file)
  writeLines(lines, con=tmp)
  close(tmp)
}

# ##############################################################################
# prepare submission scripts

# load info about datasets
motu.tasks <- read_tsv(here('parameter_space', 'data_info', 
                            'all_tasks.tsv')) %>% 
  filter(type=='mOTUs2')


for (i in seq_len(nrow(motu.tasks))){
  for (ml.method in c('enet', 'lasso', 'enet-0.5', 'randomForest')){
    prep.submission.script(motu.tasks$dataset.id[i],
                           motu.tasks$type[i],
                           motu.tasks$case[i], 
                           ml.method)
  }
}
