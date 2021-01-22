# ##############################################################################
#
##  Script to produce scripts to train the control-augmented models
#
# ##############################################################################

.libPaths(c('/g/scb2/zeller/SHARED/software/R/3.5', .libPaths()))

# packages
library("tidyverse")
library("here")

# ##############################################################################
# functions to expand parameters
# function to prepare submission scripts
prep.submission.script <- function(dataset, type, case,
                                   ctr_type, ml.method){
  
  script.file <- here("control_augmentation", "submission_temp",
                      paste0('script_', dataset, '_', type, '_', case, 
                             '_',  ctr_type, '_', ml.method, '.sh'))
  
  lines <- c(
    "#!/bin/bash", "#SBATCH -A zeller",
    "#SBATCH -t 10:00:00", "#SBATCH --mem 3G",
    paste0("#SBATCH -o /scratch/jawirbel/siamcat/job.augm.", dataset, 
           '_', type, '_', case, '_', ctr_type, '_',
           ml.method,  ".out"),
    paste0("#SBATCH -e /scratch/jawirbel/siamcat/job.augm.", dataset, 
           '_', type, '_', case, '_', ctr_type, '_',
           ml.method,  ".err"),
    "", "module load R/3.5.1-foss-2017b-X11-20171023", "",
    paste0("Rscript ../src/train_control_augmented.R --dataset ", dataset, 
           " --type ", type, " --case ", case,
           ' --ctr_type ', ctr_type, ' --ml_method ', ml.method))
  
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
    for (ctr_type in c('cohort_5', 'cohort_2', 'other_ctr', 'random')){
      for (ml.method in c('enet', 'lasso', 'enet-0.5', 'randomForest')){
        prep.submission.script(motu.tasks$dataset.id[i],
                               motu.tasks$type[i],
                               motu.tasks$case[i], ctr_type,
                               ml.method)
      }
    }
  }
}
