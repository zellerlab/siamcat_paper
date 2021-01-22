# ##############################################################################
#
##  Script to produce the configerations for all jobs to explore the
##    parameter space
#
# ##############################################################################

# packages
library("tidyverse")
library("SIAMCAT")
library("yaml")
library("progress")
library("here")

# ##############################################################################
# functions to expand parameters
expand.filtering <- function(temp, parameter.space){
  df.temp <- tibble()
  for (filt.method in parameter.space$filtering$methods){
    temp$filt.method <- filt.method
    if (filt.method == 'pass'){
      temp$filt.cutoff <- NA
      df.temp <- bind_rows(df.temp, temp)
    } else if (filt.method == 'abundance;prevalence'){
      for (f1 in parameter.space$filtering$cutoffs$abundance){
        for (f2 in parameter.space$filtering$cutoffs$prevalence){
          temp$filt.cutoff <- paste0(as.character(f1),';', as.character(f2))
          df.temp <- bind_rows(df.temp, temp)
        }
      }
    }
    for (filt.cutoff in parameter.space$filtering$cutoffs[[filt.method]]){
      temp$filt.cutoff <- as.character(filt.cutoff)
      df.temp <- bind_rows(df.temp, temp)
    }
  }
  return(df.temp)
}

expand.normalization <- function(temp, parameter.space){
  df.temp <- tibble()
  for (norm.method in parameter.space$normalization$methods){
    temp$norm.method <- norm.method
    if (norm.method %in% c('rank', 'none')){
      temp$log.n0 <- NA
      if (!str_detect(norm.method, 'std')){
        temp$sd.min.q <- NA
        df.temp <- bind_rows(df.temp, temp)
      } else {
        for (sd.min.q in parameter.space$normalization$sd.min.q){
          temp$sd.min.q <- as.character(sd.min.q)
          df.temp <- bind_rows(df.temp, temp)
        }
      }
    } else {
      for (log.n0 in parameter.space$normalization$log.n0){
        temp$log.n0 <- log.n0
        if (!str_detect(norm.method, 'std')){
          temp$sd.min.q <- NA
          df.temp <- bind_rows(df.temp, temp)
        } else {
          for (sd.min.q in parameter.space$normalization$sd.min.q){
            temp$sd.min.q <- as.character(sd.min.q)
            df.temp <- bind_rows(df.temp, temp)
          }
        }
      }
    }
  }

  return(df.temp)
}

expand.ml <- function(temp, parameter.space){
  df.temp <- tibble()
  for (ml.method in parameter.space$training$methods){
    temp$ml.method <- ml.method
    for (fs in parameter.space$training$perform.fs){
      temp$fs <- fs
      if (!fs){
        temp$fs.method <- NA
        temp$fs.cutoff <- NA
        df.temp <- bind_rows(df.temp, temp)
      } else {
        for (fs.method in parameter.space$training$fs.method){
          temp$fs.method <- fs.method
          if (fs.method == 'Wilcoxon'){
            for (fs.cutoff in parameter.space$training$fs.cutoff.wilcoxon){
              temp$fs.cutoff <- fs.cutoff
              df.temp <- bind_rows(df.temp, temp)
            }
          } else {
            for (fs.cutoff in parameter.space$training$fs.cutoff){
              temp$fs.cutoff <- fs.cutoff
              df.temp <- bind_rows(df.temp, temp)
            }
          }
        }
      }
    }
  }
  return(df.temp)
}

# function to prepare submission scripts
prep.submission.script <- function(dataset, type, case){
  
  script.file <- here("parameter_space", "submission_temp",
                      paste0('script_', dataset, '_', type, '_', case, '.sh'))
  
  lines <- c(
    "#!/bin/bash", "#SBATCH -A zeller",
    "#SBATCH -t 3-00:00:00", "#SBATCH --mem 3G",
    ifelse(type%in%c('eggNOG', 'humann2'), 
           '#SBATCH --array=1-500', 
           '#SBATCH --array=1-200'),
    paste0("#SBATCH -o /scratch/jawirbel/siamcat/job.",
           dataset, '_', type, '_', case, ".%A_%a.out"),
    paste0("#SBATCH -e /scratch/jawirbel/siamcat/job.",
           dataset, '_', type, '_', case, ".%A_%a.err"),
    "", "module load R/3.5.1-foss-2017b-X11-20171023", "",
    paste0("Rscript ../src/3_train_models.R --dataset ", dataset, " --type ",
           type, " --case ", case, ' --subset ${SLURM_ARRAY_TASK_ID}'))
  
  tmp <- file(script.file)
  writeLines(lines, con=tmp)
  close(tmp)
}

# ##############################################################################
# prepare taxonomic jobs parameter space
parameter.space <- yaml.load_file(here("parameter_space", "job_info",
                                       'parameter_space.yaml'))

temp <- tibble(test=NA_real_)
df.temp <- expand.filtering(temp, parameter.space)
df.temp2 <- tibble()
for (j in seq_len(nrow(df.temp))){
  df.temp2 <- bind_rows(df.temp2,
                        expand.normalization(df.temp[j,], parameter.space))
}
df.temp3 <- tibble()
for (j in seq_len(nrow(df.temp2))){
  df.temp3 <- bind_rows(df.temp3,
                        expand.ml(df.temp2[j,], parameter.space))
}
df.jobs <- df.temp3 %>% 
  select(-test)

# add job ids
df.jobs <- df.jobs %>%
  mutate(job.id=sample(nrow(.)))
  
write_tsv(df.jobs, here("parameter_space", 'job_info', 'job_info.tsv'))

# ##############################################################################
# prepare functional jobs parameter space
parameter.space <- yaml.load_file(here("parameter_space", "job_info",
                                       'parameter_space_func.yml'))
temp <- tibble(test=NA_real_)
df.temp <- expand.filtering(temp, parameter.space)
df.temp2 <- tibble()
for (j in seq_len(nrow(df.temp))){
  df.temp2 <- bind_rows(df.temp2,
                        expand.normalization(df.temp[j,], parameter.space))
}
df.temp3 <- tibble()
for (j in seq_len(nrow(df.temp2))){
  df.temp3 <- bind_rows(df.temp3,
                        expand.ml(df.temp2[j,], parameter.space))
}
df.jobs <- df.temp3 %>% 
  select(-test)

# add job ids
df.jobs <- df.jobs %>%
  mutate(job.id=sample(nrow(.)))

write_tsv(df.jobs, here("parameter_space", 'job_info', 'job_info_func.tsv'))

# ##############################################################################
# prepare submission scripts

# load info about datasets
all.tasks <- read_tsv(here('parameter_space', 'data_info', 'all_tasks.tsv'))

for (i in seq_len(nrow(all.tasks))){
  prep.submission.script(all.tasks$dataset.id[i],
                         all.tasks$type[i],
                         all.tasks$case[i])
}
