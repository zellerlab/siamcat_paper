# ##############################################################################
#
## Check transfer to other datasets for augmented and normal best models
#
# ##############################################################################

library("tidyverse")
library("SIAMCAT")
library("pROC")
library("here")
library("ggembl")

source(here('utils', 'cross_predictions.R'))

# ##############################################################################
# load data sets
motu.tasks <- read_tsv(here('parameter_space', 'files',
                              'auroc_all.tsv')) %>%
  filter(type=='mOTUs2')

for (ml.method in c('randomForest', 'enet', 'enet-0.5', 'lasso')){
  for (ctr_type in c('cohort_2', 'cohort_5', 'other_ctr', 'random')){
    fn.cross.pred <- here('control_augmentation', 'files',
                          paste0('cross_prediction_', ml.method,
                                 '_', ctr_type, '.tsv'))
    message("#------------------------------------------------------\n",
            ml.method, '-', ctr_type,
            "\n#------------------------------------------------------")
    if (!file.exists(fn.cross.pred)){
      # prepare cross-prediction matrix
      df.cross.predictions <- tibble(n.group=integer(0),
                                     n.pos=integer(0),
                                     frac=double(0),
                                     cp.auroc=double(0),
                                     dataset.id=character(0),
                                     Group=character(0),
                                     dataset.test=character(0),
                                     dataset.train=character(0),
                                     dataset.id.train=character(0),
                                     case.train=character(0),
                                     auroc=double(0),
                                     ext.auroc=double(0))
      for (i in seq_len(nrow(motu.tasks))){
        dataset.id <- motu.tasks$dataset.id[[i]]
        case <- motu.tasks$case[[i]]
        message(dataset.id, '-', case)
        fn.sc <- list.files(here('control_augmentation', 'models', ctr_type),
                      pattern=paste0(dataset.id, '-mOTUs2-', case),
                      full.names = TRUE) %>%
          grep(pattern=ml.method, value=TRUE)
        if (ml.method=='enet'){
          fn.sc <- grep(fn.sc, pattern = 'enet-0.5', value=TRUE, invert = TRUE)
        }
        if (length(fn.sc)==0){
          message('Problem for ', dataset.id, '-', case)
          next()
        }
        load(fn.sc)
        if (ctr_type=='random'){
          df.cross.predictions <- bind_rows(
            df.cross.predictions,
            f.cross.prediction(sc.obj, dataset.id, case, datasets))
        } else {
          df.cross.predictions <- bind_rows(
            df.cross.predictions,
            f.cross.prediction(sc.obj, dataset.id, case, ctr_type))
        }

      }
      df.cross.predictions <- df.cross.predictions %>%
        distinct()

      write_tsv(df.cross.predictions,
                file = fn.cross.pred)
    }
  }
}
