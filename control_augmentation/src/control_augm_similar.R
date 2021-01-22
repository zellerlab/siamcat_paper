# ##############################################################################
#
## Train a control-augmented model for datasets supplemented with similar data
#
# ##############################################################################

library("tidyverse")
library("optparse")
library("SIAMCAT")
library("pROC")
library("progress")
library("here")

source(here('utils', 'train_ctr_augm.R'))
source(here('utils', 'cross_predictions.R'))

for (ml.method in c('lasso', 'enet-0.5', 'enet', 'randomForest')){
  fn.cross.pred <- here('control_augmentation', 'files',
                        paste0('cross_prediction_', ml.method, 
                               '_similar.tsv'))
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
  
  ml <- ifelse(ml.method=='enet-0.5', 'enet', ml.method)
  fn.results <- here('parameter_space', 'job_info', 'full_results.tsv')
  stopifnot(file.exists(fn.results))
  motu.results <- read_tsv(fn.results) %>%
    filter(type=='mOTUs2')
  print(dim(motu.results))
  motu.tasks <- motu.results %>%
    select(dataset.id, case, type) %>%
    distinct()
  print(dim(motu.tasks))
  best.jobs <- motu.results %>%
    filter(ml.method==ml) %>%
    group_by(job.id) %>%
    summarise(n=sum(!is.na(auroc)), auc=mean(auroc)) %>%
    filter(n==nrow(motu.tasks)) %>%
    arrange(desc(auc))
  print(dim(best.jobs))
  info <- motu.results %>%
    filter(job.id==best.jobs$job.id[1]) %>%
    select(1:10) %>%
    distinct()
  info$norm.method <- 'log.std'
  info$sd.min.q <- 0
  print(info)
  
  # choose datasets to be augmented here
  augm.datasets <- c('Yu_2017', 'Jie_2017', 'He_2017', 'Qin_2012')
  for (d in augm.datasets){
    message("train augmented ", d)
    fn.results <- here('control_augmentation', 'models', 'similar',
                       paste0('sc_', d, '_', ml.method, 
                              '_augmented.RData'))
    if (!file.exists(fn.results)){
      
      # run the script same as before
      fn.sc <- list.files(here('parameter_space', 'sc'), 
                          pattern=paste0(d, '.+mOTUs2'),
                          full.names = TRUE)
      load(fn.sc)
      # filtering
      if (str_detect(info$filt.method, ';')){
        methods <- str_split(info$filt.method, ';')[[1]]
        cutoffs <- str_split(info$filt.cutoff, ';')[[1]]
        sc.obj <- filter.features(sc.obj,
                                  filter.method = methods[1],
                                  cutoff = as.numeric(cutoffs[1]),
                                  verbose=0)
        sc.obj <- filter.features(sc.obj,
                                  filter.method=methods[2],
                                  cutoff = as.numeric(cutoffs[2]),
                                  feature.type = 'filtered', verbose=0)
      } else {
        sc.obj <- filter.features(sc.obj,
                                  filter.method = info$filt.method,
                                  cutoff = as.numeric(info$filt.cutoff),
                                  verbose=0)
      }
      if (info$norm.method=='none'){
        sc.obj <- normalize.features(sc.obj, norm.method = 'pass')
      } else {
        # normalization
        sc.obj <- normalize.features(
          sc.obj,
          norm.method = info$norm.method,
          norm.param = list(log.n0=as.numeric(info$log.n0),
                            sd.min.q=as.numeric(info$sd.min.q)),verbose=0)
      }
      
      sc.obj.trained <- sc.obj
      # get control datasets
      correct.list <- list()
      for (d2 in setdiff(augm.datasets, d)){
        fn.sc <- list.files(here('parameter_space', 'sc'), pattern=d2) %>% 
          grep(pattern='mOTUs', value = TRUE) %>% sample(1)
        load(here('parameter_space', 'sc', fn.sc))
        meta <- meta(sc.obj)
        meta <- meta[meta$Group=='CTR',]
        feat <- get.orig_feat.matrix(sc.obj)
        feat <- feat[,rownames(meta)]
        correct.list[[(length(correct.list) + 1)]] <- feat
      }
      # model training
      sc.obj <- train.model.ctr(sc.obj.trained,
                                correct.list = correct.list,
                                method=ml.method,
                                perform.fs = info$fs,
                                param.fs = list(thres.fs=info$fs.cutoff,
                                                method.fs=info$fs.method,
                                                direction='absolute'),
                                n=2,
                                verbose=1)
      sc.obj <- make.predictions(sc.obj,verbose=0)
      sc.obj <- evaluate.predictions(sc.obj, verbose=0)
      
      # save model
      save(sc.obj, file=fn.results)
    } else {
      load(fn.results)
    }
    
    # calculate cross-predictions directly
    df.cross.predictions <- bind_rows(
      df.cross.predictions, 
      f.cross.prediction(sc.obj, d, names(sc.obj@label$info)[2], 
                         setdiff(augm.datasets, d)))
    
  }
  df.cross.predictions <- df.cross.predictions %>% 
    distinct()
  write_tsv(df.cross.predictions, 
            path = fn.cross.pred)
  
}



