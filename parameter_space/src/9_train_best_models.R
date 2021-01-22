# ##############################################################################
#
## Re-train best parameter set models
#
# ##############################################################################

.libPaths(c('/g/scb2/zeller/SHARED/software/R/3.5', .libPaths()))

library("tidyverse")
library("SIAMCAT")
library("optparse")
library("here")

source(here('parameter_space', 'src', 'utils.R'))

option_list <- list(
  make_option('--dataset', type='character',
              help='ID of the dataset'),
  make_option('--type', type='character',
              help='Which kind of features should be used'),
  make_option('--case', type='character',
              help='What is the positive case?'),
  make_option('--ml_method', type='character',
              help='Which type of machine learning model should be used?')
)
opt <- list(dataset='Bedarf_2017', type='mOTUs2', case='PD', ml_method='lasso')
opt <- parse_args(OptionParser(option_list=option_list))
print(opt)
stopifnot(opt$ml.method %in% c('enet', 'lasso', 'enet-0.5', 'randomForest'))
ml <- ifelse(opt$ml_method=='enet-0.5', 'enet', opt$ml_method)
stopifnot(opt$type=='mOTUs2')

# ##############################################################################
# load data sets
motu.results <- read_tsv(here('parameter_space', 
                              'job_info', 
                              'full_results.tsv')) %>% 
  filter(type=='mOTUs2')
motu.tasks <- motu.results %>% 
  select(dataset.id, case, type) %>% 
  distinct()
print(dim(motu.tasks))
# ##############################################################################
# best parameter set 
best.jobs <- motu.results %>% 
  group_by(job.id) %>% 
  summarise(n=sum(!is.na(auroc)), auc=mean(auroc)) %>% 
  filter(n==nrow(motu.tasks)) %>% 
  arrange(desc(auc))

info <- motu.results %>% 
  filter(ml.method==ml)
best.jobs.red <- best.jobs %>% 
  filter(job.id %in% info$job.id)
info <- info %>% 
  filter(job.id==best.jobs.red$job.id[1]) %>% 
  select(1:10) %>% 
  distinct()
# hard fix for other ML methods.... should be adjusted at some point, maybe?
info$norm.method <- 'log.std'
info$sd.min.q <- 0
print(info)


fn.model.trained <- here('parameter_space', 'models',
                         paste0('sc_trained_', opt$dataset, 
                                '_', opt$case, '_', 
                                ml, '.RData'))

if (!file.exists(fn.model.trained)){
  load(here('parameter_space', 'sc', 
            paste0('sc_', opt$dataset, '_', 
                   opt$case, '_', opt$type, '.RData')))
  # preprocess
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
  # normalize
  sc.obj <- 
    normalize.features(sc.obj,
                       norm.method = info$norm.method,
                       norm.param = list(log.n0=as.numeric(info$log.n0),
                                         sd.min.q=as.numeric(
                                           info$sd.min.q)),
                       verbose=0)
  # train model
  if (ml=='enet-0.5'){
    sc.obj <- train.model(sc.obj,
                          method="enet",
                          perform.fs = info$fs, 
                          param.fs = list(thres.fs=info$fs.cutoff, 
                                          method.fs=info$fs.method, 
                                          direction='absolute'),
                          param.set = list(alpha=0.5),
                          verbose=1)
  } else {
    sc.obj <- train.model(sc.obj,
                          method=opt$ml_method,
                          perform.fs = info$fs, 
                          param.fs = list(thres.fs=info$fs.cutoff, 
                                          method.fs=info$fs.method, 
                                          direction='absolute'),
                          verbose=1)
  }
  
  # apply and evaluate cross-prediction
  # determine threshold
  sc.obj <- make.predictions(sc.obj, verbose=0)
  sc.obj <- evaluate.predictions(sc.obj, verbose = 0)
  sc.obj.train <- sc.obj
  save(sc.obj.train, file = fn.model.trained)
} 

