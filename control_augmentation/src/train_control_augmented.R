# ##############################################################################
#
## Script to run the parameter space exploration
#
# ##############################################################################

.libPaths(c('/g/scb2/zeller/SHARED/software/R/3.5', .libPaths()))

library("tidyverse")
library("optparse")
library("SIAMCAT")
library("progress")
library("here")

source(here('utils', 'train_ctr_augm.R'))
source(here('utils', 'load_ctr_datasets.R'))

option_list <- list(
  make_option('--dataset', type='character',
              help='ID of the dataset'),
  make_option('--type', type='character',
              help='Which kind of features should be used'),
  make_option('--case', type='character',
              help='What is the positive case?'),
  make_option('--ctr_type', type='character',
              help='Which data should be used for control augmentation?'),
  make_option('--ml_method', type='character',
              help='Which type of machine learning model should be used?')
)
opt <- list(dataset='Bedarf_2017', type='mOTUs2', case='PD',
            ctr_type='cohort_2', ml_method='lasso')
opt <- parse_args(OptionParser(option_list=option_list))
print(opt)
stopifnot(opt$ml_method %in% c('enet', 'lasso', 'enet-0.5', 'randomForest'))
ml <- ifelse(opt$ml_method=='enet-0.5', 'enet', opt$ml_method)
stopifnot(opt$type=='mOTUs2')
stopifnot(opt$ctr_type %in% c('cohort_2', 'cohort_5', 'other_ctr', 
                              'random', 'similar'))

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
tag <- paste(opt, collapse = '-')
print(tag)

# ##############################################################################
# load
fn.sc <- here('parameter_space', 'sc', paste0(
  paste(c('sc', opt$dataset, opt$case, opt$type), collapse = '_'),
  '.RData'))
print(fn.sc)
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

# load control samples
l <- .f_load_control_data(ifelse(str_detect(opt$ctr_type, '^cohort'), 
                                 'cohort', opt$ctr_type))
correct.list <- l$feat.ctr

# model training
sc.obj <- train.model.ctr(sc.obj.trained,
                          correct.list = correct.list,
                          method=opt$ml_method,
                          perform.fs = info$fs,
                          param.fs = list(thres.fs=info$fs.cutoff,
                                          method.fs=info$fs.method,
                                          direction='absolute'),
                          n=ifelse(opt$ctr_type == 'cohort_5', 5, 2),
                          verbose=3)
sc.obj <- make.predictions(sc.obj,verbose=0)
sc.obj <- evaluate.predictions(sc.obj, verbose=0)

fn.results <- here('control_augmentation', 'models', opt$ctr_type,
                   paste0('sc_', tag, '_augmented.RData'))
if (opt$ctr_type == 'random') {
  save(sc.obj, l$datasets, file=fn.results)
} else {
  save(sc.obj, file=fn.results)
}

