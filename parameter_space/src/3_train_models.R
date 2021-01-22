# ##############################################################################
#
## Script to run the parameter space exploration
#
# ##############################################################################

.libPaths(c('/g/scb2/zeller/SHARED/software/R/3.5', .libPaths()))

library("tidyverse")
library("optparse")
library("SIAMCAT")
library("here")

option_list <- list(
  make_option('--dataset', type='character', 
              help='ID of the dataset'),
  make_option('--type', type='character',
              help='Which kind of features should be used'),
  make_option('--case', type='character',
              help='What is the positive case?'),
  make_option('--subset', type='integer',
              help='Which subset of jobs is to be used?')
)
opt <- list(dataset="Feng_2015",type="mOTUs2", case="ADA", subset=1)
opt <- parse_args(OptionParser(option_list=option_list))

job.file <- ifelse(opt$type %in% c('humann2', 'eggNOG'), 
                   'job_info_func.tsv', 
                   'job_info.tsv')

job.info <- read_tsv(
  here('parameter_space', 'job_info', job.file),
  col_types = cols(
    filt.method = col_character(),
    filt.cutoff = col_character(),
    norm.method = col_character(),
    log.n0 = col_double(),
    sd.min.q = col_double(),
    ml.method = col_character(),
    fs = col_logical(),
    fs.method = col_character(),
    fs.cutoff = col_double(),
    job.id = col_double()
  )) %>% 
  arrange(job.id)

max.jobs <- ifelse(opt$type%in%c('eggNOG', 'humann2'), 500, 200)
job.info$set <- rep(seq_len(max.jobs), 
                    each=ceiling(nrow(job.info)/max.jobs))[
                      seq_len(nrow(job.info))]
job.info <- job.info %>% 
  filter(set==opt$subset)
job.info$set <- NULL

tag <- paste(opt, collapse = '-')

# check if everything has been done or not
fn.results <- here("parameter_space/results", 
                   paste0('results_', tag, '.tsv'))
if (file.exists(fn.results)){
  stop('Every job for this dataset has run already!')
}

# check which jobs have been done already
fn.partial <- here("parameter_space/temp", 
                   paste0('partial_results_', tag, '.tsv'))
if (!file.exists(fn.partial)){
  partial.results <- job.info %>% 
    mutate(dataset.id=opt$dataset, 
           type=opt$type, 
           case=opt$case) %>% 
    mutate(processed=FALSE, auroc=NA_real_, time=NA_real_,
           ci.low=NA_real_, ci.high=NA_real_)
  write(paste0(colnames(partial.results), collapse = "\t"), 
        file = fn.partial)
} else {
  partial.results <- read_tsv(fn.partial,
                             col_types = cols(
                               dataset.id = col_character(),
                               type = col_character(),
                               case = col_character(),
                               processed = col_logical(),
                               auroc = col_double(),
                               time = col_double(),
                               ci.low = col_double(),
                               ci.high = col_double(),
                               filt.method = col_character(),
                               filt.cutoff = col_character(),
                               norm.method = col_character(),
                               log.n0 = col_double(),
                               sd.min.q = col_double(),
                               ml.method = col_character(),
                               fs = col_logical(),
                               fs.method = col_character(),
                               fs.cutoff = col_double(),
                               job.id = col_double()
                             ))
}

fn.sc <- here('parameter_space', 'sc', paste0(
  paste(c('sc', opt$dataset, opt$case, opt$type), collapse = '_'),
  '.RData'))
print(fn.sc)
load(fn.sc)

while (sum(!partial.results$processed) > 0){
  
  i <- partial.results %>% 
    filter(!processed) %>% 
    slice(1) %>% 
    pull(job.id)
  
  info <- job.info %>% filter(job.id==i)
  sc.job.current <- sc.obj
  x.start <- Sys.time()
  # more fancy with tryCatch, since some jobs can fail, which is okay
  auc.list <- tryCatch({
    # filtering
    if (str_detect(info$filt.method, ';')){
      methods <- str_split(info$filt.method, ';')[[1]]
      cutoffs <- str_split(info$filt.cutoff, ';')[[1]]
      sc.job.current <- filter.features(sc.job.current,
                                        filter.method = methods[1],
                                        cutoff = as.numeric(cutoffs[1]), 
                                        verbose=0)
      sc.job.current <- filter.features(sc.job.current,
                                        filter.method=methods[2],
                                        cutoff = as.numeric(cutoffs[2]),
                                        feature.type = 'filtered', verbose=0)
    } else {
      sc.job.current <- filter.features(sc.job.current, 
                                        filter.method = info$filt.method,
                                        cutoff = as.numeric(info$filt.cutoff),
                                        verbose=0)
    }
    if (info$norm.method=='none'){
      sc.job.current <- normalize.features(sc.job.current,
                                           norm.method = 'pass')
    } else {
      # normalization
      sc.job.current <- 
        normalize.features(
          sc.job.current,
          norm.method = info$norm.method,
          norm.param = list(log.n0=as.numeric(info$log.n0),
                            sd.min.q=as.numeric(info$sd.min.q)),
          verbose=0)
    }
    
    # model training
    sc.job.current <- train.model(sc.job.current,
                                  method=info$ml.method,
                                  perform.fs = info$fs, 
                                  param.fs = list(thres.fs=info$fs.cutoff, 
                                                  method.fs=info$fs.method, 
                                                  direction='absolute'),
                                  verbose=0)
    sc.job.current <- make.predictions(sc.job.current,verbose=0)
    sc.job.current <- evaluate.predictions(sc.job.current, verbose=0)
    
    temp <- pROC::roc(predictor=rowMeans(pred_matrix(sc.job.current)), 
                      response=sc.job.current@label$label,
                      direction = '<',
                      levels = sc.job.current@label$info,
                      ci=TRUE)
    x.end <- Sys.time()
    temp <- list(auroc=temp$auc, ci.high=temp$ci[3], ci.low=temp$ci[1],
                 time=as.numeric(difftime(x.end, x.start, units = 'secs')))},
    error=function(cond){return(list(auroc=NA_real_, 
                                     ci.low=NA_real_, 
                                     ci.high=NA_real_,
                                     time=-1))})
  
  partial.results$processed[which(partial.results$job.id==i)] <- TRUE
  partial.results$auroc[which(partial.results$job.id==i)] <- 
    as.numeric(auc.list$auroc)
  partial.results$ci.low[which(partial.results$job.id==i)] <- auc.list$ci.low
  partial.results$ci.high[which(partial.results$job.id==i)] <- auc.list$ci.high
  partial.results$time[which(partial.results$job.id==i)] <- auc.list$time
    
  
  write_tsv(partial.results %>% 
              filter(job.id==i), 
            path = fn.partial, append = TRUE)
  
  message('Finished job #', i)
}

system(paste0('mv ', fn.partial, ' ', fn.results))
