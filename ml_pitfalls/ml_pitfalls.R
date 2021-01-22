# ##############################################################################
#
## ML pitfalls, demonstrated with IBD and CRC studies
##  See also the new vignette at 
##  https://siamcat.embl.de/articles/SIAMCAT_ml_pitfalls.html
#
# ##############################################################################

library("tidyverse")
library("SIAMCAT")
library("pROC")
library("ggembl")
library("here")
library("curatedMetagenomicData")

# ##############################################################################
# Supervised feature selection with CRC studies as example
fs.cutoff <- c(20, 100, 250)

# use Thomas et al and Zeller et al as 
# training and external test dataset, respectively

## features
x <- 'ThomasAM_2018a.metaphlan_bugs_list.stool'
feat.t <- curatedMetagenomicData(x=x, dryrun=FALSE)
feat.t <- feat.t[[x]]@assayData$exprs
# clean up metaphlan profiles to contain only species-level abundances
feat.t <- feat.t[grep(x=rownames(feat.t), pattern='s__'),]
feat.t <- feat.t[grep(x=rownames(feat.t),pattern='t__', invert = TRUE),]
stopifnot(all(colSums(feat.t) != 0))
feat.t <- t(t(feat.t)/100)

x <- 'ZellerG_2014.metaphlan_bugs_list.stool'
feat.z <- curatedMetagenomicData(x=x, dryrun=FALSE)
feat.z <- feat.z[[x]]@assayData$exprs
# clean up metaphlan profiles to contain only species-level abundances
feat.z <- feat.z[grep(x=rownames(feat.z), pattern='s__'),]
feat.z <- feat.z[grep(x=rownames(feat.z),pattern='t__', invert = TRUE),]
stopifnot(all(colSums(feat.z) != 0))
feat.z <- t(t(feat.z)/100)

## metadata
meta.t <- combined_metadata %>% 
  filter(dataset_name == 'ThomasAM_2018a') %>% 
  filter(study_condition %in% c('control', 'CRC'))
rownames(meta.t) <- meta.t$sampleID
meta.z <- combined_metadata %>% 
  filter(dataset_name == 'ZellerG_2014') %>% 
  filter(study_condition %in% c('control', 'CRC'))
rownames(meta.z) <- meta.z$sampleID

## fix Metaphlan profiles
species.union <- union(rownames(feat.t), rownames(feat.z))
# add Zeller_2014-only species to the Thomas_2018 matrix
add.species <- setdiff(species.union, rownames(feat.t))
feat.t <- rbind(feat.t, 
                matrix(0, nrow=length(add.species), ncol=ncol(feat.t),
                       dimnames = list(add.species, colnames(feat.t))))

# add Thomas_2018-only species to the Zeller_2014 matrix
add.species <- setdiff(species.union, rownames(feat.z))
feat.z <- rbind(feat.z, 
                matrix(0, nrow=length(add.species), ncol=ncol(feat.z),
                       dimnames = list(add.species, colnames(feat.z))))

# prepare tibble for the results
auroc.all <- tibble(cutoff=character(0), type=character(0), 
                    study.test=character(0), AUC=double(0))

## Full model without feature selection (baseline)
sc.obj.t <- siamcat(feat=feat.t, meta=meta.t,
                    label='study_condition', case='CRC')
sc.obj.t <- filter.features(sc.obj.t, filter.method = 'prevalence',
                            cutoff = 0.01)
sc.obj.t <- normalize.features(sc.obj.t,
                               norm.method = 'log.std',
                               norm.param=list(log.n0=1e-05,
                                               sd.min.q=0))
sc.obj.t <- create.data.split(sc.obj.t,
                              num.folds = 10, num.resample = 10)
sc.obj.t <- train.model(sc.obj.t, method='lasso')
sc.obj.t <- make.predictions(sc.obj.t)
sc.obj.t <- evaluate.predictions(sc.obj.t)

auroc.all <- auroc.all %>% 
  add_row(cutoff='full', type='correct', 
          study.test='Thomas_2018', 
          AUC=as.numeric(sc.obj.t@eval_data$auroc)) %>% 
  add_row(cutoff='full', type='incorrect', 
          study.test='Thomas_2018', 
          AUC=as.numeric(sc.obj.t@eval_data$auroc)) 
# external application
sc.obj.z <- siamcat(feat=feat.z, meta=meta.z,
                    label='study_condition', case='CRC')
sc.obj.z <- make.predictions(sc.obj.t, sc.obj.z)
sc.obj.z <- evaluate.predictions(sc.obj.z)
auroc.all <- auroc.all %>% 
  add_row(cutoff='full', type='correct', 
          study.test='Zeller_2014', 
          AUC=as.numeric(sc.obj.z@eval_data$auroc)) %>% 
  add_row(cutoff='full', type='incorrect', 
          study.test='Zeller_2014', 
          AUC=as.numeric(sc.obj.z@eval_data$auroc)) 

## extract p values for feature selection
sc.obj.t <- check.associations(sc.obj.t, detect.lim = 1e-05,
                               fn.plot = 'assoc_plot.pdf')
mat.assoc <- associations(sc.obj.t)
mat.assoc$species <- rownames(mat.assoc)
# sort by p-value
mat.assoc <- mat.assoc %>% as_tibble() %>% arrange(p.val)


## train models with different feature selection cutoffs and apply
## INCORRECT PROCEDURE
for (x in fs.cutoff){
  # select x number of features based on p-value ranking
  feat.train.red <- feat.t[mat.assoc %>%
                             slice(seq_len(x)) %>%
                             pull(species),]
  sc.obj.t.fs <- siamcat(feat=feat.train.red, meta=meta.t,
                         label='study_condition', case='CRC')
  # normalize the features without filtering
  sc.obj.t.fs <- normalize.features(sc.obj.t.fs,
                                    norm.method = 'log.std',
                                    norm.param=list(log.n0=1e-05,
                                                    sd.min.q=0),
                                    feature.type = 'original')
  # take the same cross validation split as before
  data_split(sc.obj.t.fs) <- data_split(sc.obj.t)
  # train
  sc.obj.t.fs <- train.model(sc.obj.t.fs, method = 'lasso')
  # make predictions
  sc.obj.t.fs <- make.predictions(sc.obj.t.fs)
  # evaluate predictions and record the result
  sc.obj.t.fs <- evaluate.predictions(sc.obj.t.fs)
  auroc.all <- auroc.all %>% 
    add_row(cutoff=as.character(x), type='incorrect', 
            study.test='Thomas_2018',
            AUC=as.numeric(sc.obj.t.fs@eval_data$auroc))
  # apply to the external dataset and record the result
  sc.obj.z <- siamcat(feat=feat.z, meta=meta.z,
                      label='study_condition', case='CRC')
  sc.obj.z <- make.predictions(sc.obj.t.fs, sc.obj.z)
  sc.obj.z <- evaluate.predictions(sc.obj.z)
  auroc.all <- auroc.all %>% 
    add_row(cutoff=as.character(x), type='incorrect', 
            study.test='Zeller_2014', 
            AUC=as.numeric(sc.obj.z@eval_data$auroc))
}

## train models with different feature selection cutoffs and apply
## CORRECT PROCEDURE
for (x in fs.cutoff){
  # train using the original SIAMCAT object 
  # with correct version of feature selection
  sc.obj.t.fs <- train.model(sc.obj.t, method = 'lasso', 
                             perform.fs = TRUE,
                             param.fs = list(thres.fs = x,
                                             method.fs = "AUC",
                                             direction='absolute'))
  # make predictions
  sc.obj.t.fs <- make.predictions(sc.obj.t.fs)
  # evaluate predictions and record the result
  sc.obj.t.fs <- evaluate.predictions(sc.obj.t.fs)
  auroc.all <- auroc.all %>% 
    add_row(cutoff=as.character(x), type='correct', 
            study.test='Thomas_2018',
            AUC=as.numeric(sc.obj.t.fs@eval_data$auroc))
  # apply to the external dataset and record the result
  sc.obj.z <- siamcat(feat=feat.z, meta=meta.z,
                      label='study_condition', case='CRC')
  sc.obj.z <- make.predictions(sc.obj.t.fs, sc.obj.z)
  sc.obj.z <- evaluate.predictions(sc.obj.z)
  auroc.all <- auroc.all %>% 
    add_row(cutoff=as.character(x), type='correct', 
            study.test='Zeller_2014', 
            AUC=as.numeric(sc.obj.z@eval_data$auroc))
}

## Plot the results
g <- auroc.all %>%
  # facetting for plotting
  filter(study.test%in%c('Thomas_2018', 'Zeller_2014')) %>%
  mutate(split=case_when(study.test=="Thomas_2018"~
                           'Cross validation (Thomas 2018)',
                         TRUE~"External validation (Zeller 2014)")) %>%
  # convert to factor to enforce ordering
  mutate(cutoff=factor(cutoff, levels = c(fs.cutoff, 'full'))) %>%
  mutate(type=case_when(type=='full'~'correct', TRUE~type)) %>%
  ggplot(aes(x=cutoff, y=AUC, col=type)) +
    geom_point() + geom_line(aes(group=type)) +
    facet_grid(~split) +
    coord_cartesian(ylim=c(0.5, 1), expand = TRUE) +
    theme_publication(panel.grid = 'major_y') +
    xlab('Features selected') +
    ylab('AUROC') +
    scale_colour_manual(values = c('correct'='blue', 'incorrect'='red'), 
                        name='') +
    scale_y_continuous(limits = c(0.5, 1), expand = c(0,0)) + 
    theme(legend.position = 'top')

ggsave(g, filename = here('figures', 'ml_pitfalls', 'feature_selection.pdf'),
       width = 85, height = 50, units = 'mm', useDingbats=FALSE)

# ##############################################################################
# Naive dependent sample splitting with IBD studies as example

# 
data.location <- 'https://www.embl.de/download/zeller/'
# metadata
meta.all <- read_tsv(paste0(data.location, 'CD_meta/meta_all.tsv'))
# features
feat.motus <- read.table(paste0(data.location, 'CD_meta/feat_rel_filt.tsv'),
                         sep='\t', stringsAsFactors = FALSE,
                         check.names = FALSE)


datasets <- c('metaHIT', 'Lewis_2015', 'He_2017', 'Franzosa_2019', 'HMP2')

# load data
meta.all <- read_tsv(here('ibd_meta_analysis', 'data', 'meta_all.tsv'))
meta.ind <- read_tsv(here('ibd_meta_analysis', 'data', 'meta_individual.tsv'))
feat.motus <- read.table(here('ibd_meta_analysis', 'data',
                              'feat_all_rel.tsv'),
                         sep='\t', stringsAsFactors = FALSE,
                         check.names = FALSE)

# ##############################################################################
# train HMP2 with/without blocking
dataset.block <- 'HMP2'

meta.temp <- meta.all %>%
  filter(Sample_ID %in% meta.ind$Sample_ID | Study==dataset.block)
pred.block <- matrix(NA, nrow=nrow(meta.temp), ncol=2,
                     dimnames = list(meta.temp$Sample_ID,
                                     c('naive', 'blocked')))

if (!file.exists(here("ml_pitfalls", 'predictions',
                      paste0("predictions_blocking_",
                             dataset.block, ".tsv")))){
  # train
  meta.train <- meta.all %>%
    filter(Study==dataset.block) %>%
    group_by(Individual_ID) %>%
    sample_n(5, replace = TRUE) %>%
    distinct() %>%
    as.data.frame()
    

  rownames(meta.train) <- meta.train$Sample_ID
  feat.train <- feat.motus[,meta.train$Sample_ID]
  sc.obj.train <- siamcat(feat=feat.train, meta=meta.train,
                          label='Group', case='CD')
  sc.obj.train <- filter.features(sc.obj.train)
  sc.obj.train <- normalize.features(sc.obj.train, norm.method = 'log.std',
                                     norm.param=list(log.n0=1e-05,
                                                     sd.min.q=1))
  sc.obj.train.naive <- create.data.split(sc.obj.train,
                                          num.folds = 10, num.resample = 10)
  sc.obj.train.blocked <- create.data.split(sc.obj.train,
                                            num.folds = 10,
                                            num.resample = 10,
                                            inseparable = 'Individual_ID')
  sc.obj.train.naive <- train.model(sc.obj.train.naive, method='lasso')
  sc.obj.train.blocked <- train.model(sc.obj.train.blocked, method='lasso')
  save(sc.obj.train.naive, sc.obj.train.blocked,
       file = here("ml_pitfalls",
                   paste0("models_naive_vs_blocked_",
                          dataset.block, ".RData")))
  # apply
  for (i in datasets){
    if (i == dataset.block){
      
      sc.obj.train.naive <- make.predictions(sc.obj.train.naive)
      sc.obj.train.blocked <- make.predictions(sc.obj.train.blocked)

      temp <- pred_matrix(sc.obj.train.naive)
      pred.block[rownames(temp), 'naive'] <- rowMeans(temp)
      temp <- pred_matrix(sc.obj.train.blocked)
      pred.block[rownames(temp), 'blocked'] <- rowMeans(temp)
    } else {
      meta.test <- meta.all %>%
        filter(Study==i) %>%
        filter(Sample_ID %in% rownames(pred.block)) %>%
        as.data.frame()
      rownames(meta.test) <- meta.test$Sample_ID
      feat.test <- feat.motus[,meta.test$Sample_ID]
      sc.obj.test <- siamcat(feat=feat.test, meta=meta.test,
                             label='Group', case='CD')
      sc.obj.test <- make.predictions(sc.obj.train.naive, sc.obj.test)
      temp <- pred_matrix(sc.obj.test)
      pred.block[rownames(temp), 'naive'] <- rowMeans(temp)
      sc.obj.test <- make.predictions(sc.obj.train.blocked, sc.obj.test)
      temp <- pred_matrix(sc.obj.test)
      pred.block[rownames(temp), 'blocked'] <- rowMeans(temp)
    }
  }

  pred.block <- data.frame(pred.block)
  pred.block$Sample_ID <- rownames(pred.block)
  pred.block <- as_tibble(pred.block)
  write_tsv(pred.block,
            path = here('ml_pitfalls',
                        paste0('predictions_blocking_',
                               dataset.block, '.tsv')))
} else {
  pred.block <- read_tsv(here('ml_pitfalls',
                              paste0('predictions_blocking_',
                                     dataset.block, '.tsv')))
}

# evaluate
auroc.all <- tibble()
for (study.test in datasets){
    
  df.temp <- inner_join(pred.block, meta.temp) %>%
    filter(Study==study.test)
  for (type in c('naive', 'blocked')){
    temp <- roc(predictor=df.temp[[type]], response = df.temp$Group, ci=TRUE,
                levels = c('CD', 'CTR'), direction = '>')
    auroc.all <- bind_rows(auroc.all,
                           tibble(type=type,
                                  study.test=study.test,
                                  AUC=c(temp$auc)))
  }
}

# plot
g <- auroc.all %>%
  mutate(type=factor(type, levels = c('naive', 'blocked'))) %>%
  mutate(CV=case_when(study.test=='HMP2'~'CV', TRUE~'Holdout')) %>%
  ggplot(aes(x=study.test, y=AUC, fill=type)) +
  geom_bar(stat='identity', position = position_dodge(), col='black') +
  theme_publication() +
  coord_cartesian(ylim=c(0.5, 1)) +
  scale_fill_manual(values = c('blocked'='blue', 'naive'='red'), 
                      name='') +
  facet_grid(~CV, space = 'free', scales = 'free') +
  xlab('') + ylab('AUROC') +
  theme(legend.position = c(0.8, 0.8))

ggsave(g, filename = here('figures', 'ml_pitfalls', 'blocked_CV.pdf'),
       width = 85, height = 50, units = 'mm', useDingbats=FALSE)
