# ##############################################################################
#
## Script for the reproduction of the results from  Pasolli et al
##    using SIAMCAT
#
## Usage: Rscript fit_model_pasolli.R dataset
##  with dataset being the identifier of the dataset 
#   (see below for allowed values)
#
# ##############################################################################

.libPaths(c("/g/scb2/zeller/jawirbel/software/R/3.5", .libPaths()))

library("SIAMCAT")
library("tidyverse")
library("matrixStats")
library("pROC")
library("here")

args <- commandArgs(trailingOnly = TRUE)
data.set <- args[1]

# check dataset
if (!data.set %in% c('metahit', 'Zeller_fecal_colorectal_cancer',
                     'Quin_gut_liver_cirrhosis',
                     'Chatelier_gut_obesity', 'WT2D', 't2dmeta')){
  stop('Dataset not considered')
}

# ##############################################################################
# get features
data.all <- t(read.table(here('data', 'pasolli', 'pasolli.txt.bz2'), 
                         sep='\t', stringsAsFactors = FALSE,
                         quote='', comment.char = '', row.names = 1))
data.study <- data.all %>% as_tibble() %>%
  filter(str_detect(dataset_name, pattern =data.set)) %>% t() %>% as.matrix

# metadata
meta <- data.frame(t(as.matrix(data.study[
  grep(x=rownames(data.study), pattern='k__', invert = TRUE),])),
  stringsAsFactors = FALSE)
rownames(meta) <- meta$sampleID

# features
feat.all <- as.matrix(data.study[grep(x=rownames(data.study), pattern='k__'),])
colnames(feat.all) <- meta$sampleID
feat.all.species <- feat.all[grep(x=rownames(feat.all), pattern='s__'),]
feat.all.species <- feat.all.species[grep(x=rownames(feat.all.species),
                                          pattern='t__', invert = TRUE),]

# feature wrangling
feat.rel <- prop.table(apply(feat.all.species, 2, as.numeric), 2)
rownames(feat.rel) <- rownames(feat.all.species)

# ##############################################################################
# label hard-coding
if (data.set == 'Zeller_fecal_colorectal_cancer'){
  meta <- meta[meta$disease != 'large_adenoma',]
  meta$disease[meta$disease == 'cancer'] <- 'case'
  meta$disease[meta$disease == 'small_adenoma'] <- 'n'
} else if (data.set == 'metahit'){
  meta$disease[str_detect(meta$disease, 'ibd')] <- 'case'
} else if (data.set == 'Quin_gut_liver_cirrhosis'){
  meta$disease[meta$disease == 'cirrhosis'] <- 'case'
} else if (data.set == 'Chatelier_gut_obesity'){
  meta <- meta[meta$disease != 'n',]
  meta$disease[meta$disease == 'lean'] <- 'n'
  meta$disease[meta$disease == 'obesity'] <- 'case'
} else if (data.set == 'WT2D'){
  meta <- meta[meta$disease != 'impaired_glucose_tolerance',]
  meta$disease[meta$disease == 't2d'] <- 'case'
} else if (data.set == 't2dmeta'){
  meta$disease[meta$disease == 't2d'] <- 'case'
  meta <- meta[meta$disease %in% c('n', 'case'),]
}

# ##############################################################################
# SIAMCAT
siamcat <- siamcat(feat=feat.rel, meta=meta, label='disease',
                   case='case', control='n')
siamcat <- filter.features(siamcat, filter.method = 'abundance',
                           cutoff=1e-10, verbose=3)
siamcat <- create.data.split(siamcat, num.folds=10, num.resample = 20)
siamcat <- train.model(siamcat, method='randomForest', verbose = 0,
                       feature.type = 'filtered')
siamcat <- make.predictions(siamcat, verbose=0)
siamcat <- evaluate.predictions(siamcat, verbose=0)

# evaluation predictions
rocsumm = roc(response = label(siamcat)$label,
              predictor = apply(pred_matrix(siamcat), 1, 'mean'),
              ci = TRUE)
cat(data.set)
cat(c(rocsumm$ci), '\n')
