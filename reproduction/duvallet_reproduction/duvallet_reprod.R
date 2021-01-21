# ##############################################################################
#
# Script to reproduce the results from Duvallet et al., 
#   results should be more or less perfectly matching
#
## Usage: Rscript duvallet_reprod.R --metadata_in location_of_metadata_file 
##    --label_in disease_label --feat_in location_of_feature_file
#
## Results are stored in a results file on /scratch
#
# ##############################################################################

.libPaths(c("/g/scb2/zeller/jawirbel/software/R/3.5", .libPaths()))

# load packages
library('feather')
library("SIAMCAT")
library("optparse")
library('pROC')

# define arguments
option_list <- list(
  make_option('--metadata_in', type='character', 
              help='Input file containing meta-data'),
  make_option('--label_in', type='character',
              help='String for positive labe'),
  make_option('--feat_in', type='character', 
              help='Input file containing features')
)

opt <- parse_args(OptionParser(option_list=option_list))

# print parameters of the run
cat("Reproduction of the results from Duvallet et al.\n")
cat("=== Paramaters of the run:\n")
cat('metadata_in  =', opt$metadata_in, '\n')
cat('label_in     =', toupper(opt$label_in), '\n')
cat('feat_in      =', opt$feat_in, '\n')
cat('\n')

data.tag <- strsplit(
  strsplit(opt$feat_in, split='\\.')[[1]][1], 
  split='/')[[1]][4]

start.time <- proc.time()[3]

# get files as arguments
feat.file <- opt$feat_in
meta.file <- opt$metadata_in

# ##############################################################################
# load and clean data
# load features
feat <- read_feather(feat.file)
column.names <- feat$index
feat$index <- NULL
feature.names <- colnames(feat)
feat <- t(data.frame(feat))
colnames(feat) <- column.names

# load metadata
meta <- read_feather(meta.file)
meta <- data.frame(meta)
rownames(meta) <- meta[,1]

# check if label contains H or not
healthy.state <- ifelse('H' %in% unique(meta$DiseaseState), 'H', 'nonIBD')

target <- opt$label_in
if (grepl(target, pattern='_')){
  target <- strsplit(target, split='_')[[1]]
}

label <- create.label(meta=meta, label='DiseaseState', 
                      case =target, control=healthy.state)

# ##############################################################################
# start SIAMCAT workflow
siamcat <- siamcat(feat=feat, label=label, meta=meta)

# collapse at genus level
siamcat <- summarize.features(siamcat, level='g__')

# prepare cross-validation 
siamcat <- create.data.split(siamcat, num.folds = 5, 
                             num.resample = 5, 
                             stratify = TRUE)

# train random forest classifier
siamcat <- train.model(siamcat, method='randomForest', 
                       param.set=list('ntree'=c(1000, 1000)), 
                       feature.type = 'original')

# make predictions
siamcat <- make.predictions(siamcat)

# evaluation predictions
rocsumm = roc(response = label(siamcat)$label, 
              predictor = apply(pred_matrix(siamcat), 1, 'mean'),
              ci = TRUE)

cat('\n', feat.file, '\n')
cat(c(rocsumm$ci), '\n')

save(siamcat, file=paste0('/scratch/jawirbel/duvallet_reprod/files/', 
                          data.tag, '_', opt$label_in, '.RData'))

line <- paste(c(paste0(data.tag, '_', opt$label_in), rocsumm$ci), 
              collapse='\t')
write(line,file="/scratch/jawirbel/duvallet_reprod/files/results.txt",
      append=TRUE)

cat('\nFinished training and evaluation the model in ', 
    proc.time()[3] - start.time, 's...\n', sep='')
