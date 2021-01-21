# ##############################################################################
#
## Ocean metagenomes for classification
#
# ##############################################################################

library("tidyverse")
library("SIAMCAT")
library("readxl")
library("ggembl")
library("here")

# ##############################################################################
# get data

# metadata
meta <- read_excel(here('data', 'ocean', 'ocean.xlsx'),
                   sheet = 8) %>%
  filter(Layer=='SRF') %>%
  separate(col=`PANGAEA sample id`, into=c('metaG', 'metaT'), sep='-')

# features
fn.feat <- here('data', 'ocean',
                'mitags_tab_genus.tsv.gz')
feat <- read.table(fn.feat,
                   stringsAsFactors = FALSE, check.names = FALSE,
                   quote = '', sep='\t', skip = 1)
rownames(feat) <- feat$V1
feat$V1 <- NULL
feat <- as.matrix(feat)
feat <- t(feat)
# adjust feature names
otu.names <- read_lines(fn.feat, n_max = 1)
otu.names <- str_split(otu.names, '\t')[[1]]
otu.names <- otu.names[-1]
if (str_detect(fn.feat, 'otu')){
  otu.names.clean <- str_split_fixed(otu.names, ' ', n = 2)
  rownames(feat) <- otu.names.clean[,1]
} else {
  rownames(feat) <- otu.names
  # otu.names.clean <- otu.names # str_split_fixed(otu.names, ';', n=6)
}

feat.rel <- prop.table(feat, 2)

# ##############################################################################
# train polar vs non-polar
meta.mG <- meta %>%
  as.data.frame()
rownames(meta.mG) <- meta.mG$metaG

sc.obj <- siamcat(feat=feat.rel, meta=meta.mG, label='polar', case='polar')
hist(log10(matrixStats::rowMaxs(get.orig_feat.matrix(sc.obj))), 100)
hist(rowMeans(get.orig_feat.matrix(sc.obj) != 0), 100)
sc.obj <- filter.features(sc.obj, filter.method = 'prevalence',
                          cutoff=0.05, verbose=3)
sc.obj <- check.associations(sc.obj, detect.lim = 1e-06, prompt = FALSE)
## the data split is a bit more comlicated, since we want to block by
##  ocean region, but this information is folded into the label as well, since 
##  only artic ocean regions are polar. So we have to split the polar regions 
##  normally, but the other ocean regions in a blocked way and then combine the
##  data splits
meta.p <- meta.mG[meta.mG$polar=='polar',]
meta.np <- meta.mG[meta.mG$polar=='non polar',]
assign.folds.ocean <- function(meta, num.folds, inseparable = NULL){
  
  foldid <- rep(0, nrow(meta))
  
  # If stratify is not TRUE, make sure that num.sample is not
  # bigger than number.folds
  if (nrow(meta) <= num.folds) {
    warning("+++ num.samples is exceeding number of folds, setting
            CV to (k-1) unstratified CV")
    num.folds <- nrow(meta) - 1
  }
  if (!is.null(inseparable)) {
    strata <- unique(meta[[inseparable]])
    sid <-
      sample(rep(seq_len(num.folds), length.out =
                   length(strata)))
    for (s in seq_along(strata)) {
      idx <- which(meta[[inseparable]] == strata[s])
      foldid[idx] <- sid[s]
    }
    stopifnot(all(!is.na(foldid)))
  } else {
    foldid <- sample(rep(seq_len(num.folds),
                         length.out = nrow(meta)))
  }
  names(foldid) <- rownames(meta)
  return(foldid)
}

train.list <- list(NULL)
test.list <- list(NULL)
for (r in seq_len(5)){
  train.temp <- list(NULL)
  test.temp <- list(NULL)
  foldid.p <- assign.folds.ocean(meta=meta.p, num.folds=5)
  foldid.np <- assign.folds.ocean(meta=meta.np, num.folds=5, 
                                   inseparable='Ocean.region')
  for (f in seq_len(5)){
    # select test examples
    test.idx <- c(names(which(foldid.p == f)), names(which(foldid.np == f)))
    train.idx <- c(names(which(foldid.p != f)), names(which(foldid.np != f)))
    
    train.temp[f] <- list(train.idx)
    test.temp[f] <- list(test.idx)
    stopifnot(length(intersect(train.idx, test.idx)) == 0)
  }
  train.list[[r]] <- train.temp
  test.list[[r]] <- test.temp
}

data_split(sc.obj) <- list(
  training.folds = train.list,
  test.folds = test.list,
  num.resample = 5,
  num.folds = 5
)

sc.obj <- normalize.features(sc.obj, norm.method = 'log.std',
                             norm.param = list(log.n0=1e-06, sd.min.q=0))
sc.obj <- train.model(sc.obj, method = 'enet')
sc.obj <- make.predictions(sc.obj)
sc.obj <- evaluate.predictions(sc.obj)
model.evaluation.plot(sc.obj,
                      fn.plot = here('figures', 'ocean', 'eval_plot.pdf'))
model.interpretation.plot(sc.obj, consens.thres = 0.8,
                          fn.plot = here('figures',
                                         'ocean', 'interpret_plot.pdf'))
df.plot <- tibble(metaG=rownames(pred_matrix(sc.obj)),
                  pred=rowMeans(pred_matrix(sc.obj)))
df.plot <- left_join(df.plot, meta) %>%
  select(polar, pred) %>%
  mutate(type='metaG') %>%
  group_by(polar) %>%
  mutate(x_value=rank(pred)/n()) %>%
  ungroup()

# show that this also works on metaT samples
meta.mT <- meta %>%
  as.data.frame()
rownames(meta.mT) <- meta.mT$metaT

sc.obj.mT <- siamcat(feat=feat.rel, meta=meta.mT)
sc.obj.mT <- make.predictions(sc.obj, sc.obj.mT)
pred.mT <- tibble(metaT=rownames(pred_matrix(sc.obj.mT)),
                  pred=rowMeans(pred_matrix(sc.obj.mT)))
df.plot.mT <- left_join(pred.mT, meta) %>%
  select(polar, pred) %>%
  mutate(type='metaT') %>%
  group_by(polar) %>%
  mutate(x_value=rank(pred)/n()) %>%
  ungroup()

df.plot.all <- bind_rows(df.plot, df.plot.mT)
thres <- sc.obj@eval_data$roc$thresholds[
  which(sc.obj@eval_data$roc$specificities > 0.9)[1]]
g <- df.plot.all %>%
  ggplot(aes(x=x_value, y=pred)) +
  geom_point() +
  facet_grid(type~polar) +
  geom_hline(yintercept = thres) +
  theme_publication(panel.grid = 'major') +
  ylab('Model prediction') +
  xlab('Relative rank') +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0))
ggsave(g, filename = here('figures', 'ocean', 'pred_to_mGmT.pdf'),
       width = 100, height = 70, units = 'mm', useDingbats=FALSE)
save(sc.obj, file = here('ocean', 'model_polar.RData'))


# ##############################################################################
# train high vs low iron

meta.mG <- meta %>%
  mutate(label=case_when(Iron.5m > 10^(-3.7) ~'high', TRUE~'low')) %>%
  as.data.frame()
rownames(meta.mG) <- meta.mG$metaG

sc.obj.i <- siamcat(feat=feat.rel, meta=meta.mG, label='label', case='low')
hist(log10(matrixStats::rowMaxs(get.orig_feat.matrix(sc.obj.i))), 100)
hist(rowMeans(get.orig_feat.matrix(sc.obj.i) != 0), 100)
sc.obj.i <- filter.features(sc.obj.i, filter.method = 'prevalence',
                            cutoff=0.05, verbose=3)
sc.obj.i <- check.associations(sc.obj.i, detect.lim = 1e-06, prompt = FALSE)

meta.p <- meta.mG[meta.mG$label=='low',]
meta.np <- meta.mG[meta.mG$label=='high',]

train.list <- list(NULL)
test.list <- list(NULL)
for (r in seq_len(5)){
  train.temp <- list(NULL)
  test.temp <- list(NULL)
  foldid.p <- assign.folds.ocean(meta=meta.p, num.folds=5)
  foldid.np <- assign.folds.ocean(meta=meta.np, num.folds=5, 
                                  inseparable='Ocean.region')
  for (f in seq_len(5)){
    # select test examples
    test.idx <- c(names(which(foldid.p == f)), names(which(foldid.np == f)))
    train.idx <- c(names(which(foldid.p != f)), names(which(foldid.np != f)))
    
    train.temp[f] <- list(train.idx)
    test.temp[f] <- list(test.idx)
    stopifnot(length(intersect(train.idx, test.idx)) == 0)
  }
  train.list[[r]] <- train.temp
  test.list[[r]] <- test.temp
}

data_split(sc.obj.i) <- list(
  training.folds = train.list,
  test.folds = test.list,
  num.resample = 5,
  num.folds = 5
)



sc.obj.i <- normalize.features(sc.obj.i, norm.method = 'log.std',
                               norm.param = list(log.n0=1e-06, sd.min.q=0))
sc.obj.i <- train.model(sc.obj.i, method = 'enet')
sc.obj.i <- make.predictions(sc.obj.i)
sc.obj.i <- evaluate.predictions(sc.obj.i)
model.evaluation.plot(sc.obj.i, fn.plot = here('figures',
                                               'ocean', 'eval_plot_iron.pdf'))
model.interpretation.plot(sc.obj.i, consens.thres = 0.5,
                          fn.plot = here('figures',
                                         'ocean', 'interpret_plot_iron.pdf'))
df.plot <- tibble(metaG=rownames(pred_matrix(sc.obj.i)),
                  pred=rowMeans(pred_matrix(sc.obj.i)))
df.plot <- left_join(df.plot, meta) %>%
  mutate(label=case_when(Iron.5m > 10^(-3.5) ~'high', TRUE~'low')) %>%
  select(label, pred) %>%
  mutate(type='metaG') %>%
  group_by(label) %>%
  mutate(x_value=rank(pred)/n()) %>%
  ungroup()

# show that this also works on metaT samples
meta.mT <- meta %>%
  as.data.frame()
rownames(meta.mT) <- meta.mT$metaT

sc.obj.mT <- siamcat(feat=feat.rel, meta=meta.mT)
sc.obj.mT <- make.predictions(sc.obj.i, sc.obj.mT)
pred.mT <- tibble(metaT=rownames(pred_matrix(sc.obj.mT)),
                  pred=rowMeans(pred_matrix(sc.obj.mT)))
df.plot.mT <- left_join(pred.mT, meta) %>%
  mutate(label=case_when(Iron.5m > 10^(-3.5) ~'high', TRUE~'low')) %>%
  select(label, pred) %>%
  mutate(type='metaT') %>%
  group_by(label) %>%
  mutate(x_value=rank(pred)/n()) %>%
  ungroup()

df.plot.all <- bind_rows(df.plot, df.plot.mT)
thres <- sc.obj.i@eval_data$roc$thresholds[
  which(sc.obj.i@eval_data$roc$specificities > 0.9)[1]]
g <- df.plot.all %>%
  ggplot(aes(x=x_value, y=pred)) +
  geom_point() +
  facet_grid(type~label) +
  geom_hline(yintercept = thres) +
  theme_publication(panel.grid = 'major') +
  ylab('Model prediction') +
  xlab('Relative rank') +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0))
ggsave(g, filename = here('figures', 'ocean', 'pred_to_mGmT_iron.pdf'),
       width = 100, height = 70, units = 'mm', useDingbats=FALSE)

save(sc.obj.i, file = here('ocean', 'model_iron.RData'))