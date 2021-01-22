# ##############################################################################
#
## Check how similar the datasets are towards each other
#
# ##############################################################################

library("tidyverse")
library("here")
library("SIAMCAT")
library("ggembl")

# load info about datasets
motu.tasks <- read_tsv(here('parameter_space', 'data_info', 
                            'all_tasks.tsv')) %>% 
  filter(type=='mOTUs2') %>% 
  select(dataset.id) %>% 
  distinct()

fn.res <- here('control_augmentation', 'files', 'ctr_auroc.tsv')
if (!file.exists(fn.res)){
  df.res <- tibble(Study.1=character(0), Study.2=character(0), AUC=double(0))

  for (i in seq_len((nrow(motu.tasks)-1))){
    dataset.1 <- motu.tasks$dataset.id[i]
    x <- list.files(here('parameter_space', 'sc'), pattern=dataset.1) %>% 
      grep(pattern='mOTUs', value=TRUE) %>% 
      sample(size=1)
    load(here('parameter_space', 'sc', x))
    feat.1 <- get.orig_feat.matrix(sc.obj)
    meta.1 <- as.data.frame(meta(sc.obj))
    meta.1$Study <- dataset.1
    meta.1 <- meta.1[meta.1$Group=='CTR',]
    meta.1 <- meta.1[,c('Study', 'Group')]
    
    for (j in seq(from=(i+1), to=nrow(motu.tasks))){
      dataset.2 <- motu.tasks$dataset.id[j]
      message(dataset.1, '-', dataset.2)
      y <- list.files(here('parameter_space', 'sc'), pattern=dataset.2) %>% 
        grep(pattern='mOTUs', value=TRUE) %>% 
        sample(size=1)
      load(here('parameter_space', 'sc', y))
      feat.2 <- get.orig_feat.matrix(sc.obj)
      meta.2 <- as.data.frame(meta(sc.obj))
      meta.2$Study <- dataset.2
      meta.2 <- meta.2[meta.2$Group=='CTR',]
      meta.2 <- meta.2[,c('Study', 'Group')]
      
      # train classifier
      feat.all <- cbind(feat.1, feat.2)
      meta.all <- rbind(meta.1, meta.2)
      
      sc.obj <- siamcat(feat=feat.all, meta=meta.all, 
                        label='Study', case=dataset.1,
                        verbose=0)
      sc.obj <- filter.features(sc.obj, filter.method = 'abundance', 
                                cutoff=0.001, verbose=0)
      sc.obj <- normalize.features(sc.obj, norm.method='log.std', 
                                   norm.param = list(log.n0=1e-05, sd.min.q=0),
                                   verbose=0)
      sc.obj <- create.data.split(sc.obj, num.folds = 5, num.resample = 5, 
                                  verbose=0)
      sc.obj <- train.model(sc.obj, method='lasso',verbose=0)
      sc.obj <- make.predictions(sc.obj, verbose=0)
      sc.obj <- evaluate.predictions(sc.obj, verbose=0)
      df.res <- df.res %>% 
        add_row(Study.1=dataset.1, Study.2=dataset.2, 
                AUC=as.numeric(eval_data(sc.obj)$auroc))
    }
  }

  write_tsv(df.res, fn.res)
} else {
  df.res <- read_tsv(fn.res)
}
  

g <- df.res %>% 
  mutate(AUC=case_when(AUC > 0.5~AUC, TRUE~1-AUC)) %>% 
  ggplot(aes(x=Study.1, y=Study.2, fill=AUC)) + 
  geom_tile() + theme_poster() + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) + 
  geom_text(aes(label=sprintf('%.2f', AUC))) + 
  scale_fill_continuous_embl(limits=c(0.5, 1)) + 
  xlab('Study 1') + 
  ylab("Study 2")
ggsave(g, filename = here('figures', 'cross_prediction', 
                          'study_similarity.pdf'), 
       width = 12, height = 9, useDingbats=FALSE)
