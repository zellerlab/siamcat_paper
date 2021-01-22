# ##############################################################################
#
## IBD meta-analysis and ML pitfalls
#
# ##############################################################################

library("tidyverse")
library("SIAMCAT")
library("pROC")
library("ggembl")
library("progress")
library("here")

# functions needed
source(here('utils', 'train_ctr_augm.R'))
source(here('utils', 'load_ctr_datasets.R'))

# ##############################################################################
# prepare dataset 
datasets <- c('metaHIT', 'Lewis_2015', 'He_2017', 'Franzosa_2019', 'HMP2')
tax.info <- read_tsv(here('data', 'motus_taxonomy.tsv'))
if (!file.exists(here("ibd_meta_analysis", "data", "feat_genus.tsv"))){
    
  meta <- list()
  feat <- list()
  
  for (d in datasets){
    # metadata
    fn.meta <- here('data', 'meta', paste0('meta_', d, '.tsv'))
    meta.temp <- read_tsv(fn.meta) 
    if ("Country" %in% colnames(meta.temp)){
      if ('spanish' %in% meta.temp$Country){
        meta.temp <- meta.temp %>% 
          filter(Country=='spanish')
      }  
    }
    if ('Subject' %in% colnames(meta.temp)){
      meta.temp <- meta.temp %>% 
        mutate(Individual_ID=paste0("ID_", Subject), Timepoint=Time)
    }
    if ("Individual_ID" %in% colnames(meta.temp)){
      if ('Timepoint' %in% colnames(meta.temp)){
        meta.temp <- meta.temp %>% 
          select(Sample_ID, Group, Library_Size, Individual_ID, Timepoint)
      } else {
        meta.temp <- meta.temp %>% 
          mutate(Timepoint=Sampling_day) %>% 
          select(Sample_ID, Group, Library_Size, Individual_ID, Timepoint)
      }
    } else {
      meta.temp <- meta.temp %>% 
        select(Sample_ID, Group, Library_Size) %>% 
        mutate(Individual_ID=Sample_ID, Timepoint=0)
    }
    if ('nonIBD' %in% meta.temp$Group){
      meta.temp <- meta.temp %>% 
        mutate(Group=case_when(Group=='nonIBD'~'CTR',
                               TRUE~Group))
    }
    
    meta.temp <- meta.temp %>% 
      filter(Group %in% c('CTR', 'CD')) %>% 
      mutate(Study=d)
    meta.temp$Sample_ID <- make.names(meta.temp$Sample_ID)
    meta[[d]] <- meta.temp
  
    # features
    fn.feat <- here('data', 'features', 'motus', paste0(d, '_motus.tsv'))
    feat.temp <- read.table(fn.feat, sep='\t', stringsAsFactors = FALSE,
                            check.names = FALSE, quote = '', comment.char = '')
    stopifnot(nrow(feat.temp) == 14213)
    feat.temp <- feat.temp[,colSums(feat.temp) > 100]
    feat.temp <- prop.table(as.matrix(feat.temp), 2)
    
    # remove unmapped 
    feat.temp <- feat.temp[-which(rownames(feat.temp)=='-1'),]
    
    # change rownames to mOTU ids
    rownames(feat.temp) <- str_extract(rownames(feat.temp),
                                       '(ref|meta)_mOTU_v25_[0-9]{5}')
    feat[[d]] <- feat.temp
  }
  
  meta.all <- bind_rows(meta)
  feat.all <- do.call(cbind, feat)
  
  meta.all <- meta.all %>% 
    filter(Sample_ID %in% colnames(feat.all))
  
  # filter by individual for holdout testing
  meta.ind <- meta.all %>% 
    group_by(Individual_ID) %>% 
    filter(Timepoint==min(Timepoint)) %>% 
    ungroup()
  
  feat.all <- feat.all[,meta.all$Sample_ID]
  
  # convert to genus level abundances
  tax.info <- tax.info %>% 
    filter(!str_detect(family, '^NA'))
  bins.unique <- unique(tax.info$genus)
  feat.genus <- matrix(NA, nrow=length(bins.unique), ncol=ncol(feat.all),
                            dimnames=list(bins.unique, colnames(feat.all)))
  pb <- progress_bar$new(total=length(bins.unique))
  for (x in bins.unique){
    feat.genus[x,] <- colSums(feat.all[tax.info %>% 
                                              filter(genus==x) %>% 
                                              pull(mOTU_ID),,drop=FALSE])
    pb$tick()
  }
  
  # save feature and metadata
  write_tsv(meta.ind, path = here('ibd_meta_analysis', 
                                  'data', 'meta_individual.tsv'))
  write_tsv(meta.all, path = here('ibd_meta_analysis', 
                                  'data', 'meta_all.tsv'))
  write.table(feat.all, file = here('ibd_meta_analysis', 'data', 
                                    'feat_all_rel.tsv'), sep='\t', 
              quote = FALSE)
  write.table(feat.genus, file = here('ibd_meta_analysis', 'data', 
                                      'feat_genus.tsv'), sep='\t', 
              quote = FALSE)
} else {
  feat.genus <- read.table(here('ibd_meta_analysis', 'data', 'feat_genus.tsv'),
                           sep='\t', quote = '', comment.char = '', 
                           stringsAsFactors = FALSE, check.names = FALSE)
  meta.all <- read_tsv(here('ibd_meta_analysis', 'data', 'meta_all.tsv'))
  meta.ind <- read_tsv(here('ibd_meta_analysis', 'data', 
                            'meta_individual.tsv'))
}

ml.method <- 'enet-0.5'
stopifnot(ml.method %in% c('lasso', 'enet', 'enet-0.5'))

# ##############################################################################
# filter by prevalence
prev.per.study <- vapply(rownames(feat.genus), function(x){
  vapply(datasets, function(y){
    mean(feat.genus[x,meta.ind %>% 
                      filter(Study==y) %>% 
                      pull(Sample_ID)] != 0)
  }, FUN.VALUE = double(1))
}, FUN.VALUE = double(length(datasets)))
prev.per.study <- t(prev.per.study)
f.idx <- which(rowMeans(prev.per.study > 0.05) > 0.7)

feat.genus.filt <- feat.genus[f.idx,]

# ##############################################################################
# load control datasets
feat.ctr <- .f_load_control_data('cohort')$feat.ctr
# collapse at genus level
for (x in seq_len(length(feat.ctr))){
  ctr.data <- feat.ctr[[x]]
  rownames(ctr.data) <- str_extract(rownames(ctr.data), 
                                    '(ref|meta)_mOTU_v25_[0-9]{5}')
  feat.g <- matrix(NA, nrow=length(bins.unique), 
                   ncol=ncol(ctr.data),
                   dimnames=list(bins.unique, colnames(ctr.data)))
  pb <- progress_bar$new(total=length(bins.unique))
  for (g in bins.unique){
    feat.g[g,] <- colSums(ctr.data[tax.info %>% 
                                     filter(genus==g) %>% 
                                     pull(mOTU_ID),,drop=FALSE])
    pb$tick()
  }
  feat.ctr[[x]] <- feat.g
}

# ##############################################################################
# perform machine learning meta-analysis
fn.pred <- here("ibd_meta_analysis", "predictions", 
                paste0("predictions_", ml.method, ".tsv"))
if (!file.exists(fn.pred)){
  
  pred.matrix <- matrix(NA, nrow=nrow(meta.all), ncol=length(datasets),
                          dimnames = list(meta.all$Sample_ID, datasets))
  for (i in datasets){
    meta.train <- meta.all %>% 
      filter(Study==i) %>% 
      as.data.frame()
    rownames(meta.train) <- meta.train$Sample_ID
    
    insp <- NULL
    if (i %in% c('Lewis_2015', 'HMP2')){
      insp <- 'Individual_ID'
      if (i == 'HMP2'){
        meta.train <- meta.all %>% 
          filter(Study=='HMP2') %>% 
          group_by(Individual_ID) %>% 
          sample_n(5, replace = TRUE) %>% 
          distinct() %>% 
          as.data.frame()
        rownames(meta.train) <- meta.train$Sample_ID
      }
    }

    feat.train <- feat.genus.filt[,meta.train$Sample_ID]
    sc.obj.train <- siamcat(feat=feat.train, meta=meta.train, 
                            label='Group', case='CD', verbose=0)
    sc.obj.train <- normalize.features(sc.obj.train, norm.method = 'log.std',
                                       norm.param=list(log.n0=1e-05, 
                                                       sd.min.q=0),
                                       feature.type = 'original', verbose=0)
    sc.obj.train <- create.data.split(sc.obj.train,
                                      num.folds = 10, num.resample = 10,
                                      inseparable = insp, verbose=0)
    if (ml.method=='enet-0.5'){
      sc.obj.train <- train.model.ctr(sc.obj.train, method=ml.method, 
                                      correct.list = feat.ctr)
    } else {
      sc.obj.train <- train.model(sc.obj.train, method=ml.method)
    }
    message("+ Trained model for ", i)
    for (i2 in datasets){
      if (i == i2){
        sc.obj.train <- make.predictions(sc.obj.train, verbose=0)
        temp <- pred_matrix(sc.obj.train)
        pred.matrix[rownames(temp), i] <- rowMeans(temp)
      } else {
        meta.test <- meta.ind %>%
          filter(Study==i2) %>%
          as.data.frame()
        rownames(meta.test) <- meta.test$Sample_ID
        feat.test <- feat.genus[,meta.test$Sample_ID]
        sc.obj.test <- siamcat(feat=feat.test, meta=meta.test,
                                label='Group', case='CD', verbose=0)
        sc.obj.test <- make.predictions(sc.obj.train, sc.obj.test, verbose=0)
        temp <- pred_matrix(sc.obj.test)
        pred.matrix[rownames(temp), i] <- rowMeans(temp)
      }
    }
    message("+ Finished predictions for ", i)
    save(sc.obj.train,
         file = here("ibd_meta_analysis", "models",
                     paste0("trained_model_", i, "_",
                            ml.method, "_augmented.RData")))
  }
  
  pred.matrix <- as.data.frame(pred.matrix)
  pred.matrix$Sample_ID <- rownames(pred.matrix)
  pred.matrix <- as_tibble(pred.matrix)
  write_tsv(pred.matrix, path = fn.pred)
} else {
  pred.matrix <- read_tsv(fn.pred)
}

# ##############################################################################
# evaluate predictions
auroc.all <- tibble()

for (study.test in datasets){
  temp <- meta.ind %>% 
    filter(Study==study.test)
  df.temp <- inner_join(pred.matrix, temp, by='Sample_ID')
  for (study.train in c(datasets)){
    temp <- roc(predictor=df.temp[[study.train]],
                response=df.temp$Group,
                levels = c('CD', 'CTR'), direction = '>')
    auroc.all <- bind_rows(auroc.all, 
                           tibble(study.train=study.train, 
                                  study.test=study.test,
                                  AUC=c(temp$auc)))
  }
}

for (j in c('HMP2', 'metaHIT', 'Lewis_2015')){
  load(here('ibd_meta_analysis', 'models', 
            paste0('trained_model_', j, '_', ml.method, '_augmented.RData')))
  sc.obj.train <- evaluate.predictions(sc.obj.train)
  auroc.all <- auroc.all %>% 
    mutate(AUC=case_when(study.train==j & study.test==j~
                           as.numeric(sc.obj.train@eval_data$auroc),
                         TRUE~AUC))
}
  

test.average <- auroc.all %>% 
  filter(study.train!=study.test) %>% 
  group_by(study.test) %>% 
  summarise(AUC=mean(AUC)) %>% 
  ungroup() %>% 
  mutate(study.train="Average")


g <- bind_rows(auroc.all, test.average) %>% 
  mutate(CV=study.train == study.test) %>%
  mutate(split=case_when(study.train=='LOSO'~'LOSO', 
                         study.train=='Average'~'Average', 
                         TRUE~'none')) %>% 
  mutate(split=factor(split, levels = c('none', 'Average', 'LOSO'))) %>% 
  mutate(study.train=factor(study.train, levels = c(datasets, 
                                                    'Average', 'LOSO'))) %>% 
  mutate(study.test=factor(study.test, levels = c(rev(datasets), 
                                                  'Average', 'LOSO'))) %>% 
  ggplot(aes(y=study.test, x=study.train, fill=AUC, size=CV, color=CV)) +
    geom_tile() + theme_publication() +
    # test in tiles
    geom_text(aes_string(label="format(AUC, digits=2)"), 
              col='white', size=2)+
    # color scheme
    scale_fill_gradientn(colours=rev(c('darkgreen','forestgreen', 
                                       'chartreuse3','lawngreen', 
                                       'yellow')), limits=c(0.5, 1)) +
    # axis position/remove boxes/ticks/facet background/etc.
    scale_x_discrete(position='top') + 
    theme(axis.line=element_blank(), 
          axis.ticks = element_blank(), 
          axis.text.x.top = element_text(angle=45, hjust=.1), 
          panel.grid=element_blank(), 
          panel.border=element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_blank()) + 
    xlab('Training Set') + ylab('Test Set') + 
    scale_color_manual(values=c('#FFFFFF00', 'grey'), guide=FALSE) + 
    scale_size_manual(values=c(0, 1), guide=FALSE) + 
    facet_grid(~split, scales = 'free', space = 'free')

ggsave(g, filename = here('figures', 'ibd_meta_analysis', 
                          paste0("performance_", ml.method,".pdf")),
       width = 120, height = 65, units = 'mm', useDingbats=FALSE)

# ##############################################################################
# association plot heatmap
assoc.list <- list()
weight.list <- list()
for (d in datasets){
  meta.train <- meta.ind %>% 
    filter(Study==d) %>% 
    as.data.frame()
  rownames(meta.train) <- meta.train$Sample_ID
  feat.train <- feat.genus.filt[,meta.train$Sample_ID]
  sc.obj.train <- siamcat(feat=feat.train, meta=meta.train, 
                          label='Group', case='CD')
  # associations
  sc.obj.train <- check.associations(sc.obj.train, detect.lim = 1e-04, 
                                     prompt = FALSE,
                                     feature.type = 'original', 
                                     sort.by = 'p.val')
  temp <- associations(sc.obj.train)
  temp$species <- rownames(temp)
  
  assoc.list[[d]] <- temp %>% 
    select(species, fc, auc, p.adj) %>% 
    mutate(Study=d)

  load(here('ibd_meta_analysis', 'models',
            paste0('trained_model_', d, '_', ml.method, '_augmented.RData')))

  temp <- feature_weights(sc.obj.train)
  temp$species <- rownames(temp)
  weight.list[[d]] <- temp %>%
    select(species, median.rel.weight, mean.rel.weight, percentage) %>%
    mutate(Study=d) %>%
    mutate(r.med=rank(-abs(median.rel.weight)),
           r.mean=rank(-abs(mean.rel.weight)))

}
df.assoc <- bind_rows(assoc.list)
df.assoc <- df.assoc %>% filter(species!='unclassified')

df.weights <- bind_rows(weight.list)
df.weights <- df.weights %>% filter(species!='unclassified')

x <- df.assoc %>% 
  group_by(species) %>% 
  summarise(m=mean(auc), n.filt=any(auc < 0.25 | auc > 0.75)) %>% 
  filter(n.filt) %>% 
  arrange(m) %>% 
  filter(species!='unclassified')


g <- df.assoc %>% 
  filter(species %in% x$species) %>% 
  mutate(species=factor(species, levels = rev(x$species))) %>% 
  mutate(Study=factor(Study, levels = datasets)) %>% 
  mutate(l=case_when(p.adj < 0.01~'*', TRUE~'')) %>% 
  ggplot(aes(y=species, x=Study, fill=fc)) + 
    geom_tile() + 
    scale_fill_gradient2_embl(palette = 'Blue-Red', limits=c(-2.57, 2.57)) + 
    theme_minimal() + 
    geom_text(aes(label=l)) +
    theme(panel.grid = element_blank()) + 
    xlab('') + ylab('') +
    theme(axis.text = element_text(size=6))
ggsave(g, filename = here('figures', 'ibd_meta_analysis', 'heatmap_ibd.pdf'),
       width = 120, height = 80, units = 'mm', useDingbats=FALSE)

abs.weights <- df.weights %>% 
  group_by(Study) %>% 
  summarise(s.med=sum(abs(median.rel.weight)),
            s.mean=sum(abs(mean.rel.weight)))

g1 <- df.weights %>% 
  full_join(abs.weights) %>% 
  mutate(median.rel.weight=median.rel.weight/s.med,
         mean.rel.weight=mean.rel.weight/s.mean) %>% 
  mutate(r.mean=case_when(r.mean > 20~NA_real_, TRUE~r.mean)) %>%
  mutate(r.med=case_when(r.med > 20~NA_real_, TRUE~r.med)) %>%
  mutate(incl=species %in% x$species) %>% 
  group_by(Study, incl) %>% 
  summarise(w=sum(abs(median.rel.weight))) %>% 
  filter(incl) %>% 
  ungroup() %>% 
  mutate(Study=factor(Study, levels = datasets)) %>% 
  ggplot(aes(x=Study, y=w)) + 
    geom_bar(stat="identity") + 
    theme_publication() + 
    xlab('') + 
    ylab('Proportion of model shown')


g2 <- df.weights %>% 
  full_join(abs.weights) %>% 
  mutate(median.rel.weight=median.rel.weight/s.med,
         mean.rel.weight=mean.rel.weight/s.mean) %>% 
  mutate(r.mean=case_when(r.mean > 20~NA_real_, TRUE~r.mean)) %>%
  mutate(r.med=case_when(r.med > 20~NA_real_, TRUE~r.med)) %>%
  filter(species %in% x$species) %>% 
  mutate(species=factor(species, levels = rev(x$species))) %>% 
  mutate(Study=factor(Study, levels = datasets)) %>% 
  ggplot(aes(y=species, x=Study, fill=median.rel.weight)) + 
    geom_tile() + 
    scale_fill_gradientn(colours=rev(c('#007A53', '#009F4D', "#6CC24A", 'white',
                                   "#EFC06E", "#FFA300", '#BE5400')), 
                         limits=c(-0.15, 0.15)) +
    theme_minimal() + 
    geom_text(aes(label=r.med), col='black', size= 2) +
    theme(panel.grid = element_blank()) + 
    xlab('') + ylab('') +
    theme(axis.text = element_text(size=6))

ggsave(g2, filename = here('figures', 'ibd_meta_analysis', 
                           paste0('weights_heatmap_', ml.method ,'.pdf')),
       width = 150, height = 80, units = 'mm', useDingbats=FALSE)

# order of the genera
tax.info %>% 
  filter(genus %in% x$species) %>% 
  select(genus, order) %>% 
  distinct() %>% 
  mutate(genus=factor(genus, levels = rev(x$species))) %>% 
  ggplot(aes(x=genus, y=1, fill=order)) + 
    geom_tile() + 
    coord_flip()
  
# ##############################################################################
# confounder plot for study/disease?

feat.conf <- feat.genus[,meta.ind$Sample_ID]
label.conf <- meta.ind$Group
names(label.conf) <- meta.ind$Sample_ID
study <- as.factor(meta.ind$Study)
names(study) <- meta.ind$Sample_ID

var.label <- vapply(rownames(feat.conf), FUN=function(x){
  x <- feat.conf[x,]
  x <- rank(x)/length(x)
  ss.tot <- sum((x - mean(x))^2)/length(x)
  ss.o.i <- sum(vapply(unique(label.conf), function(s){
    sum((x[label.conf==s] - mean(x[label.conf==s]))^2)
  }, FUN.VALUE = double(1)))/length(x)
  return(1-ss.o.i/ss.tot)
}, FUN.VALUE = double(1))
if (any(is.infinite(var.label))){
  var.label[is.infinite(var.label)] <- NA
}
var.batch <- vapply(rownames(feat.conf), FUN=function(x){
  x <- feat.conf[x,names(study)]
  x <- rank(x)/length(x)
  ss.tot <- sum((x - mean(x))^2)/length(x)
  ss.o.i <- sum(vapply(levels(study), function(s){
    sum((x[study==s] - mean(x[study==s]))^2)
  }, FUN.VALUE = double(1)))/length(x)
  return(1-ss.o.i/ss.tot)
}, FUN.VALUE = double(1))
if (any(is.infinite(var.batch))){
  var.batch[is.infinite(var.batch)] <- NA
}
df.plot <- tibble(label=var.label, batch=var.batch, species=names(var.label))
df.plot$mean <- rowMeans(log10(feat.conf+1e-05))

g.conf <- df.plot %>% 
  mutate(selection=case_when(species %in% c(x %>% filter(m < 0.5) %>% 
                                              pull(species))~'blue',
                             species %in% c(x %>% filter(m > 0.5) %>% 
                                              pull(species))~'red',
                             TRUE~'grey')) %>%
  ggplot(aes(x=label, y=batch, size=mean, col=selection)) + 
  geom_point(stroke=0) + 
  theme_publication()  +
  theme(aspect.ratio = 1) + 
  scale_size_continuous(range = c(0.2, 4)) + 
  xlim(0, 0.4) + ylim(0, 0.4) + 
  scale_colour_manual(values=c('blue'='#3B6FB660',
                               'grey'='#70737230', 
                               'red'='#D4164560')) +
  geom_abline(intercept = 0, slope=1, col='darkgrey', linetype=3)
ggsave(g.conf, filename = here("figures", "ibd_meta_analysis", 
                               "conf_plot_ibd.pdf"),
       width = 80, height = 80, useDingbats=FALSE, units = 'mm')
