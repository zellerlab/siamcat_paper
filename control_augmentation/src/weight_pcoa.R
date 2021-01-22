# ##############################################################################
#
## Re-train best parameter set models and extract feature weights
##    and viz feature weights
#
# ##############################################################################

library("here")
library("tidyverse")
library("SIAMCAT")
library("vegan")
library("labdsv")
library("ggembl")
library("yaml")
library("pheatmap")
library("ape")
library("ggrepel")
library("matrixStats")
library("cowplot")

colours <- yaml.load_file(here("parameter_space", "data_info", "colours.yaml"))
disease_colours <- unlist(colours$disease)
motus.taxonomy <- read_tsv(here("data", "motus_taxonomy.tsv"))
motus.tree <- read.tree('~/Documents/utils/general_data/mOTUS_tree_2.5.nwx')

# ##############################################################################
# load data sets
ml.method <- 'enet-0.5'
stopifnot(ml.method %in% c('enet-0.5', 'enet', 'lasso'))
motu.tasks <- read_tsv(here("parameter_space", "files", "cross_predictions",
                            paste0("cross_prediction_", ml.method,
                                   ".tsv"))) %>% 
  select(dataset.id.train, case.train, auroc, augmented) %>% 
  distinct() %>% 
  filter(augmented, auroc > 0.7) 

# ##############################################################################
# loop through datasets
# extract feature weights of best and augmented models
# perform differential abundance testing and extract the results
if (!file.exists(here("parameter_space", "files", "weights", "weights.tsv"))){
  
  weight.list <- list()
  weight.list.augm <- list()
  assoc.list <- list()
  assoc.list.g <- list()
  assoc.list.full <- list()
  assoc.list.g.full <- list()
  for (i in seq_len(nrow(motu.tasks))){
    for (ml.method in c('enet', 'lasso', 'enet-0.5')){
      tag <- paste0(motu.tasks$dataset.id.train[i], '_', 
                    motu.tasks$case.train[i],
                    '_', ml.method)
      message(tag)
      # normal models
      fn.model.trained <- here('parameter_space', 'models', 'naive',
                               paste0('sc_trained_', tag,  '.RData'))
      if (!file.exists(fn.model.trained)){
        message("No model found for ",
                paste0(motu.tasks$dataset.id.train[i], '_', 
                       motu.tasks$case.train[i]))
        next()
      }
      load(fn.model.trained)
      # extract features
      weight.mat <- weight_matrix(sc.obj.train)
      colnames(weight.mat) <- paste0(motu.tasks$dataset.id.train[i], '-',
                                     motu.tasks$case.train[i], '-original-',
                                     ifelse(ml.method=='enet-0.5',
                                            'enet0.5', ml.method), '-',
                                     seq_len(ncol(weight.mat)))
      weight.mat <- as_tibble(weight.mat, rownames='mOTU_ID') %>%
        pivot_longer(cols=-1, names_to='model', values_to='weight')
      weight.list[[tag]] <- weight.mat
    }
    
    # associations
    tag.assoc <- paste0(motu.tasks$dataset.id.train[i], '-', 
                        motu.tasks$case.train[i])
    load(here('parameter_space', 'models', "naive", 
              paste0('sc_trained_', motu.tasks$dataset.id.train[i],
                     '_', motu.tasks$case.train[i], "_enet.RData")))
    sc.obj.train <- check.associations(sc.obj.train, 
                                       fn.plot = here('temp.pdf'),
                                       verbose=0)
    assoc.mat <- as_tibble(associations(sc.obj.train), rownames='mOTU_ID') %>% 
      mutate(Study=tag.assoc)
    assoc.list.full[[tag.assoc]] <- assoc.mat %>% select(fc, mOTU_ID, Study)
    assoc.list[[tag.assoc]] <- assoc.mat %>% 
      filter(p.adj < 0.2) %>% 
      select(fc, mOTU_ID, Study)
    
    # genus level associations
    included.motus <- enframe(rownames(get.filt_feat.matrix(sc.obj.train))) %>% 
      transmute(mOTU_name=value) %>% 
      mutate(mOTU_ID=str_extract(mOTU_name, '(ref|meta)_mOTU_v25_[0-9]{5}')) %>% 
      left_join(motus.taxonomy) %>% 
      select(-mOTU_ID, -specI_cluster, -mOTU) %>% 
      as.data.frame()
    rownames(included.motus) <- included.motus$mOTU_name
    included.motus$mOTU_name <- NULL
    included.motus <- as.matrix(included.motus)
    
    sc.test <- siamcat(feat=get.filt_feat.matrix(sc.obj.train), 
                       label = sc.obj.train@label, verbose = 0)
    tax_table(physeq(sc.test)) <- included.motus
    sc.test <- summarize.features(sc.test, level = 'genus', verbose=0)
    sc.test <- check.associations(sc.test, 
                                  fn.plot = here('temp.pdf'),
                                  verbose=0, alpha = 0.2,
                                  feature.type = 'original')
    assoc.mat <- as_tibble(associations(sc.test), rownames='mOTU_ID') %>% 
      mutate(Study=tag.assoc)
      
    assoc.list.g.full[[tag.assoc]] <- assoc.mat %>% select(fc, mOTU_ID, Study)
    assoc.list.g[[tag.assoc]] <- assoc.mat %>% 
      filter(p.adj < 0.2) %>% 
      select(fc, mOTU_ID, Study)
    
    # augmented models
    for (ml.method in c('enet', 'lasso', 'enet-0.5')){
      fraction <- 5
      fn.augm.model <- here('parameter_space', 'models', 'augmented',
                            paste0('sc_', motu.tasks$dataset.id.train[i],
                                   '-mOTUs2-', motu.tasks$case.train[i], '-',
                                   fraction, '-FALSE-', ml.method, '-FALSE',
                                   '_augmented.RData'))
      if (file.exists(fn.augm.model)){
        load(fn.augm.model)
        # extract features
        weight.mat <- weight_matrix(sc.obj)
        colnames(weight.mat) <- paste0(motu.tasks$dataset.id.train[i], '-',
                                       motu.tasks$case.train[i], '-augm', 
                                       fraction,
                                       '-',
                                       ifelse(ml.method=='enet-0.5',
                                              'enet0.5', ml.method), '-',
                                       seq_len(ncol(weight.mat)))
        weight.mat <- as_tibble(weight.mat, rownames='mOTU_ID') %>%
          pivot_longer(cols=-1, names_to='model', values_to='weight')
        weight.list.augm[[paste0(tag, '-augm', fraction,
                                 '-', ml.method)]] <- weight.mat
      } else {
        message("Problem with: ", motu.tasks$dataset.id.train[i], '-',
                motu.tasks$case.train[i], '-', ml.method)
      }
    }
  }
  
  weights.all <- bind_rows(weight.list) %>% 
    pivot_wider(names_from=model, values_from=weight, 
                values_fill=list('weight'=0))
  write_tsv(weights.all, path = here("parameter_space", "files", 
                                     "weights_naive.tsv"))
  weights.all <- bind_rows(weight.list.augm) %>% 
    pivot_wider(names_from=model, values_from=weight, 
                values_fill=list('weight'=0))
  write_tsv(weights.all, path = here("parameter_space", "files", 
                                     "weights_augm.tsv"))
  
  assoc.mat <- bind_rows(assoc.list) %>% 
    pivot_wider(names_from = Study, values_from = fc, 
                values_fill = list('fc'=0))
  write_tsv(assoc.mat, path = here("parameter_space", "files", 
                                   "assoc_mat.tsv"))
  assoc.mat.full <- bind_rows(assoc.list.full) %>% 
    pivot_wider(names_from = Study, values_from = fc, 
                values_fill = list('fc'=0))
  write_tsv(assoc.mat.full, path = here("parameter_space", "files", 
                                        "assoc_mat_full.tsv"))
  assoc.mat.g <- bind_rows(assoc.list.g) %>% 
    pivot_wider(names_from = Study, values_from = fc, 
                values_fill = list('fc'=0))
  write_tsv(assoc.mat.g, path = here("parameter_space", "files", 
                                   "assoc_mat_genus.tsv"))
  assoc.mat.g.full <- bind_rows(assoc.list.g.full) %>% 
    pivot_wider(names_from = Study, values_from = fc, 
                values_fill = list('fc'=0))
  write_tsv(assoc.mat.g.full, path = here("parameter_space", "files", 
                                          "assoc_mat_genus_full.tsv"))
  unlink(here('temp.pdf'))  

} else {
  weights.mat.naive <- read_tsv(here("parameter_space", "files", "weights",
                                     "weights_naive.tsv"))
  weights.mat.augm <- read_tsv(here("parameter_space", "files", "weights",
                                     "weights_augm.tsv"))
  assoc.mat <- read_tsv(here("parameter_space", "files", "assoc_mat.tsv"))
  assoc.mat.full <- read_tsv(here("parameter_space", "files",
                                  "assoc_mat_full.tsv"))
  assoc.mat.g <- read_tsv(here("parameter_space", "files", 
                               "assoc_mat_genus.tsv"))
  assoc.mat.g.full <- read_tsv(here("parameter_space", "files", 
                                    "assoc_mat_genus_full.tsv"))
}

# ##############################################################################
# General trafos of weights

# normalize for model size
weights.mat.naive.rel <- weights.mat.naive %>% 
  mutate_if(is_double, function(x){x/sum(abs(x))})
weights.mat.augm.rel <- weights.mat.augm %>% 
  mutate_if(is_double, function(x){x/sum(abs(x))})

# summarize as mean across models
weights.rel.naive.ind <- weights.mat.naive.rel %>% 
  pivot_longer(cols=-mOTU_ID) %>% 
  mutate(name=str_remove(name, '-[0-9]*$')) %>% 
  group_by(mOTU_ID, name) %>% 
  summarise(weight=mean(value)) %>% 
  pivot_wider(names_from = name, values_from = weight,
              values_fill = list('weight'=NA)) %>% 
  ungroup
weights.rel.augm.ind <- weights.mat.augm.rel %>% 
  pivot_longer(cols=-mOTU_ID) %>% 
  mutate(name=str_remove(name, '-[0-9]*$')) %>% 
  group_by(mOTU_ID, name) %>% 
  summarise(weight=mean(value)) %>% 
  pivot_wider(names_from = name, values_from = weight,
              values_fill = list('weight'=NA)) %>% 
  ungroup

# ##############################################################################
# PCoA of models
ml.method <- ifelse(ml.method == 'enet-0.5', 'enet0.5', ml.method)
weight.pco <- function(mat, id, dist.method='mod.jacc', indiv=TRUE){
  stopifnot(dist.method %in% c('mod.jacc', 'canberra'))
  stopifnot(id %in% c('lasso', 'enet-', 'enet0.5', '-'))
  if (class(indiv)=='logical'){
    if (id == 'enet-' && !indiv){
      id <- 'enet$'
    }
  }
  
  stopifnot('mOTU_ID' %in% colnames(mat))
  
  sel.cols <- which(str_detect(colnames(mat), id))
  mat.red <- mat %>% 
    select(c(mOTU_ID, sel.cols))
  # select only models that can actually manage something
  sel.models <- paste0(motu.tasks$dataset.id.train, '-', 
                       motu.tasks$case.train)
  sel.cols.2 <- unlist(map(sel.models, 
                           .f = function(x){which(
                             str_detect(colnames(mat.red), x))}))
  mat.red <- mat.red %>% 
    select(c(mOTU_ID, sel.cols.2))
  
  if (dist.method=='canberra'){
    temp <- mat.red %>% 
      as.data.frame()
    rownames(temp) <- temp$mOTU_ID
    temp$mOTU_ID <- NULL
    temp <- as.matrix(temp)
    dist.mat <- vegdist(t(temp), method='canberra')
  } else {
    list.d <- list()
    for (x in c('pos', 'neg')){
      if (x=='pos'){
        .f <- function(x){(ifelse(x < 0, 0, x))}
      } else {
        .f <- function(x){(ifelse(x > 0, 0, x))}
      }
      temp <- mat.red %>%
        select_if(is_double) %>% 
        mutate_all(.f) %>%
        as.data.frame() %>% as.matrix
      d <- vegdist(t(abs(temp)), method='jaccard', binary=FALSE)
      d <- as.matrix(d) %>%
        as_tibble(rownames = "Dataset_1") %>%
        pivot_longer(-Dataset_1, names_to = 'Dataset_2', values_to = x) %>%
        filter(!is.na(!!sym(x)))
      list.d[[x]] <- d
    }
    all.w <- full_join(list.d$pos, list.d$neg) %>%
      mutate(pos=case_when(is.na(pos)~1,TRUE~pos)) %>%
      mutate(neg=case_when(is.na(neg)~1,TRUE~neg))
    d.all <- rowMeans(all.w %>% select(pos, neg))
    all.w$all <- d.all
    dist.mat <- all.w %>%
      select(Dataset_1, Dataset_2, all) %>%
      pivot_wider(names_from = Dataset_2, values_from = all) %>%
      as.data.frame
    rownames(dist.mat) <- dist.mat$Dataset_1
    dist.mat$Dataset_1 <- NULL
    dist.mat <- as.matrix(dist.mat)
  }
  
  pco.res <- pco(dist.mat)
  var.explained <- sprintf(fmt='%.2f', 
                           pco.res$eig[1:2]/sum(pco.res$eig)*100)
  if (class(indiv) == 'logical' && indiv){
    df.plot <- pco.res$points %>% 
      as_tibble(rownames = 'ID', .name_repair=make.names) %>% 
      separate(ID, into = c('dataset.id', 'case', 'type', 
                            'ml.method', 'no'), sep = '-') %>% 
      mutate(group=paste0(dataset.id, '-', case))
    
  } else {
    df.plot <- pco.res$points %>% 
      as_tibble(rownames = 'ID', .name_repair=make.names)
    if (id=='-'){
      df.plot <- df.plot %>% 
        separate(ID, into = c('dataset.id', 'case'), sep = '-') %>% 
        mutate(group=paste0(dataset.id, '-', case))
    } else {
      df.plot <- df.plot %>% 
        separate(ID, into = c('dataset.id', 'case', 'type', 
                              'ml.method'), sep = '-') %>% 
        mutate(group=paste0(dataset.id, '-', case))
    }
  }
  g <- df.plot %>% 
    ggplot(aes(x=X, y=X.1, col=case)) + 
    geom_point() +
    scale_colour_manual(values=disease_colours) +
    theme_publication() + 
    xlab(paste0('PCo1 [', var.explained[1],'%]')) +
    ylab(paste0('PCo2 [', var.explained[2],'%]')) +
    NULL
  if (indiv == 'cross'){
    g <- df.plot %>% 
      group_by(group, dataset.id, case) %>% 
      summarise(x=mean(X), y=mean(X.1), s.x=sd(X), s.y=sd(X.1)) %>% 
      ggplot(aes(x=x, y=y)) + 
        geom_errorbarh(aes(col=case, y=y, xmin=x-s.x, xmax=x+s.x)) +
        geom_errorbar(aes(col=case, x=x, ymin=y-s.y, ymax=y+s.y)) +
        geom_point(aes(fill=case),shape=23) + 
        scale_fill_manual(values=disease_colours) +
        scale_colour_manual(values=disease_colours) +
        theme_publication() + 
        geom_text_repel(aes(label=dataset.id, col=case)) + 
        NULL
  } else if (indiv=='90percent'){
    g <- df.plot %>% 
      group_by(group, dataset.id, case) %>% 
      summarise(x=mean(X), y=mean(X.1), 
                s.x.mx=quantile(X, probs = 0.95), 
                s.x.mn=quantile(X, probs = 0.05), 
                s.y.mx=quantile(X.1, probs = 0.95), 
                s.y.mn=quantile(X.1, probs = 0.05)) %>% 
      ggplot(aes(x=x, y=y)) + 
        geom_errorbarh(aes(col=case, y=y, xmin=s.x.mn, xmax=s.x.mx)) +
        geom_errorbar(aes(col=case, x=x, ymin=s.y.mn, ymax=s.y.mx)) +
        geom_point(aes(fill=case),shape=23) + 
        scale_fill_manual(values=disease_colours) +
        scale_colour_manual(values=disease_colours) +
        theme_publication() + 
        geom_text_repel(aes(label=dataset.id, col=case)) + 
        NULL
  } else if (indiv){
    g <- g + stat_ellipse(aes(group=group), col='#00000060')
  } else {
    g <- g + geom_text_repel(aes(label=dataset.id))
  }
  return(g)
}

# of all models
pdf(here('figures', 'weights', paste0('pco_weights_', ml.method, '.pdf')), 
    width = 6, height = 5, useDingbats = FALSE)
for (distance in c('canberra')){
  g <- weight.pco(weights.mat.naive.rel, ml.method, distance)
  print(g + ggtitle(paste0('Distance: ', distance, ', ML: ', 
                           ml.method, ', Type: Naive')))
  g <- weight.pco(weights.mat.augm.rel, ml.method, distance)
  print(g + ggtitle(paste0('Distance: ', distance, ', ML: ', 
                           ml.method, ', Type: Augm')))
  g <- weight.pco(weights.mat.naive.rel, ml.method, distance, indiv='cross')
  print(g + ggtitle(paste0('Distance: ', distance, ', ML: ', 
                           ml.method, ', Type: Naive')))
  g <- weight.pco(weights.mat.augm.rel, ml.method, distance, indiv='cross')
  print(g + ggtitle(paste0('Distance: ', distance, ', ML: ', 
                           ml.method, ', Type: Augm')))
  g <- weight.pco(weights.mat.naive.rel, ml.method, distance, indiv='90percent')
  print(g + ggtitle(paste0('Distance: ', distance, ', ML: ', 
                           ml.method, ', Type: Naive')))
  g <- weight.pco(weights.mat.augm.rel, ml.method, distance, indiv='90percent')
  print(g + ggtitle(paste0('Distance: ', distance, ', ML: ', 
                           ml.method, ', Type: Augm')))
  g <- weight.pco(weights.rel.naive.ind, ml.method, distance, indiv = FALSE)
  print(g + ggtitle(paste0('Distance: ', distance, ', ML: ', 
                           ml.method, ', Type: Naive')))
  g <- weight.pco(weights.rel.augm.ind, ml.method, distance, indiv=FALSE)
  print(g + ggtitle(paste0('Distance: ', distance, ', ML: ', 
                           ml.method, ', Type: Augm')))
}
dev.off()


# ##############################################################################
# PCoA of associations
pdf(here('figures', 'weights', 'pco_assoc.pdf'), 
    width = 6, height = 5, useDingbats = FALSE)
for (distance in c('canberra')){
  # motus level
  g <- weight.pco(assoc.mat, id='-', dist.method = distance, indiv = FALSE)
  print(g + ggtitle(paste0("Distance: ", distance, 
                           ", Level: motus, Type: filtered")))
  g <- weight.pco(assoc.mat.full, id='-', dist.method = distance, indiv = FALSE)
  print(g + ggtitle(paste0("Distance: ", distance, 
                           ", Level: motus, Type: all")))
  # genus level
  g <- weight.pco(assoc.mat.g, id='-', dist.method = distance, indiv = FALSE)
  print(g + ggtitle(paste0("Distance: ", distance, 
                           ", Level: genus, Type: filtered")))
  g <- weight.pco(assoc.mat.g.full, id='-', dist.method = distance, 
                  indiv = FALSE)
  print(g + ggtitle(paste0("Distance: ", distance, 
                           ", Level: genus, Type: all")))
}
dev.off()

# ##############################################################################
# compute importance of genera for weights
col.order <- c('Wen_2017-AS',
               'Zhang_2015-RA',
               'Zeller_2014-CRC','Yu_2017-CRC', 'Wirbel_2019-CRC', 
               #'Vogtmann_2016-CRC',
               'Thomas_2019-CRC','Feng_2015-CRC',
               #'Zeller_2014-ADA', 'Feng_2015-ADA',
               'Thomas_2019-ADA',
               'Jie_2017-ACVD',
               'metaHIT-CD', 'Lewis_2015-CD', 'HMP2-CD', 'He_2017-CD','Franzosa_2019-CD',
               'metaHIT-UC','HMP2-UC','Franzosa_2019-UC',
               'Qin_2014-LIV',
               #'Loomba_2017-NAFLD','Hoyles_2018-NAFLD',
               'Bedarf_2017-PD')


sel.cols <- which(str_detect(colnames(weights.rel.augm.ind), 'enet0.5'))
disease.weights <- weights.rel.augm.ind %>%
  select(c(mOTU_ID, sel.cols)) %>% 
  pivot_longer(cols=-mOTU_ID) %>%
  mutate(name=str_remove(name, '-(augm5|original)-.*')) %>% 
  filter(name%in%col.order) %>%
  separate(col=name, into=c('dataset', 'case'), sep='-') %>% 
  group_by(mOTU_ID, case) %>% 
  summarise(x=mean(value)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = case, values_from = x) 
sel.motus.mat <- disease.weights %>% 
  mutate_if(is_double, .funs = list(function(x){x/sum(abs(x))})) %>% 
  pivot_longer(-mOTU_ID, names_to = 'case', values_to = 'weight') %>% 
  group_by(case) %>%
  arrange(desc(abs(weight))) %>% 
  mutate(cum.size=cumsum(abs(weight))) %>% 
  filter(cum.size < 0.5) %>% 
  ungroup() %>% 
  select(mOTU_ID, case, weight) %>% 
  pivot_wider(names_from = case, values_from = weight, 
              values_fill=list("weight"=0)) %>% 
  as.data.frame()
rownames(sel.motus.mat) <- sel.motus.mat$mOTU_ID
sel.motus.mat$mOTU_ID <- NULL
sel.motus.mat <- as.matrix(sel.motus.mat)

df.plot <- enframe(rowSums(sel.motus.mat > 0), name='mOTU', value='No_pos') %>% 
  full_join(enframe(rowSums(sel.motus.mat < 0), name='mOTU', value='No_neg'),
            by='mOTU') %>% 
  mutate(max.value=rowMaxs(sel.motus.mat)) %>% 
  mutate(min.value=rowMins(sel.motus.mat)) %>% 
  mutate(mOTU_ID=str_extract(mOTU, '(ref|meta)_mOTU_v25_[0-9]{5}')) %>% 
  left_join(motus.taxonomy %>% select(mOTU_ID, genus))

temp <- df.plot %>% 
  mutate(No_neg=-No_neg) %>% 
  mutate(n.rev=case_when(abs(No_neg) > No_pos~No_neg, 
                         No_pos > abs(No_neg)~No_pos,
                         abs(min.value) > max.value~No_neg,
                         TRUE~No_pos)) %>% 
  mutate(weight=case_when(abs(No_neg) > No_pos~min.value, 
                          No_pos > abs(No_neg)~max.value,
                          abs(min.value) > max.value~min.value,
                          TRUE~max.value)) %>% 
  arrange(n.rev, weight)

temp.pos <- df.plot %>% 
  group_by(genus) %>% 
  summarise(mn=min(min.value), mx=max(max.value)) %>% 
  filter(!str_detect(genus, '^NA')) %>% 
  filter(mx > 0.015) %>% arrange(mx)
temp.neg <- df.plot %>% 
  group_by(genus) %>% 
  summarise(mn=min(min.value), mx=max(max.value)) %>% 
  filter(!str_detect(genus, '^NA')) %>% 
  filter(mn < -0.015) %>% arrange(mn)

x.neg <- df.plot %>% 
  filter(genus %in% temp.neg$genus) %>% 
  group_by(genus) %>% 
  summarise(n=n(), n.dir=sum(No_neg),
            mx=min(min.value)) %>% 
  arrange(desc(n.dir), mx) %>% 
  mutate(genus=factor(genus, levels = genus))
x.pos <- df.plot %>% 
  filter(genus %in% temp.pos$genus) %>% 
  group_by(genus) %>% 
  summarise(n=n(), n.dir=sum(No_pos),
            mx=max(max.value)) %>% 
  arrange(desc(n.dir), desc(mx)) %>% 
  mutate(genus=factor(genus, levels = rev(genus)))

n.max <- max(x.pos$n, x.neg$n)

g.2.neg <- sel.motus.mat %>% 
  as_tibble(rownames = 'mOTU') %>% 
  pivot_longer(-mOTU, names_to = 'disease', values_to = 'weight') %>% 
  filter(mOTU %in% (df.plot %>% 
                      filter(genus %in% temp.neg$genus) %>% 
                      pull(mOTU))) %>% 
  mutate(mOTU_ID=str_extract(mOTU, '(ref|meta)_mOTU_v25_[0-9]{5}')) %>% 
  left_join(motus.taxonomy %>% select(mOTU_ID, genus)) %>% 
  mutate(genus=factor(genus, levels = levels(x.neg$genus))) %>% 
  mutate(weight=case_when((weight > 0)~0, TRUE~-weight)) %>% 
  filter(weight> 0) %>% 
  arrange(weight) %>% 
  ggplot(aes(x=genus, y=weight, fill=disease)) + 
    geom_jitter(width = 0.08, pch=23) + 
    xlab('') + ylab('Weight') + 
    theme_publication() + 
    theme(axis.text.x = element_blank(), axis.ticks.x=element_blank())  +
    scale_fill_manual(values=disease_colours, guide=FALSE) + 
    ylim(0, 0.061) +
    NULL
g.2.pos <- sel.motus.mat %>% 
  as_tibble(rownames = 'mOTU') %>% 
  pivot_longer(-mOTU, names_to = 'disease', values_to = 'weight') %>% 
  filter(mOTU %in% (df.plot %>% 
                      filter(genus %in% temp.pos$genus) %>% 
                      pull(mOTU))) %>% 
  mutate(mOTU_ID=str_extract(mOTU, '(ref|meta)_mOTU_v25_[0-9]{5}')) %>% 
  left_join(motus.taxonomy %>% select(mOTU_ID, genus)) %>% 
  mutate(genus=factor(genus, levels = levels(x.pos$genus))) %>% 
  mutate(weight=case_when((weight < 0)~0, TRUE~weight)) %>% 
  filter(weight> 0) %>% 
  arrange(weight) %>% 
  ggplot(aes(x=genus, y=weight, fill=disease)) + 
    geom_jitter(width = 0.08, pch=23, colour='black') + 
    xlab('') + ylab('Weight') + 
    theme_publication() + 
    theme(axis.text.x = element_blank(), axis.ticks.x=element_blank())  +
    scale_fill_manual(values=disease_colours, guide=FALSE) + 
    ylim(0, 0.061) +
    NULL
g.4.neg <- sel.motus.mat %>% 
  as_tibble(rownames = 'mOTU') %>% 
  pivot_longer(-mOTU, names_to = 'disease', values_to = 'weight') %>% 
  filter(mOTU %in% (df.plot %>% 
                      filter(genus %in% temp.neg$genus) %>% 
                      pull(mOTU))) %>% 
  mutate(mOTU_ID=str_extract(mOTU, '(ref|meta)_mOTU_v25_[0-9]{5}')) %>% 
  left_join(motus.taxonomy %>% select(mOTU_ID, genus)) %>% 
  mutate(genus=factor(genus, levels = levels(x.neg$genus))) %>% 
  mutate(weight=case_when((weight > 0)~0, TRUE~weight)) %>% 
  filter(weight < 0) %>% 
  ggplot(aes(x=genus, fill=disease)) + 
  geom_bar(stat="count") + 
  theme_publication() +
  scale_fill_manual(values=disease_colours, guide=FALSE) + 
  ylim(0, 15) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=1),
        panel.grid.major.y = element_line(colour='lightgrey', linetype=3)) + 
  NULL
g.4.pos <- sel.motus.mat %>% 
  as_tibble(rownames = 'mOTU') %>% 
  pivot_longer(-mOTU, names_to = 'disease', values_to = 'weight') %>% 
  filter(mOTU %in% (df.plot %>% 
                      filter(genus %in% temp.pos$genus) %>% 
                      pull(mOTU))) %>% 
  mutate(mOTU_ID=str_extract(mOTU, '(ref|meta)_mOTU_v25_[0-9]{5}')) %>% 
  left_join(motus.taxonomy %>% select(mOTU_ID, genus)) %>% 
  mutate(genus=factor(genus, levels = levels(x.pos$genus))) %>% 
  mutate(weight=case_when((weight < 0)~0, TRUE~weight)) %>% 
  filter(weight> 0) %>% 
  ggplot(aes(x=genus, fill=disease)) + 
    geom_bar(stat="count") + 
    theme_publication() +
    scale_fill_manual(values=disease_colours, guide=FALSE) + 
    ylim(0, 15) +
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=1),
          panel.grid.major.y = element_line(colour='lightgrey', linetype=3)) + 
    NULL
g.pos <- plot_grid(g.2.pos, g.4.pos, ncol = 1,
                   rel_heights = c(0.4, 0.6))
g.neg <- plot_grid(g.2.neg, g.4.neg, ncol = 1,
                   rel_heights = c(0.4, 0.6))
ggsave(g.pos, filename = here('figures', 'weights', 'genus_plot_pos.pdf'),
       width = 4, height = 5, useDingbats=FALSE)
ggsave(g.neg, filename = here('figures', 'weights', 'genus_plot_neg.pdf'),
       width = 6, height = 5, useDingbats=FALSE)

# heatmap of odds ratios for the selected genera
# negative first
prep <- sel.motus.mat %>% 
  as_tibble(rownames = 'mOTU') %>% 
  pivot_longer(-mOTU, names_to = 'disease', values_to = 'weight') %>% 
  filter(mOTU %in% (df.plot %>% 
                      filter(genus %in% as.character(x.neg$genus)) %>% 
                      pull(mOTU))) %>% 
  mutate(mOTU_ID=str_extract(mOTU, '(ref|meta)_mOTU_v25_[0-9]{5}')) %>% 
  left_join(motus.taxonomy %>% select(mOTU_ID, genus)) %>% 
  mutate(weight=as.numeric(weight!=0)) %>% 
  filter(weight != 0)
test <- vapply(unique(prep$disease), FUN=function(x){
  vapply(unique(prep$genus), FUN=function(y){
    temp <- prep %>% 
      mutate(disease=case_when(disease==x~x, TRUE~'other')) %>% 
      mutate(genus=case_when(genus==y~y, TRUE~'other'))
    t <- fisher.test(table(temp$disease, temp$genus))
    if (is.infinite(t$estimate)){
      if (table(temp$disease, temp$genus)[x,y] > 0){
        return(Inf)
      } else {
        return(NA_real_)
      }
    } else if (t$estimate==0){
      return(NA_real_)
    } else {
      return(t$estimate)
    }
    
  }, FUN.VALUE = double(1))}, 
  FUN.VALUE=double(length(unique(prep$genus))))
test %>% 
  as_tibble(rownames = 'genus') %>% 
  pivot_longer(-genus, values_to = 'odds_ratio', names_to = 'disease') %>% 
  mutate(odds_ratio=case_when(is.infinite(odds_ratio)~13,
                              TRUE~odds_ratio)) %>% 
  mutate(genus=factor(genus, levels = levels(x.neg$genus))) %>% 
  mutate(disease=factor(disease, levels=names(disease_colours))) %>% 
  ggplot(aes(x=genus, y=disease, fill=odds_ratio)) + 
    geom_tile() + 
    scale_fill_continuous_embl(na.value = 'white') +
    theme_publication() +
    theme(axis.text.x=element_text(angle=90, hjust=1))

# ##############################################################################
# Viz genera heatmaps

# function to plot weight heatmaps for different genera
heatmap.genus <- function(mat, id, g){
  stopifnot(id %in% c('lasso', 'enet$', 'enet0.5'))
  stopifnot('mOTU_ID' %in% colnames(mat))

  sel.cols <- which(str_detect(colnames(mat), id))
  mat.red <- mat %>% 
    select(c(mOTU_ID, sel.cols))
  
  # select motus that make up more than 50% of relative weight
  sel.motus <- mat.red %>%
    pivot_longer(cols=-mOTU_ID) %>%
    mutate(name=str_remove(name, '-(augm5|original)-.*')) %>% 
    filter(name%in%col.order) %>%
    group_by(name) %>%
    arrange(desc(abs(value))) %>%
    mutate(cum.size=cumsum(abs(value))) %>%
    filter(cum.size < 0.5) %>%
    ungroup()
  mat.red.filt <- mat.red %>% 
    filter(mOTU_ID %in% sel.motus$mOTU_ID) %>% 
    mutate(mOTU_ID=str_extract(mOTU_ID, '(ref|meta)_mOTU_v25_[0-9]{5}')) %>% 
    left_join(motus.taxonomy %>% select(mOTU_ID, genus), by='mOTU_ID') %>% 
    filter(str_detect(genus, g)) %>% 
    as.data.frame()
  mat.red.filt$genus <- NULL
  rownames(mat.red.filt) <- mat.red.filt$mOTU_ID
  mat.red.filt$mOTU_ID <- NULL
  
  colnames(mat.red.filt) <- str_remove(colnames(mat.red.filt), 
                                       paste0('-(augm5|original)-', id))
  # rename rows with nice names
  motus.names <- sel.motus %>% 
    transmute(full.name=mOTU_ID, 
              mOTU_ID=str_extract(mOTU_ID, '(ref|meta)_mOTU_v25_[0-9]{5}')) %>% 
    mutate(nice.name=str_remove(full.name, 
                                ' \\[(ref|meta)_mOTU_v25_[0-9]{5}\\]')) %>% 
    mutate(nice.name=paste0(nice.name, '-', str_extract(mOTU_ID, '[0-9]{5}'))) %>% 
    distinct() %>% 
    left_join(motus.taxonomy %>% select(mOTU_ID, genus), by='mOTU_ID') %>% 
    filter(str_detect(genus, g))
  
  
  # re-order cols according to taxonomy
  sub.tree <- keep.tip(motus.tree, motus.names$mOTU_ID)
  
  
  mat.red.filt <- -mat.red.filt[sub.tree$tip.label,
                                intersect(col.order, colnames(mat.red.filt))]
  rownames(mat.red.filt) <- motus.names$nice.name[match(rownames(mat.red.filt), 
                                                        motus.names$mOTU_ID)]
  # temp <- hclust(vegdist(mat.red.filt, method = 'manhattan'),method = 'ward.D2')
  lim <- 0.07# max(abs(mat.red.filt))
  g <- pheatmap(mat.red.filt, cluster_cols = FALSE, cluster_rows = FALSE,
                color=colorRampPalette(
                  embl.palette.data$diverging$`Blue-Red`$value)(100),
                breaks = seq(-lim, lim, length.out = 100))
  return(g)
}

pdf(here('figures', 'weights', 'genus_heatmaps.pdf'), 
    useDingbats = FALSE, width = 8, height = 6)
for (g in unique(c(as.character(x.pos$genus), 
                   as.character(x.neg$genus)))){
  print(heatmap.genus(weights.rel.augm.ind, 'enet0.5', g)) 
}
dev.off()

print(heatmap.genus(weights.rel.augm.ind, 'enet0.5', 'Fuso'))
print(heatmap.genus(weights.rel.augm.ind, 'enet0.5', 'Lactobacillus'))
print(heatmap.genus(weights.rel.augm.ind, 'enet0.5', 'Bacteroides'))
print(heatmap.genus(weights.rel.augm.ind, 'enet0.5', 'Solobacterium'))
print(heatmap.genus(weights.rel.augm.ind, 'enet0.5', 'Veillonella'))
print(heatmap.genus(weights.rel.augm.ind, 'enet0.5', 'Neisseria'))
print(heatmap.genus(weights.rel.augm.ind, 'enet0.5', 'Actinomyces'))
print(heatmap.genus(weights.rel.augm.ind, 'enet0.5', 'Oscillibacter'))
print(heatmap.genus(weights.rel.augm.ind, 'enet0.5', 'Gemella'))
print(heatmap.genus(weights.rel.augm.ind, 'enet0.5', 'Staphylococcus'))
print(heatmap.genus(weights.rel.augm.ind, 'enet0.5', 'Romboutsia'))
print(heatmap.genus(weights.rel.augm.ind, 'enet0.5', 'Anaerostipes'))
print(heatmap.genus(weights.rel.augm.ind, 'enet0.5', 'Bifidobacterium'))
print(heatmap.genus(weights.rel.augm.ind, 'enet0.5', 'Eubacterium'))
dev.off()

# # ##############################################################################
# motus association heatmaps as well

heatmap.assoc.genus <- function(mat, g){
  temp <- mat %>% 
    pivot_longer(-mOTU_ID, values_to = 'gFC', names_to = 'dataset') %>% 
    mutate(mOTU_name=mOTU_ID) %>% 
    mutate(mOTU_ID=str_extract(mOTU_name, '(ref|meta)_mOTU_v25_[0-9]{5}')) %>% 
    left_join(motus.taxonomy %>% select(mOTU_ID, genus)) %>% 
    filter(str_detect(genus, g)) %>% 
    select(mOTU_name, gFC, dataset) %>% 
    filter(dataset%in%col.order) %>% 
    pivot_wider(values_from = gFC, names_from = dataset) %>% 
    as.data.frame
  rownames(temp) <- temp$mOTU_name
  temp$mOTU_name <- NULL
  temp <- as.matrix(temp)
  
  motus.names <- enframe(rownames(temp), value = 'full.name') %>% 
    mutate(mOTU_ID=str_extract(full.name,
                               '(ref|meta)_mOTU_v25_[0-9]{5}')) %>% 
    mutate(nice.name=str_remove(full.name, 
                                ' \\[(ref|meta)_mOTU_v25_[0-9]{5}\\]')) %>% 
    mutate(nice.name=paste0(nice.name, '-', 
                            str_extract(mOTU_ID, '[0-9]{5}'))) %>% 
    distinct() %>% 
    left_join(motus.taxonomy %>% select(mOTU_ID, genus), by='mOTU_ID') %>% 
    filter(str_detect(genus, g))
  sub.tree <- keep.tip(motus.tree, motus.names$mOTU_ID)
  
  rownames(temp) <- str_extract(rownames(temp), '(ref|meta)_mOTU_v25_[0-9]{5}')
  temp <- temp[sub.tree$tip.label,
                                intersect(col.order, colnames(temp))]
  rownames(temp) <- motus.names$nice.name[match(rownames(temp), 
                                                motus.names$mOTU_ID)]
  # temp <- hclust(vegdist(mat.red.filt, method = 'manhattan'),method = 'ward.D2')
  lim <- 1.5# max(abs(mat.red.filt))
  g <- pheatmap(temp, cluster_cols = FALSE, cluster_rows = FALSE,
                color=colorRampPalette(
                  embl.palette.data$diverging$`Blue-Red`$value)(100),
                breaks = seq(-lim, lim, length.out = 100))
}

pdf('~/Documents/google_drive/SIAMCAT/WIP/assoc_heatmaps.pdf',
    useDingbats = FALSE, width = 8, height = 6)
heatmap.assoc.genus(assoc.mat, 'Fuso')
heatmap.assoc.genus(assoc.mat, 'Rumi')
heatmap.assoc.genus(assoc.mat, 'Lachno')
heatmap.assoc.genus(assoc.mat, 'Blautia')
heatmap.assoc.genus(assoc.mat, 'Eubac')
heatmap.assoc.genus(assoc.mat, 'Lact')
heatmap.assoc.genus(assoc.mat, 'Bifido')
dev.off()

# # ##############################################################################
# correlation between associations and weights




df.sim <- full_join(weights.rel.augm.ind %>%
                      pivot_longer(-mOTU_ID, names_to = 'Model',
                                   values_to='Weight') %>% 
                      separate(col = Model,
                               into=c('dataset.id', 'case', 'type', 
                                      'model.type'),
                               sep='-') %>% 
                      mutate(Dataset=paste0(dataset.id, '-', case)) %>% 
                      filter(model.type=='enet0.5') %>%
                      transmute(Dataset, w.augm=Weight, mOTU_ID),
                    weights.rel.naive.ind %>%
                      pivot_longer(-mOTU_ID, names_to = 'Model',
                                   values_to='Weight') %>% 
                      separate(col = Model,
                               into=c('dataset.id', 'case', 'type', 
                                      'model.type'),
                               sep='-') %>% 
                      mutate(Dataset=paste0(dataset.id, '-', case)) %>% 
                      filter(model.type=='enet0.5') %>%
                      transmute(Dataset, w.naive=Weight, mOTU_ID)) %>% 
  full_join(assoc.mat.full %>% 
              pivot_longer(-mOTU_ID, names_to = 'Dataset', 
                           values_to = 'gFC'))
df.sim %>%
  group_by(Dataset) %>%
  summarise(c.w=cor(w.naive, w.augm, use='pairwise.complete.obs',
                  method='spearman'),
            c.n=cor(w.naive, -gFC, use='pairwise.complete.obs',
                    method='spearman'),
            c.a=cor(w.augm, -gFC, use='pairwise.complete.obs',
                    method='spearman')) %>% 
  pivot_longer(-Dataset, names_to = 'type', values_to='corr') %>% 
  ggplot(aes(x=type, y=corr)) +
    geom_boxplot()

g <- df.sim %>% 
  filter(Dataset %in% col.order) %>% 
  filter(!(w.naive==0 & w.augm == 0)) %>% 
  ggplot(aes(x=w.naive, y=w.augm)) + 
    geom_point() + 
    facet_wrap(~Dataset, ncol=4, nrow = 5,
                scales = 'free')
ggsave(g, filename = '~/Desktop/weights.pdf', width = 12, height = 12)
