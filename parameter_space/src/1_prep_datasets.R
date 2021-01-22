# ##############################################################################
#
##  Script to prepare datasets for the parameter space exploration
#
# ##############################################################################

library("tidyverse")
library("SIAMCAT")
library("here")
library("matrixStats")
library("vegan")
library("ggembl")
library("ggpubr")
library("curatedMetagenomicData")
library("feather")

# ##############################################################################
# mOTUs2 and eggNOG datasets

datasets <- c("Bedarf_2017", "Feng_2015", "Franzosa_2019", "He_2017",
              "HMP2", "Hoyles_2018", "Jie_2017", "Karlsson_2013", 
              "Kushugolova_2018", "LeChatelier_2013", "Lewis_2015", 
              "Li_2017", "Loomba_2017", "metaHIT", "Qin_2012", "Qin_2014", 
              "Thomas_2019", "Vogtmann_2016", "Wen_2017", "Wirbel_2019", 
              "Yachida_2019", "Yu_2017", "Zeller_2014", "Zhang_2015")

data.info <- tibble(dataset.id=character(0), type=character(0),
                    case=character(0), n.ctr=integer(0), 
                    n.case=integer(0))
df.list <- list()
metadata.list <- list()
for (d in datasets){
  message(d)
  # metadata
  fn.meta <- here('data', 'meta', paste0('meta_', d, '.tsv'))
  meta <- read_tsv(fn.meta, col_types = cols()) %>% 
    mutate(Sample_ID=make.names(Sample_ID))
  
  # features
  fn.feat <- paste0(here('data', 'features', 'motus', paste0(d, '_motus.tsv')))
  feat.temp <- read.table(fn.feat, sep='\t', stringsAsFactors = FALSE,
                          check.names = FALSE, quote = '', comment.char = '')
  feat.temp <- as.matrix(feat.temp)
  colnames(feat.temp) <- make.names(colnames(feat.temp))
  # check that these are mOTUs 2.5 profiles, not mOTUs 2.0
  # and that the import worked well (weird taxa names)
  stopifnot(nrow(feat.temp) == 14213)
  
  # remove samples with very low sequencing depth
  feat.temp <- feat.temp[,colSums(feat.temp) > 100]
  feat.rel <- prop.table(feat.temp, 2)
  
  # get overlap in metadata and features
  # since it not always fits between metadata and features
  samples <- intersect(meta$Sample_ID, colnames(feat.temp))
  stopifnot(length(samples) > 0)
  message(length(samples))
  meta <- meta %>% 
    filter(Sample_ID %in% samples)
  feat.rel <- feat.rel[,samples]
  
  # for datasets with multiple samples per indiviudal, 
  #   compute distances and take the first sample for each indiviual
  if (d %in% c('HMP2', 'Lewis_2015', 'metaHIT', 
               'Kushugolova_2018', 'Zhang_2015')){
    if (d == 'metaHIT'){
      meta <- meta %>% 
        filter(Country!='danish') %>% 
        mutate(Timepoint=Sampling_day)
      feat.temp <- feat.temp[,meta$Sample_ID]
      feat.rel <- feat.rel[,meta$Sample_ID]
    }
    if (d == 'Lewis_2015'){
      meta <- meta %>% 
        mutate(Individual_ID=paste0('ID_', Subject)) %>% 
        mutate(Timepoint=Time)
    }
    
    feat.filt <- feat.temp[rowMaxs(feat.rel) > 1e-03,]
    
    dist.mat <- as.matrix(vegdist(t(feat.filt)))
    diag(dist.mat) <- NA
    dist.mat[lower.tri(dist.mat)] <- NA
    temp <- as_tibble(dist.mat) %>% 
      mutate(Sample=rownames(dist.mat)) %>% 
      pivot_longer(-Sample) %>% 
      filter(!is.na(value))
    df.plot <- left_join(temp, meta %>% 
                           transmute(Sample=Sample_ID, ID1=Individual_ID)) %>% 
      left_join(meta %>% transmute(name=Sample_ID, ID2=Individual_ID)) %>% 
      mutate(type=ID1==ID2) %>% 
      mutate(Study=d)
    df.list[[d]] <- df.plot
    
    meta <- meta %>% 
      group_by(Individual_ID) %>% 
      filter(Timepoint==min(Timepoint)) %>% 
      ungroup()
  }
  metadata.list[[d]] <- meta %>% 
    mutate(Study=d) %>% 
    select(Sample_ID, Group, Study)
  meta <- as.data.frame(meta)
  rownames(meta) <- meta$Sample_ID
  message(nrow(meta))
  
  
  fn.feat.eggNOG <- here('data', 'features', 'eggNOG', 
                         paste0('eggNOG_', d, '.tsv'))
  feat.eggNOG <- read.table(fn.feat.eggNOG, sep='\t', stringsAsFactors = FALSE,
                            check.names = FALSE, row.names = 1, header = TRUE)
  stopifnot(all(rownames(meta) %in% colnames(feat.eggNOG)))
  feat.eggNOG.rel <- prop.table(as.matrix(feat.eggNOG), 2)

  # loop through groups
  groups <- setdiff(unique(meta$Group), 'CTR')
  for (g in groups){
    # create a label with SIAMCAT
    label <- create.label(label='Group', case=g, 
                          control='CTR', meta=meta, verbose=0)
    
    # siamcat
    sc.obj <- siamcat(feat=feat.rel, label=label, meta=meta, verbose = 0)
    n.ctr <- sum(sc.obj@label$label==-1)
    n.case <- sum(sc.obj@label$label==1)
    data.info <- data.info %>%
      add_row(dataset.id=d, case=g, type='mOTUs2',
              n.ctr=n.ctr, n.case=n.case)
    sc.obj <- create.data.split(sc.obj, num.folds = 5, num.resample = 5,
                                verbose=0)
    # save object
    save(sc.obj, file = here('parameter_space/sc/',
                             paste0('sc_', d,'_', g, '_mOTUs2.RData')))
    
    # create an eggNOG sc object as well and save
    sc.obj.eggNOG <- siamcat(feat=feat.eggNOG.rel, 
                             label=label, meta=meta, verbose = 0)
    data_split(sc.obj.eggNOG) <- data_split(sc.obj)
    data.info <- data.info %>%
      add_row(dataset.id=d, case=g, type='eggNOG',
              n.ctr=n.ctr, n.case=n.case)
    sc.obj <- sc.obj.eggNOG
    save(sc.obj, file = here('parameter_space/sc/',
                             paste0('sc_', d,'_', g, '_eggNOG.RData')))
  }
}


# plot similarity between and across individuals
df.plot.all <- bind_rows(df.list)

g <- df.plot.all %>% 
  mutate(type=as.character(type)) %>% 
  mutate(type=case_when(type=='TRUE'~'Within subjects', 
                        TRUE~'Across subjects')) %>% 
  ggplot(aes(x=Study, y=value, fill=type)) + 
    geom_boxplot(outlier.alpha = 0.1, outlier.stroke = 0) + 
    theme_publication() + 
    xlab('') + ylab('Bray-Curtis dissimilarity') +
    scale_fill_manual(values=c('#70737285', '#D4164585'), 
                      labels=c('Across subjects', 'Within subjects'),
                      name='') + 
    theme(axis.text.x=element_text(angle=45, hjust=1)) + 
    stat_compare_means(method='wilcox.test') + 
    NULL
ggsave(g, filename = here('figures', 'parameter_space', 
                          'subject_similarity.pdf'), 
       useDingbats=FALSE, units = 'mm', 
       width = 180, height = 80)


# ##############################################################################
# curatedMetagenomics

all.metadata <- combined_metadata %>% 
  filter(body_site=='stool')
for (i in unique(all.metadata$dataset_name)){
  temp <- all.metadata %>%
    filter(dataset_name == i)
  info <- temp %>%
    group_by(study_condition) %>%
    summarise(n=n()) %>%
    filter(n > 20)
  if (nrow(info) == 1 | !('control' %in% info$study_condition)){
    next()
  }
  case <- setdiff(info$study_condition, 'control')
  
  if (i == 'NielsenHB_2014'){
    temp <- temp %>% 
      filter(country!='DNK')
    temp <- temp %>% 
      mutate(timepoint=str_extract(sampleID, '_[0-9]+$')) %>% 
      mutate(timepoint=str_remove(timepoint, '_')) %>% 
      mutate(timepoint=as.numeric(timepoint)) %>% 
      group_by(subjectID) %>% 
      filter(timepoint==min(timepoint)) %>% 
      ungroup
  } else if (i %in% c('Heitz-BuschartA_2016', 'KosticAD_2015', 
                      'VatanenT_2016', 'RaymondF_2016', 'VincentC_2016')){
    next()
  }
  
  meta.temp <- temp %>%
    filter(study_condition %in% info$study_condition) %>%
    select(sampleID, subjectID, study_condition, visit_number,
           days_from_first_collection)
  metadata.list[[i]] <- meta.temp %>% 
    transmute(Sample_ID=sampleID, Group=study_condition, Study=i)
  
  meta.temp <- as.data.frame(meta.temp)
  rownames(meta.temp) <- meta.temp$sampleID
  
  x <- paste0(i, '.metaphlan_bugs_list.stool')
  feat <- curatedMetagenomicData(x=x, dryrun=FALSE)
  feat <- feat[[x]]@assayData$exprs
  # clean up metaphlan data
  feat <- feat[grep(x=rownames(feat), pattern='s__'),]
  feat <- feat[grep(x=rownames(feat),pattern='t__', invert = TRUE),]
  stopifnot(all(colSums(feat) != 0))
  feat <- t(t(feat)/100)
  
  
  # get humann2 data
  x <- paste0(i, '.pathabundance_relab.stool')
  feat.humann <- curatedMetagenomicData(x=x, dryrun=FALSE)
  feat.humann <- feat.humann[[x]]@assayData$exprs
  # clean up humann2 data
  sel.feat <- c(1, which(str_detect(rownames(feat.humann), '\\|')))
  feat.red <- feat.humann[sel.feat,]
  if (any(colSums(feat.red) != 0)){
    feat.red <- prop.table(feat.red, 2)
    feat.red[is.na(feat.red)] <- 0
  } else {
    feat.red <- prop.table(feat.red, 2)
  }
  
  
  for (x in case){
    label <- create.label(label='study_condition', case=x, control='control',
                          meta=meta.temp)
    sc.obj <- siamcat(feat=feat, label=label, meta=meta.temp, verbose=0)
    n.ctr <- sum(sc.obj@label$label==-1)
    n.case <- sum(sc.obj@label$label==1)
    data.info <- data.info %>%
      add_row(dataset.id=i, case=x, type='metaphlan',
              n.ctr=n.ctr, n.case=n.case)
    sc.obj <- create.data.split(sc.obj, num.folds = 5, num.resample = 5, 
                                verbose=0)
    save(sc.obj, file = here("parameter_space", "sc", 
                             paste0('sc_', i, '_', x, '_metaphlan.RData')))
    
    # get humann2 datasets as well
    sc.obj.humann <- siamcat(feat=feat.red, label=label, 
                             meta=meta.temp, verbose=0)
    data_split(sc.obj.humann) <- data_split(sc.obj)
    sc.obj <- sc.obj.humann
    save(sc.obj, file = here("parameter_space", "sc", 
                             paste0('sc_', i, '_', x, '_humann2.RData')))
    data.info <- data.info %>%
      add_row(dataset.id=i, case=x, type='humann2',
              n.ctr=n.ctr, n.case=n.case)
  }
}

# ##############################################################################
# 16S from Duvallet et al.
data.location <- 'data/duvallet/'

feat.files <- list.files(here(data.location), pattern='otu_table')
meta.files <- list.files(here(data.location), pattern='metadata')

# create dataset info table
cases <- toupper(str_extract(feat.files, '^[a-z0-9]+'))
dataset.ids <- str_remove(str_remove(feat.files, '^[a-z0-9]+_'),
                          '.otu_table.clean.feather')

# prep files
for (id in seq_along(dataset.ids)){
  dataset.id <- dataset.ids[id]
  case <- tolower(cases[id])
  message(paste(c(case, dataset.id), collapse = '_'))
  
  # metadata
  meta.file <- paste0(case, '_',
                      dataset.id, '.metadata.clean.feather')
  meta <- read_feather(here(data.location, meta.file))
  Sample_ID <- meta %>% pull(1)
  meta <- meta %>%
    mutate(Sample_ID = Sample_ID) %>%
    mutate(DiseaseState = toupper(DiseaseState)) %>%
    select(Sample_ID, DiseaseState)
  metadata.list[[dataset.id]] <- meta %>% 
    transmute(Sample_ID=Sample_ID, Group=DiseaseState, Study=dataset.id)
  
  # label preprocessing
  if (dataset.id == 'scher'){
    meta <- meta %>%
      mutate(DiseaseState = ifelse(DiseaseState == 'RA', 'ART', DiseaseState))
  } else if (dataset.id == 'zhang'){
    meta <- meta %>%
      mutate(DiseaseState = ifelse(DiseaseState == 'CIRR', 'LIV',
                                   DiseaseState))
  } else if (dataset.id %in% c('willing', 'morgan')){
    meta <- meta %>%
      mutate(DiseaseState = ifelse(DiseaseState %in% c('CD', "UC"), 'IBD',
                                   DiseaseState))
  } else if (dataset.id == 'papa'){
    meta <- meta %>%
      mutate(DiseaseState = ifelse(DiseaseState %in% c('CD', "UC"), 'IBD',
                                   DiseaseState)) %>%
      mutate(DiseaseState = ifelse(DiseaseState == 'NONIBD', 'H',
                                   DiseaseState))
  } else if (dataset.id == 'gevers'){
    meta <- meta %>%
      mutate(DiseaseState = ifelse(DiseaseState == 'NONIBD', 'H', 'IBD'))
  }
  
  if (meta %>% filter(DiseaseState == 'H') %>% nrow < 10){
    print(dataset.id)
    next()
  }
  
  meta <- as.data.frame(meta)
  rownames(meta) <- meta$Sample_ID
  
  # features and preprocessing
  feat.file <- paste0(case, '_',
                      dataset.id, '.otu_table.clean.feather')
  feat <- read_feather(here(data.location, feat.file))
  column.names <- feat$index
  feat$index <- NULL
  feature.names <- colnames(feat)
  feat <- t(data.frame(feat))
  colnames(feat) <- column.names
  feat.rel <- prop.table(feat, 2)
  
  # aggregate on genus level
  label <- create.label(label='DiseaseState', case=toupper(case),
                        control='H', meta=meta)
  sc.obj <- siamcat(feat=feat.rel, label=label)
  sc.obj <- summarize.features(sc.obj, level='g__')
  sc.obj <- create.data.split(sc.obj, num.folds = 5, num.resample = 5)
  n.ctr <- sum(sc.obj@label$label == -1)
  n.case <- sum(sc.obj@label$label == 1)
  data.info <- data.info %>%
    add_row(dataset.id=dataset.id, type='RDP', 
            case=toupper(case),
            n.ctr=n.ctr, n.case=n.case)
  
  # save features
  save(sc.obj, file = here("parameter_space/sc",
                           paste0('sc_', dataset.id, '_',
                                  toupper(case), '_RDP.RData')))
}


# ##############################################################################
# save data info
write_tsv(data.info, path = here('parameter_space', 'data_info', 
                                 'all_tasks.tsv'))

# some basic (impressive?) info
meta.all <- bind_rows(metadata.list)
meta.all <- meta.all %>% distinct()
write_tsv(meta.all, path=here('parameter_space', 'data_info', 'meta_full.tsv'))
