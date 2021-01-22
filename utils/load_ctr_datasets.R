# ##############################################################################
#
## Load control datasets for the control-augmentation
#
# ##############################################################################

.f_load_control_data <- function(type){
  datasets <- NULL
  stopifnot(type %in% c('cohort', 'random', 'other_ctr'))
  
  if (type == 'cohort'){
    
    ## Zeevi
    # metadata
    meta.zeevi <- read_tsv(here("data", 'meta', 'meta_Zeevi.tsv')) %>% 
      group_by(Individual_ID) %>% 
      filter(Timepoint==min(Timepoint)) %>% 
      ungroup()
    # features
    feat.zeevi <- read.table(
      here('data', 'features', 'motus', 'Zeevi_2014_motus.tsv'),
      stringsAsFactors = FALSE, check.names = FALSE, sep='\t',
      row.names = 1, header = TRUE, comment.char = '', quote = '')
    feat.zeevi <- feat.zeevi[,colSums(feat.zeevi) > 100]
    feat.zeevi <- feat.zeevi[,intersect(colnames(feat.zeevi), 
                                        meta.zeevi$Sample_ID)]
    feat.zeevi <- prop.table(as.matrix(feat.zeevi), 2)
    message('Loaded Zeevi data!')
    ## TwinsUK
    feat.twinsuk <- read.table(
      here('data', 'features', 'motus', 'TwinsUK_motus.tsv'),
      stringsAsFactors = FALSE, check.names = FALSE, sep='\t',
      row.names = 1, header = TRUE, comment.char = '', quote = '')
    feat.twinsuk <- feat.twinsuk[,colSums(feat.twinsuk) > 100]
    feat.twinsuk <- prop.table(as.matrix(feat.twinsuk), 2)
    message('Loaded TwinsUK data!')
    ## Schirmer
    feat.schirmer <- read.table(
      here('data', 'features', 'motus', 'Schirmer_2016_motus.tsv'),
      stringsAsFactors = FALSE, check.names = FALSE, sep='\t',
      row.names = 1, header = TRUE, comment.char = '', quote = '')
    feat.schirmer <- feat.schirmer[,colSums(feat.schirmer) > 100]
    feat.schirmer <- prop.table(as.matrix(feat.schirmer), 2)
    message('Loaded Schirmer data!')
    feat.ctr <- list(feat.zeevi, 
                     feat.twinsuk, 
                     feat.schirmer)
  } else if (type == 'other_ctr'){
    # metaHIT DNK samples
    meta.metaHIT <- read_tsv(here('data', 'meta', 'meta_metaHIT.tsv'))
    meta.metaHIT <- meta.metaHIT %>% filter(Country=='danish')
    feat.metaHIT <- read.table(here('data', 'features', 'motus', 
                                    'metaHIT_motus.tsv'),
                               sep='\t', stringsAsFactors = FALSE, 
                               check.names = FALSE, quote = '', 
                               comment.char = '')
    feat.metaHIT <- feat.metaHIT[,meta.metaHIT$Sample_ID]
    feat.metaHIT <- feat.metaHIT[,colSums(feat.metaHIT) > 100]
    feat.metaHIT <- prop.table(as.matrix(feat.metaHIT), 2)
    message('Loaded metaHIT data!')
    # Backhed mothers
    meta.backhed <- read_tsv(here('data', 'meta', 'meta_Backhed_2015.tsv')) %>% 
      filter(Description=='Mother')
    feat.backhed <- read.table(here('data', 'features', 
                                    'motus', 'Backhed_2015_motus.tsv'),
      stringsAsFactors = FALSE, check.names = FALSE, sep='\t',
      row.names = 1, header = TRUE, comment.char = '', quote = '')
    feat.backhed <- feat.backhed[,colSums(feat.backhed) > 100]
    feat.backhed <- feat.backhed[,intersect(colnames(feat.backhed), 
                                            meta.backhed$Sample_ID)]
    feat.backhed <- prop.table(as.matrix(feat.backhed), 2)
    message('Loaded Backhed data!')
    # Poyet
    meta.poyet <- read_tsv(here('data', 'meta', 'meta_Poyet_2019.tsv')) %>% 
      group_by(Individual_ID) %>% 
      filter(Timepoint==min(Timepoint)) %>% 
      ungroup()
    feat.poyet <- read.table(here('data', 'features', 
                                  'motus', 'Poyet_2019_motus.tsv'),
      stringsAsFactors = FALSE, check.names = FALSE, sep='\t',
      row.names = 1, header = TRUE, comment.char = '', quote = '')
    feat.poyet <- feat.poyet[,colSums(feat.poyet) > 100]
    feat.poyet <- feat.poyet[,intersect(colnames(feat.poyet), 
                                        meta.poyet$Sample_ID)]
    feat.poyet <- prop.table(as.matrix(feat.poyet), 2)
    message('Loaded Poyet data!')
    # Zhu
    meta.zhu <- read_tsv(here('data', 'meta', 'meta_Zhu_2018.tsv')) %>% 
      filter(Group=='CTR')
    feat.zhu <- read.table(here('data', 'features', 
                                'motus', 'Zhu_2018_motus.tsv'),
      stringsAsFactors = FALSE, check.names = FALSE,  sep='\t',
      row.names = 1, header = TRUE, comment.char = '', quote = '')
    feat.zhu <- feat.zhu[,colSums(feat.zhu) > 100]
    feat.zhu <- feat.zhu[,intersect(colnames(feat.zhu), 
                                    str_remove(meta.zhu$Sample_ID, '[DB]+$'))]
    feat.zhu <- prop.table(as.matrix(feat.zhu), 2)
    message('Loaded Zhu data!')
    # Vincent
    meta.vincent <- read_tsv(here('data', 'meta', 'meta_Vincent_2016.tsv')) %>% 
      filter(Group=='CTR') %>% 
      group_by(Individual_ID) %>% 
      filter(Timepoint==min(Timepoint)) %>% 
      ungroup()
    feat.vincent <- read.table(here('data', 'features', 
                                    'motus', 'Vincent_2016_motus.tsv'),
      stringsAsFactors = FALSE, check.names = FALSE,  sep='\t',
      row.names = 1, header = TRUE, comment.char = '', quote = '')
    feat.vincent <- feat.vincent[,colSums(feat.vincent) > 100]
    feat.vincent <- feat.vincent[,intersect(colnames(feat.vincent), 
                                            meta.vincent$Sample_ID)]
    feat.vincent <- prop.table(as.matrix(feat.vincent), 2)
    message('Loaded Vincent data!')
    # combine
    feat.ctr <- list(feat.metaHIT, feat.poyet, feat.zhu, 
                         feat.backhed, feat.vincent)
  } else if (type == 'random'){
    datasets <- motu.tasks %>% 
      select(dataset.id) %>% 
      distinct() %>% 
      filter(dataset.id != opt$dataset) %>% 
      pull(dataset.id) %>% 
      sample(4)
    print(datasets)
    feat.ctr <- list()
    for (d in datasets){
      fn.sc <- list.files(here('parameter_space', 'sc'), pattern=d) %>% 
        grep(pattern='mOTUs', value = TRUE) %>% sample(1)
      load(here('parameter_space', 'sc', fn.sc))
      feat <- get.orig_feat.matrix(sc.obj)
      feat.ctr[[(length(feat.ctr) + 1)]] <- feat
    }
  }
  return(list('feat.ctr'=feat.ctr, 'datasets'=datasets))
}


