# ##############################################################################
#
## Check transfer to other datasets for normal best models
#
# ##############################################################################

library("tidyverse")
library("SIAMCAT")
library("pROC")
library("here")
library("ggembl")
library("yaml")

colours <- yaml.load_file(here("parameter_space", "data_info", "colours.yaml"))
disease_colours <- unlist(colours$disease)
# ##############################################################################
# function to get cross-predictions

source(here('utils', 'cross_predictions.R'))

# ##############################################################################
# load data sets

motu.tasks <- read_tsv(here('parameter_space', 'files', 
                              'auroc_all.tsv')) %>% 
  filter(type=='mOTUs2')
for (ml.method in c('enet', 'enet-0.5', 'lasso', 'randomForest')){
  fn.cross.pred <- here('parameter_space', 'files',
                        paste0('cross_prediction_', ml.method, '.tsv'))

  if (!file.exists(fn.cross.pred)){
    # prepare cross-prediction matrix
    df.cross.predictions <- tibble(n.group=integer(0),
                                   n.pos=integer(0),
                                   frac=double(0),
                                   cp.auroc=double(0),
                                   dataset.id=character(0),
                                   Group=character(0),
                                   dataset.test=character(0),
                                   dataset.train=character(0),
                                   dataset.id.train=character(0),
                                   case.train=character(0),
                                   auroc=double(0),
                                   ext.auroc=double(0))
    
    for (i in seq_len(nrow(motu.tasks))){
      dataset.id <- motu.tasks$dataset.id[[i]]
      case <- motu.tasks$case[[i]]
      message(dataset.id, '-', case)
      fn.sc <- here('parameter_space', 'models',
                    paste0('sc_trained_', dataset.id, '_', case, '_', 
                           ml.method,'.RData'))
      if (!file.exists(fn.sc)){
        message('Problem for ', dataset.id, '-', case)
        next()
      }
      load(fn.sc)
      df.cross.predictions <- bind_rows(
        df.cross.predictions, 
        f.cross.prediction(sc.obj.train, dataset.id, case))
    }
    
    df.cross.predictions <- df.cross.predictions %>% 
      distinct()
    
    write_tsv(df.cross.predictions, 
              path = fn.cross.pred)
  } else {
    df.cross.predictions <- read_tsv(fn.cross.pred)
  }

  # ############################################################################
  # plot batch effects and cross predictions

  # BATCH EFFECT
  # calculate batch effect (BE) measure
  be <- df.cross.predictions %>% 
    filter(Group=='CTR') %>% 
    filter(dataset.id!=dataset.id.train) %>% 
    mutate(cp.auroc=case_when(cp.auroc < 0.5~0.5, TRUE~cp.auroc)) %>% 
    mutate(batch_effect=abs(0.5-cp.auroc)*2)

  g <- be %>% 
    mutate(lab=sprintf(fmt='%.2f', batch_effect)) %>% 
    mutate(lab=case_when(batch_effect < 0.75~lab, TRUE~'')) %>% 
    mutate(case.train=factor(case.train, levels = names(disease_colours))) %>%
    mutate(batch_effect=case_when(batch_effect<0~0, TRUE~batch_effect)) %>% 
    mutate(case.test=str_extract(dataset.test, '_[A-Z0-9a-z]*$')) %>%
    mutate(case.test=str_remove(case.test, '_')) %>% 
    mutate(case.test=factor(case.test, levels = names(disease_colours))) %>% 
    arrange(case.test) %>% 
    mutate(dataset.id=factor(dataset.id, levels = unique(dataset.id))) %>% 
    ggplot(aes(x=dataset.train, y=dataset.id, fill=batch_effect)) + 
    geom_tile() + 
    scale_fill_gradientn(colours =rev(c('white', 
                                        embl.palette.data$sequential$Blue$value)), 
                         limits=c(0, 1), name='Batch effect')  +
    geom_text(aes(label=lab), colour='white', size=1) +
    theme_publication() +
    theme(panel.grid = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust=1),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = 'lightgrey', colour = NA),
          panel.border = element_blank()) +
    facet_grid(~case.train, scales = 'free', space = 'free') +
    xlab('Training set') + ylab('Test set') + 
    NULL
  ggsave(g, filename = here('figures', 'cross_prediction', 
                            paste0('batch_effect_heatmap_', 
                                   ml.method, '.pdf')), 
         width = 170, height = 90, units = 'mm')
  
  # DISEASE CROSS-PREDICTION PROPENSITY
  # calculate disease cross-prediction
  dcp <- df.cross.predictions %>%
    filter(Group!='CTR') %>%
    mutate(dcp=frac)
  g <- dcp %>%
    mutate(lab=sprintf(fmt='%.2f', dcp)) %>%
    mutate(lab=case_when(dcp < 0.10~'', TRUE~lab)) %>%
    mutate(Group=factor(Group, levels = names(disease_colours))) %>%
    mutate(case.train=factor(case.train,
                             levels = names(disease_colours))) %>%
    ggplot(aes(x=dataset.train, y=dataset.test, fill=dcp)) +
    geom_tile() +
    facet_grid(Group~case.train, scales='free', space = 'free') +
    geom_text(aes(label=lab), colour='white', size=1) +
    theme_publication() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust=1),
          axis.ticks = element_blank(),
          panel.border = element_blank(),
          panel.background = element_rect(fill = 'lightgrey', colour = NA)) +
    xlab('Training set') + ylab('Test set') +
    NULL
  ggsave(width = 210, height = 140, units = 'mm', 
         g + scale_fill_gradientn(
           colours = c('lightgrey', embl.palette.data$sequential$Red$value),
           name='Cross prediction'),
         filename = here('figures', 'cross_prediction',
                         paste0('cross_pred_heatmap_r_', ml.method, '.pdf')))
  ggsave(width = 210, height = 140, units = 'mm', 
         g + scale_fill_gradientn(
           colours = c('lightgrey', embl.palette.data$sequential$Green$value),
           name='Cross prediction'),
         filename = here('figures', 'cross_prediction',
                         paste0('cross_pred_heatmap_g_', ml.method, '.pdf')))
  
  
  # FPR as boxplots
  g <- df.cross.predictions %>% 
    select(dataset.id, Group, dataset.train, frac) %>% 
    distinct() %>% 
    mutate(Group=factor(Group, levels = c('CTR', names(disease_colours)))) %>% 
    ggplot(aes(x=Group, y=frac, fill=Group)) + 
      geom_boxplot(outlier.shape = NA) + 
      geom_jitter(width = 0.08) + 
      scale_fill_manual(values = c('CTR'='grey', disease_colours), 
                        name='Group', guide=FALSE) + 
      xlab('') + 
      ylab('False positive rate') + 
      theme_publication()
  
  # External AUROC for data transfer for CRC/UC/CD
  df <- df.cross.predictions %>% 
    select(dataset.id, dataset.id.train, Group, case.train, 
           dataset.train, ext.auroc, auroc) %>% 
    filter(Group %in% c('CRC', 'UC', "CD")) %>% 
    filter(case.train %in% c('CRC', 'UC', "CD")) %>% 
    filter(Group==case.train)
    
  original.auc <- df %>% 
    select(dataset.train, case.train, auroc) %>% 
    distinct() %>% 
    arrange(auroc)
  
  
  g <- df %>% 
    group_by(dataset.train, case.train) %>% 
    summarise(m=mean(ext.auroc, na.rm=TRUE), .groups = 'drop') %>% 
    mutate(dataset.train=
             factor(dataset.train, 
                    levels = original.auc$dataset.train)) %>% 
    ggplot(aes(x=dataset.train, y=m, fill=case.train))  +
      geom_bar(stat='identity') + 
      facet_grid(~case.train, scales = 'free', space = 'free') + 
      coord_cartesian(ylim=c(0.5, 1)) + 
      geom_point(data=original.auc, aes(x=dataset.train, y=auroc), pch=8,
                 colour='#E40046') +
      geom_jitter(data=df, aes(x=dataset.train, y=ext.auroc), width = 0.08,
                  col='#54585A') + 
      theme_publication() + 
      scale_fill_manual(values=disease_colours)
}
