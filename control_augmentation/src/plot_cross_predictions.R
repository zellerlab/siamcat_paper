# ##############################################################################
#
## Plot results from the control augmentation
#
# ##############################################################################

library("here")
library("tidyverse")
library("ggembl")
library("pROC")
library("yaml")
library("SIAMCAT")

colour.list <- yaml.load_file(here('parameter_space', 'data_info', 
                                   'colours.yaml'))
disease_colours <- unlist(colour.list$disease)

# loop through the machine learning methods, i guess
ml.method <- 'enet-0.5'

# get all the data
# naive
df.naive <- read_tsv(here('parameter_space', 'files', 
                          paste0('cross_prediction_', ml.method, '.tsv'))) %>% 
  mutate(type='naive')
# augmented
df.augm <- map(c('cohort_2',  'cohort_5', 'other_ctr', 'random', 'similar'),
               .f=function(x){read_tsv(here('control_augmentation', 'files',
                                            paste0('cross_prediction_', 
                                                   ml.method, '_', x, 
                                                   '.tsv'))) %>% 
                   mutate(ctr_type=x)}) %>% 
  bind_rows() %>% 
  mutate(type='augmented')

df.all <- bind_rows(df.naive, df.augm)

# ##############################################################################
# plot everything

ctr_type_levels <- c('naive', 'similar',  'cohort_2', 'cohort_5', 
                     'other_ctr', 'random')
ctr_type_colours <- c('naive'='#009F4D', 'similar'='#FFA300', 
                      'random'='#E40046', 'other_ctr'='#8246AF', 
                      'cohort_2'='#307FE2', 'cohort_5'='#307FE2')
# auc values
df.all %>% 
  group_by(dataset.train, type, ctr_type) %>% 
  summarize(auc=unique(auroc), .groups='drop') %>% 
  mutate(ctr_type=case_when(is.na(ctr_type)~'naive', TRUE~ctr_type)) %>% 
  select(-type) %>% 
  pivot_wider(names_from = ctr_type, values_from=auc) %>% 
  pivot_longer(cols = -c(dataset.train, naive), names_to='ctr_type', 
               values_to='auc') %>% 
  mutate(ctr_type=factor(ctr_type, levels = ctr_type_levels)) %>% 
  ggplot(aes(x=naive, y=auc)) + 
    geom_abline(slope = 1, intercept = 0, col='darkgrey') + 
    geom_point() + 
    facet_grid(~ctr_type) + 
    xlab('AUROC for naive models') + 
    ylab("AUROC for control-augmented models") + 
    theme_publication()
    
df.all %>% 
  group_by(dataset.train, type, ctr_type) %>% 
  summarize(auc=unique(auroc), .groups='drop') %>% 
  mutate(ctr_type=case_when(is.na(ctr_type)~'naive', TRUE~ctr_type)) %>% 
  select(-type) %>% 
  pivot_wider(names_from = ctr_type, values_from=auc) %>% 
  pivot_longer(cols = -c(dataset.train, naive), names_to='ctr_type', 
               values_to='auc') %>% 
  mutate(diff=auc-naive) %>% 
  mutate(ctr_type=factor(ctr_type, levels = ctr_type_levels)) %>% 
  ggplot(aes(x=ctr_type, y=diff)) + 
    geom_hline(yintercept = 0, col='darkgrey') + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(width = 0.08) + 
    xlab('') + 
    ylab("AUROC improvement") + 
    theme_publication()

# FPR values
g <- df.all %>%
  filter(ctr_type %in% c(NA_character_, 'cohort_5')) %>% 
  select(dataset.id, dataset.test, Group, dataset.train, 
         case.train, frac, type) %>% 
  filter(case.train!=Group) %>% 
  distinct() %>% 
  mutate(is_ctr=Group=='CTR') %>% 
  mutate(Group=str_extract(dataset.test, '_[a-zA-Z0-9]*$')) %>% 
  mutate(Group=str_remove(Group, '_')) %>% 
  mutate(Group=factor(Group, levels = c(names(disease_colours)))) %>% 
  mutate(type=factor(type, levels = c('naive', 'augmented'))) %>% 
  ggplot(aes(x=Group, y=frac, fill=type)) + 
    geom_boxplot(outlier.shape = NA, colour='black') + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.08),
                pch=16, colour='#00000050') + 
    scale_fill_manual(values = c('#307FE2', '#E40046'), name='Type') + 
    xlab('') + 
    ylab('False positive rate') + 
    theme_publication() + 
    geom_hline(yintercept = 0.1) + 
    facet_grid(is_ctr~.)
ggsave(g, filename = here('figures', 'cross_prediction', 'fpr_plot.pdf'),
       useDingbats=FALSE, width = 190, height = 100, units = 'mm')

# external AUC VALUES
# External AUROC for data transfer for CRC/UC/CD
df <- df.all %>% 
  filter(ctr_type %in% c(NA_character_, 'cohort_5')) %>% 
  select(dataset.test, dataset.id, dataset.id.train, Group, case.train, 
         dataset.train, ext.auroc, auroc, type) %>% 
  filter(Group %in% c('CRC', 'UC', "CD")) %>% 
  filter(case.train %in% c('CRC', 'UC', "CD")) %>% 
  filter(Group==case.train) %>% 
  mutate(type=factor(type, levels = c('naive', 'augmented'))) %>% 
  mutate(case.train=factor(case.train, levels = c('CRC', 'CD', 'UC')))

order <- df %>% 
  group_by(dataset.test, case.train) %>% 
  summarise(m=mean(ext.auroc, na.rm=TRUE), .groups = 'drop') %>% 
  arrange(m)


g <- df %>% 
  group_by(dataset.test, case.train, type) %>% 
  summarise(m=mean(ext.auroc, na.rm=TRUE), 
            s=sd(ext.auroc, na.rm=TRUE),
            .groups = 'drop') %>% 
  mutate(up=m+s, low=m-s) %>% 
  mutate(dataset.test=factor(dataset.test, levels = order$dataset.test)) %>% 
  ggplot(aes(x=dataset.test, y=m, fill=type))  +
    geom_bar(stat='identity', position = position_dodge()) + 
    geom_errorbar(aes(ymax=up, ymin=low), position = position_dodge()) +
    geom_jitter(data=df, aes(x=dataset.test, y=ext.auroc), 
                position = position_jitterdodge(jitter.width = 0.08),
                pch=16, colour='#00000050') + 
    facet_grid(~case.train, scales = 'free', space = 'free') + 
    coord_cartesian(ylim=c(0.5, 1)) + 
    scale_fill_manual(values = c('#307FE2', '#E40046'), name='Type') + 
    theme_publication() + 
    xlab('Test set') + 
    ylab('Model transfer AUROC') +
    NULL
ggsave(g, filename = here('figures', 'cross_prediction', 'external_auc.pdf'),
       useDingbats=FALSE, width = 190, height = 70, units = 'mm')
# test?
temp <- df %>% 
  select(dataset.test, dataset.id.train,case.train, ext.auroc, type) %>% 
  pivot_wider(names_from = type, values_from = ext.auroc) %>% 
  filter(complete.cases(.)) 
wilcox.test(temp %>% filter(case.train=='CRC') %>% pull(naive), 
            temp %>% filter(case.train=='CRC') %>% pull(augmented),  
            paired = TRUE)
wilcox.test(temp %>% filter(case.train=='CD') %>% pull(naive), 
            temp %>% filter(case.train=='CD') %>% pull(augmented),  
            paired = TRUE)
wilcox.test(temp %>% filter(case.train=='UC') %>% pull(naive), 
            temp %>% filter(case.train=='UC') %>% pull(augmented),  
            paired = TRUE)
  
motu.tasks.sel <- df.all %>%
  select(dataset.id.train, case.train, auroc, ctr_type) %>%
  distinct() %>%
  filter(ctr_type=='cohort_5') %>% 
  filter(auroc > 0.7) %>%
  mutate(id=paste0(dataset.id.train, '_', case.train))


# batch effects for all ctr-augmentations
be <- df.all %>% 
  filter(Group=='CTR') %>% 
  filter(dataset.id!=dataset.id.train) %>% 
  mutate(cp.auroc=case_when(cp.auroc < 0.5~0.5, TRUE~cp.auroc)) %>% 
  mutate(batch_effect=abs(0.5-cp.auroc)*2) %>% 
  filter(dataset.train%in%motu.tasks.sel$id)
g3 <- be %>% 
  mutate(ctr_type=case_when(is.na(ctr_type)~'naive', TRUE~ctr_type)) %>% 
  mutate(ctr_type=factor(ctr_type, levels = ctr_type_levels)) %>% 
  ggplot(aes(x=ctr_type, y=batch_effect, fill=ctr_type))  +
    geom_boxplot() + 
    theme_publication() + 
    xlab('') + ylab('Cross-study portability') + 
    scale_fill_manual(values = ctr_type_colours, guide=FALSE)

# try out
g.final.a <- be %>% 
  select(dataset.id, dataset.train, ctr_type, batch_effect) %>% 
  mutate(ctr_type=case_when(is.na(ctr_type)~'naive', TRUE~ctr_type)) %>% 
  distinct() %>% 
  pivot_wider(names_from=ctr_type, values_from=batch_effect) %>% 
  # select(-naive) %>% 
  pivot_longer(cols=-c(dataset.id, dataset.train, cohort_5), 
               names_to='ctr_type',
               values_to='batch_effect') %>% 
  mutate(ctr_type=factor(ctr_type, levels = ctr_type_levels)) %>% 
  ggplot(aes(x=cohort_5, y=batch_effect)) + 
    geom_abline(slope = 1, intercept = 0) + 
    geom_point(alpha=0.5) + 
    facet_grid(~ctr_type) + 
    theme_publication() + 
    xlab(paste0('Cross-study portability after augmentation with five ', 
                'times control samples')) + 
    ylab('Cross-study portability') + 
    xlim(0,1 ) + ylim(0,1)
# quantify delta?
be %>% 
  select(dataset.id, dataset.train, ctr_type, batch_effect) %>% 
  mutate(ctr_type=case_when(is.na(ctr_type)~'naive', TRUE~ctr_type)) %>% 
  distinct() %>% 
  pivot_wider(names_from=ctr_type, values_from=batch_effect) %>% 
  pivot_longer(cols=-c(dataset.id, dataset.train, cohort_5), 
               names_to='ctr_type',
               values_to='batch_effect') %>% 
  mutate(diff=cohort_5-batch_effect) %>% 
  group_by(ctr_type) %>% 
  summarise(m=mean(diff, na.rm = TRUE))


# cross-predictions for all ctr-augmentations
dcp <- df.all %>%
  filter(Group!='CTR') %>%
  mutate(dcp=frac) %>% 
  filter(dataset.train%in%motu.tasks.sel$id)
g4 <- dcp %>% 
  mutate(ctr_type=case_when(is.na(ctr_type)~'naive', TRUE~ctr_type)) %>% 
  mutate(ctr_type=factor(ctr_type, levels = ctr_type_levels)) %>% 
  mutate(same=case.train!=Group) %>% 
  ggplot(aes(x=ctr_type, y=dcp, fill=ctr_type))  +
  geom_boxplot() + 
  theme_publication() + 
  xlab('') + ylab('Prediction rate\nfor other diseases') + 
  scale_fill_manual(values = ctr_type_colours, guide=FALSE) + 
  facet_grid(~same) + 
  NULL
g5 <- dcp %>% 
  mutate(ctr_type=case_when(is.na(ctr_type)~'naive', TRUE~ctr_type)) %>% 
  mutate(ctr_type=factor(ctr_type, levels = ctr_type_levels)) %>% 
  filter(case.train!=Group) %>% 
  ggplot(aes(x=case.train, y=dcp, fill=case.train)) + 
    geom_boxplot() + 
    facet_grid(~ctr_type) + 
    theme_publication() + 
    xlab('') + 
    ylab('Prediction rate\nfor other diseases') + 
    scale_fill_manual(values=disease_colours, name='Disease') + 
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

g.final.b <- dcp %>% 
  select(dataset.test, dataset.train, ctr_type, dcp, case.train, Group) %>% 
  mutate(ctr_type=case_when(is.na(ctr_type)~'naive', TRUE~ctr_type)) %>% 
  distinct() %>% 
  pivot_wider(names_from=ctr_type, values_from=dcp) %>% 
  pivot_longer(cols=-c(dataset.test, dataset.train, cohort_5, Group, case.train), 
               names_to='ctr_type', values_to='dcp') %>%
  mutate(ctr_type=factor(ctr_type, levels = ctr_type_levels)) %>% 
  mutate(same=case.train==Group) %>% 
  # filter(same) %>%
  ggplot(aes(x=cohort_5, y=dcp)) + 
    geom_abline(slope = 1, intercept = 0) + 
    geom_point(alpha=0.5) + 
    facet_grid(same~ctr_type) + 
    theme_publication() + 
    xlab(paste0('Prediction rate of models augmented with five ', 
                'times control samples')) + 
    ylab('Prediction rate') + 
    xlim(0,1 ) + ylim(0,1)
dcp %>% 
  select(dataset.test, dataset.train, ctr_type, dcp, case.train, Group) %>% 
  mutate(ctr_type=case_when(is.na(ctr_type)~'naive', TRUE~ctr_type)) %>% 
  distinct() %>% 
  pivot_wider(names_from=ctr_type, values_from=dcp) %>% 
  pivot_longer(cols=-c(dataset.test, dataset.train, cohort_5, 
                       Group, case.train), 
               names_to='ctr_type', values_to='dcp') %>% 
  mutate(diff=cohort_5-dcp) %>% 
  mutate(same=case.train==Group) %>% 
  group_by(ctr_type, same) %>% 
  summarise(m=mean(diff, na.rm=TRUE))


# save final figures
ggsave(g.final.a, filename = here('figures', 'cross_prediction', 
                                  paste0('batch_effect_comp_', ml.method, 
                                         '.pdf')),
       width = 170, height = 55, units = 'mm', useDingbats=FALSE)
ggsave(g.final.b, filename = here('figures', 'cross_prediction', 
                                paste0('cross_pred_comp_', ml.method, 
                                       '.pdf')),
       width = 170, height = 80, units = 'mm', useDingbats=FALSE)


# save figures
ggsave(g, filename = here('figures', 'cross_prediction', 
                          paste0("auroc_ctr_augm_", 
                                 ml.method, ".pdf")),
       width = 5, height = 4, useDingbats=FALSE)
ggsave(g2, filename = here('figures', 'cross_prediction', 
                          paste0("auroc_change_ctr_augm_", 
                                 ml.method, ".pdf")),
       width = 4, height = 3, useDingbats=FALSE)
ggsave(g3, filename = paste0('~/Desktop/',
                           paste0("batch_effect_improvement_", 
                                  ml.method, ".png")),
       width = 4, height = 3)
ggsave(g4, filename = paste0('~/Desktop/',
                           paste0("cross_pred_improvement_2_", 
                                  ml.method, ".png")),
       width = 4, height = 3)
ggsave(g5, filename = here('figures', 'cross_prediction', 
                           paste0("cross_pred_improvement_box_", 
                                  ml.method, ".pdf")),
       width = 7, height = 3, useDingbats=FALSE)

# loop through the ctr augmentation types
for (x in ctr_type_levels){
  if (x=='naive') next()
  
  # batch effect heatmap for each augmentation
  g <- be %>% 
    filter(ctr_type==x) %>% 
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
                            paste0('batch_effect_heatmap_', ml.method, 
                                   '_', x, '.pdf')),
         width = 120, height = 80, units = 'mm')
  # cross-prediction for each augmentation
  g <- dcp %>%
    filter(ctr_type==x) %>% 
    mutate(lab=sprintf(fmt='%.2f', dcp)) %>%
    mutate(lab=case_when(dcp < 0.1~'', TRUE~lab)) %>%
    mutate(case.train=factor(case.train, levels = names(disease_colours))) %>% 
    mutate(Group=factor(Group, levels = names(disease_colours))) %>% 
    mutate(dataset.test=factor(dataset.test)) %>% 
    mutate(dataset.train=factor(dataset.train, 
                                levels = rev(unique(dataset.train)))) %>% 
    ggplot(aes(x=dataset.train, y=dataset.test, fill=dcp)) +
    geom_tile() +
    geom_text(aes(label=lab), colour='white', size=1) +
    theme_publication() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust=1),
          axis.ticks = element_blank(),
          panel.border = element_blank(),
          panel.background = element_rect(fill = 'lightgrey', colour = NA)) +
    xlab('Training set') + ylab('Test set') +
    facet_grid(Group~case.train, scales = 'free', space = 'free') +
    NULL
  ggsave(width = 170, height = 140, units = 'mm', 
         g + scale_fill_gradientn(
           colours = c('lightgrey', embl.palette.data$sequential$Red$value),
           name='Cross prediction'),
         filename = here('figures', 'cross_prediction',
                         paste0('cross_pred_heatmap_ctr_r_', ml.method, 
                                '_', x, '.pdf')))
  ggsave(width = 170, height = 140, units = 'mm', 
         g + scale_fill_gradientn(
           colours = c('lightgrey', embl.palette.data$sequential$Green$value),
           name='Cross prediction'),
         filename = here('figures', 'cross_prediction',
                         paste0('cross_pred_heatmap_ctr_g_', ml.method, 
                                '_', x, '.pdf')))
}

# ##############################################################################
# plot some actual examples

motu.tasks <- df.all %>%
  select(dataset.id.train, case.train) %>%
  distinct() %>%
  mutate(id=paste0(dataset.id.train, '_', case.train))


pdf(here('figures', 'cross_prediction', 'examples.pdf'), 
    width = 5, height = 3, useDingbats = FALSE)
for (ref.study in motu.tasks$id){
  message("#----------------------------------------------------\n", 
          ref.study, '\n#----------------------------------------------------')
  fn.model.naive <- list.files(here('parameter_space','models'), 
                               pattern=ref.study, full.names = TRUE) %>% 
    grep(pattern=ml.method, value=TRUE)
  load(fn.model.naive)

  df.plot <- enframe(rowMeans(pred_matrix(sc.obj.train)), 
                     name='Sample_ID', value='prediction') %>% 
    left_join(enframe(label(sc.obj.train)$label, name='Sample_ID', 
                      value='label'),
              by='Sample_ID') %>% 
    mutate(label=case_when(label==-1~names(label(sc.obj.train)$info)[1],
                           label==1~names(label(sc.obj.train)$info)[2])) %>% 
    mutate(study=ref.study)
  # 10% threshold
  t <- eval_data(sc.obj.train)$roc$threshold[which(
    eval_data(sc.obj.train)$roc$specificities > 0.9)[1]]
  
  
  for (ext.study in motu.tasks %>% filter(id!=ref.study) %>% 
       pull(dataset.id.train) %>% unique()){
    message(ext.study)
    fn.ext <- list.files(here('parameter_space', 'sc'), 
                         pattern=ext.study, full.names=TRUE) %>% 
      grep(pattern='mOTUs2', value=TRUE) %>% 
      sample(1)
    load(fn.ext)
    sc.obj <- make.predictions(sc.obj.train, sc.obj, verbose=0)
    df.plot.ext <- enframe(rowMeans(pred_matrix(sc.obj)), 
                           name='Sample_ID', value='prediction') %>% 
      left_join(enframe(label(sc.obj)$label, name='Sample_ID', 
                        value='label'),
                by='Sample_ID') %>% 
      mutate(label=case_when(label==-1~names(label(sc.obj)$info)[1],
                             label==1~names(label(sc.obj)$info)[2])) %>% 
      mutate(study=ext.study)
    # calculate AUROC
    x <- auc(cases=df.plot %>% filter(label!='CTR') %>% pull(prediction),
             controls=df.plot.ext %>% filter(label=='CTR') %>% pull(prediction),
             direction='<')
                   
    # calculate prediction rate for other disease
    dcp <- df.plot.ext %>% 
      mutate(x=prediction>=t) %>% 
      group_by(label) %>% 
      summarise(m=mean(x), .groups='drop') %>% filter(label!='CTR') %>% 
      pull(m) %>% sprintf(fmt='%.2f')
    g <- bind_rows(df.plot, df.plot.ext) %>% 
      mutate(study=factor(study, levels = c(ref.study, ext.study))) %>% 
      mutate(label=factor(label, levels = c('CTR', names(disease_colours)))) %>% 
      ggplot(aes(x=study, y=prediction, fill=label)) + 
      geom_boxplot() + 
      scale_fill_manual(values=c("CTR"='darkgrey', disease_colours), 
                        guide=FALSE) + 
      theme_publication(panel.grid = 'major_y') + 
      ylab("Model prediction") + 
      xlab('') + ggtitle(paste0('cross prediction: ', dcp, 
                                ', AUC: ', sprintf(fmt='%.2f', x)))
    print(g)
  }
}
dev.off()
