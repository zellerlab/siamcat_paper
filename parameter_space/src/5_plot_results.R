# ##############################################################################
#
## Plot results
#
# ##############################################################################

library("tidyverse")
library("ggembl")
library("cowplot")
library("here")
library("SIAMCAT")
library("yaml")


colour.list <- yaml.load_file(here('parameter_space', 'data_info', 
                                   'colours.yaml'))
disease_colours <- unlist(colour.list$disease)
type_colours <- unlist(colour.list$type)

# all.results
df.everything <- read_tsv(here("parameter_space", "job_info", 
                               "full_results.tsv"),
                          col_types = cols(
                            .default = col_character(),
                            log.n0 = col_double(),
                            sd.min.q = col_double(),
                            fs = col_logical(),
                            fs.cutoff = col_double(),
                            job.id = col_double(),
                            processed = col_logical(),
                            auroc = col_double(),
                            time = col_double(),
                            ci.low = col_double(),
                            ci.high = col_double()
                          ))
data.info <- read_tsv(here("parameter_space", "data_info", "all_tasks.tsv"))

# ##############################################################################
# plot best-performing parameter set and export data

auc.list <- list()
max.size <- data.info %>% 
  group_by(type) %>% 
  summarise(n=n()) %>% 
  pull(n) %>% 
  max %>% +.5
for (data.type in unique(df.everything$type)){
  
  best.jobs <- df.everything %>% 
    filter(type==data.type) %>% 
    group_by(job.id) %>% 
    summarise(n=sum(!is.na(auroc)), auc=mean(auroc)) %>% 
    filter(!is.na(auc)) %>%
    arrange(desc(auc))
    
  info <- df.everything %>% 
    filter(type==data.type) %>% 
    filter(job.id==best.jobs$job.id[1]) %>% 
    select(1:10) %>% 
    distinct()

  best.points <- df.everything %>% 
    filter(type==data.type) %>%
    filter(job.id==best.jobs$job.id[1])

  # ############################################################################
  # extract feature matrix and metadata
  best.points$n.feat <- NA_real_
  for (x in seq_len(nrow(best.points))){
    temp <- best.points[x,]
    fn <- paste0(temp$dataset.id, '_', temp$case, '_', temp$type)
    message(fn)
    load(here('parameter_space', 'sc', paste0('sc_', fn, '.RData')))
    
    # perform feature filtering
    if (str_detect(info$filt.method, ';')){
      methods <- str_split(info$filt.method, ';')[[1]]
      cutoffs <- str_split(info$filt.cutoff, ';')[[1]]
      sc.obj <- filter.features(sc.obj,
                                filter.method = methods[1],
                                cutoff = as.numeric(cutoffs[1]), 
                                verbose=0)
      sc.obj <- filter.features(sc.obj,
                                filter.method=methods[2],
                                cutoff = as.numeric(cutoffs[2]),
                                feature.type = 'filtered', verbose=0)
    } else {
      sc.obj <- filter.features(sc.obj, 
                                filter.method = info$filt.method,
                                cutoff = as.numeric(info$filt.cutoff),
                                verbose=0)
    }
    feat.mat <- get.filt_feat.matrix(sc.obj)
    best.points$n.feat[x] <- nrow(feat.mat)
    write.table(feat.mat, file = here('parameter_space', 'files', 
                                      'feat', paste0('feat_rel_', fn, '.tsv')),
                sep='\t', quote = FALSE)
  }
  
  best.points <- best.points %>% 
    left_join(data.info)
  
  # plot results
  g1 <- best.points %>% 
    mutate(n.all=n.case+n.ctr) %>% 
    mutate(case.split=factor(case.split, 
                             names(sort(-table(case.split))))) %>% 
    arrange(desc(case.split), case, n.all) %>% 
    mutate(label.nice=factor(label.nice, levels = label.nice)) %>% 
    ggplot(aes(y=label.nice, x=auroc, col=case.nice)) + 
      theme_publication(panel.grid = 'major_x') + 
      # theme(panel.grid.major.y = element_blank()) +
      geom_segment(aes(y=label.nice, yend=label.nice,
                       x = ci.low, xend=ci.high), size=0.75) + 
      geom_point(shape=23, size=1.5, col='black', aes(fill=case)) +
      scale_colour_manual(values=disease_colours, guide=FALSE) +
      scale_fill_manual(values=disease_colours, guide=FALSE) +
      coord_cartesian(xlim=c(0.5, 1), ylim=c(0.5, max.size), expand = FALSE) + 
      xlab('AUROC') + 
      ylab('')
  ymax <- ifelse(data.type=='RDP', 150, 200)
  g2 <- best.points %>% 
    mutate(n.all=n.case+n.ctr) %>% 
    mutate(case.split=factor(case.split, 
                             names(sort(-table(case.split))))) %>% 
    arrange(desc(case.split), case, n.all) %>% 
    mutate(label.nice=factor(label.nice, levels = label.nice)) %>% 
    ggplot(aes(x=label.nice, y=n.all, col=case.nice)) +
      geom_bar(stat='identity', fill='white') + 
      theme_publication(panel.grid = 'major_x') + 
      theme(panel.grid.major.y = element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) +
      coord_flip(ylim=c(0, ymax)) + 
      scale_color_manual(values=disease_colours, guide=FALSE) + 
      xlab('') + ylab('N samples') + 
      geom_text(aes(label=n.all, y=155), col='black')
  
  g.all <- plot_grid(g1, g2, rel_widths = c(0.8, 0.2))
  ggsave(g.all, filename = here('figures', 'parameter_space',
                                paste0('best_jobs_', data.type,'.pdf')),
         useDingbats=FALSE, width = 110, height = 75, units = 'mm')
  
  # save AUCs
  auc.list[[data.type]] <- best.points %>% 
    select(dataset.id, type, case, auroc, time, ci.low, ci.high, 
           case.split, dataset.nice, case.nice, 
           label.nice, n.feat, n.ctr, n.case)
}

auroc.all <- bind_rows(auc.list)
write_tsv(auroc.all, path = here('parameter_space', 'files', 'auroc_all.tsv'))

# ##############################################################################
# plot func vs tax
df.temp <- auroc.all %>% 
  mutate(tax=case_when(type%in%c('mOTUs2', 'metaphlan')~auroc, 
                       TRUE~NA_real_)) %>% 
  mutate(func=case_when(type%in%c('eggNOG', 'humann2')~auroc, 
                        TRUE~NA_real_)) %>%
  filter(!is.na(tax) | !is.na(func)) %>% 
  mutate(type=case_when(type%in%c('metaphlan', 'humann2')~'meta', 
                        TRUE~'motus')) %>% 
  group_by(dataset.nice, type, case.nice, n.ctr, n.case) %>% 
  summarise(tax=na.omit(unique(tax))[1], func=na.omit(unique(func))[1])

cor(df.temp$tax, df.temp$func)
cor.test(df.temp$tax, df.temp$func)

g <- df.temp %>% 
  mutate(n.all=n.ctr+n.case) %>% 
  ggplot(aes(x=tax, y=func, shape=type, col=case.nice, size=n.all)) + 
    geom_point(alpha=0.8) + 
    geom_abline(intercept = 0, slope = 1, col='darkgrey') +
    theme_publication(panel.grid = 'major') + 
    coord_cartesian(xlim=c(0.5, 1), ylim=c(0.5, 1), expand = FALSE) + 
    xlab('AUROC based on taxonomic profiles') + 
    ylab('AUROC based on functional profiles') + 
    scale_shape_manual(values=c(17, 19), guide=FALSE) + 
    scale_colour_manual(values = disease_colours, guide=FALSE) +
    NULL
ggsave(g, filename = here('figures', 'parameter_space', 'cor_tax_func.pdf'),
       width = 100, height = 85, units = 'mm', useDingbats=FALSE)

# ##############################################################################
# plot dataset size vs ci
g <- auroc.all %>% 
  mutate(n.all=n.ctr+n.case) %>% 
  mutate(ci=ci.high-ci.low) %>% 
  ggplot(aes(x=n.all, y=ci, col=type)) + 
    geom_point(alpha=0.8) + 
    theme_publication(panel.grid = 'major') + 
    xlab('Dataset size') + ylab('Range of 95% confidence interval') +
    scale_colour_manual(values = type_colours, guide=FALSE) +
    scale_shape_manual(values=c(8, 2, 1, 17, 19)) + 
    scale_y_continuous(limits = c(0, 0.4), expand = c(0, 0)) + 
    scale_x_log10() +
    NULL
ggsave(g, filename = here('figures', 'parameter_space', 'cor_size_ci.pdf'),
       width = 80, height = 80, units = 'mm', useDingbats=FALSE)

# plot dataset size vs AUROC
g <- auroc.all %>% 
  mutate(n.all=n.ctr+n.case) %>% 
  ggplot(aes(x=n.all, y=auroc, col=type)) + 
  geom_point(alpha=0.8) + 
  theme_publication(panel.grid = 'major') + 
  xlab('Dataset size') + ylab('Best AUROC') +
  scale_colour_manual(values = type_colours, guide=FALSE) +
  scale_shape_manual(values=c(8, 2, 1, 17, 19)) + 
  scale_y_continuous(limits = c(0.5, 1), expand = c(0, 0)) +
  scale_x_log10() +
  NULL
ggsave(g, filename = here('figures', 'parameter_space', 'cor_size_auc.pdf'),
       width = 80, height = 80, units = 'mm', useDingbats=FALSE)

auroc.all %>% 
  mutate(n.all=n.ctr+n.case) %>% 
  mutate(size=n.all >= 100) %>% 
  group_by(size) %>% 
  summarise(n.over=sum(auroc >= 0.75), n.total=n()) %>% 
  mutate(prob=n.over/n.total)

# ##############################################################################
# plot mOTUs and metaphlan
df.temp <- auroc.all %>% 
  filter(type%in%c('mOTUs2', 'metaphlan')) %>% 
  select(type, dataset.nice, case.nice, auroc, label.nice) %>% 
  spread(type, auroc) %>% 
  drop_na() %>% 
  mutate(m=(metaphlan+mOTUs2)/2) %>% 
  arrange(case.nice, m)

wilcox.test(df.temp$metaphlan, df.temp$mOTUs2, paired = TRUE)

g <- df.everything %>% 
  filter(type%in%c('mOTUs2', 'metaphlan')) %>% 
  filter(dataset.nice%in%df.temp$dataset.nice) %>% 
  filter(case.nice%in%df.temp$case.nice) %>% 
  mutate(label.nice=factor(label.nice, levels = df.temp$label.nice)) %>% 
  ggplot(aes(x=label.nice, y=auroc, fill=type)) + 
    geom_boxplot(outlier.alpha = 0.1) + 
    theme_publication(panel.grid = 'major_y') + 
    xlab('') + ylab('AUROC') + 
    scale_y_continuous(limits = c(0.5, 1), expand = c(0,0)) +
    scale_fill_manual(values = type_colours, guide=FALSE) + 
    geom_point(data=df.temp %>% 
                 mutate(label.nice=factor(label.nice, levels = label.nice)), 
               aes(x=label.nice, y=metaphlan), inherit.aes = FALSE, shape=17) +
    geom_point(data=df.temp %>% 
                 mutate(label.nice=factor(label.nice, levels = label.nice)), 
               aes(x=label.nice, y=mOTUs2), inherit.aes = FALSE, shape=19) +
    theme(axis.text.x = element_text(angle=45, hjust=1)) + 
    NULL
ggsave(g, filename = here('figures', 'parameter_space', 'motus_metaphlan.pdf'),
       width = 100, height = 85, units = 'mm', useDingbats=FALSE)

# ##############################################################################
# parameter space figures

# machine learning method
x <- df.everything %>% 
  filter(!is.na(auroc)) %>% 
  filter(!is.infinite(auroc)) %>% 
  group_by(label.nice, type, ml.method) %>% 
  summarise(m=max(auroc)) %>% 
  ungroup() %>% 
  mutate(ml.method=
           factor(ml.method,
                  levels = rev(c('randomForest', 'lasso', 'enet')))) %>%
  mutate(type=factor(type, 
                     levels = c('metaphlan', 'humann2', 'mOTUs2', 
                                'eggNOG', 'RDP'))) %>% 
  filter(!is.na(ml.method))

df.test <- x %>%
  pivot_wider(names_from = ml.method, values_from = m)
# lasso randomForest
wilcox.test(df.test$lasso, df.test$randomForest, paired = TRUE)
# lasso enet
wilcox.test(df.test$lasso, df.test$enet, paired = TRUE)
# enet randomForest
wilcox.test(df.test$enet, df.test$randomForest, paired = TRUE)


g.ml.method <- x %>% 
  ggplot(aes(x=ml.method, y=m, fill=ml.method)) +
    geom_hline(yintercept = 0.5, colour='darkgrey') + 
    # geom_line(aes(group=id),col='lightgrey', size=.7) +
    geom_boxplot() + 
    theme_publication() + 
    scale_fill_embl(guide=FALSE) +
    # facet_grid(type~.) + 
    xlab('') + ylab('Best AUROC') + ylim(0.2, 1) + 
    coord_flip(ylim = c(0,1), expand = FALSE) + 
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

## similar stuff, but other parameters
# norm.method
g.norm.method <- df.everything %>% 
  filter(!is.na(auroc)) %>% 
  filter(!is.infinite(auroc)) %>% 
  mutate(id=paste0(dataset.id, case.nice, ml.method)) %>% 
  mutate(norm.method=
           factor(norm.method,levels = c('none', 'rank.std',
                                         'log.std', 'log.clr'))) %>%
  mutate(type=factor(type, 
                     levels = c('metaphlan', 'humann2', 'mOTUs2', 
                                'eggNOG', 'RDP'))) %>% 
  mutate(ml.method=
           factor(ml.method,
                  levels = c('randomForest', 'lasso', 'enet'))) %>%
  ggplot(aes(x=ml.method, y=auroc, fill=norm.method)) +
    geom_hline(yintercept = 0.5, colour='darkgrey') +
    geom_boxplot() + 
    theme_publication() + 
    xlab('') + ylab('Best AUROC') + ylim(0.2, 1)

# filtering.method
g.filt.method <- df.everything %>% 
  filter(!is.na(auroc)) %>% 
  filter(!is.infinite(auroc)) %>% 
  group_by(dataset.id, type, case.nice, filt.method) %>% 
  summarise(m=max(auroc)) %>% 
  ungroup() %>% 
  mutate(id=paste0(dataset.id, case.nice)) %>% 
  mutate(filt.method=factor(filt.method,
                          levels = c('pass', 'abundance', 
                                     'prevalence', 
                                     'abundance;prevalence'))) %>%
  mutate(type=factor(type, 
                     levels = c('metaphlan', 'humann2', 'mOTUs2', 
                                'eggNOG', 'RDP'))) %>% 
  ggplot(aes(x=filt.method, y=m, fill=filt.method)) +
  geom_hline(yintercept = 0.5, colour='darkgrey') + 
  geom_line(aes(group=id),col='lightgrey', size=.7) +
  geom_boxplot() + 
  theme_publication() + 
  facet_grid(type~.) + 
  xlab('') + ylab('Best AUROC') + ylim(0.2, 1) + 
  coord_flip(ylim = c(0,1), expand = FALSE)

# feature selection
g.fs.cutoff <- df.everything %>% 
  filter(!is.na(auroc)) %>% 
  filter(!is.infinite(auroc)) %>% 
  mutate(fs.cutoff=case_when(is.na(fs.cutoff)~'None',
                             TRUE~as.character(fs.cutoff))) %>% 
  # group_by(dataset.id, type, case, fs.cutoff) %>%
  # summarise(m=max(auroc)) %>%
  # ungroup() %>%
  mutate(id=paste0(dataset.id, case)) %>% 
  mutate(fs.cutoff=factor(fs.cutoff,
                            levels = c('10', '25', '50', "100", '200', '400', 
                                       '500', '1000', '2000','None'))) %>%
  mutate(type=factor(type, 
                     levels = c('metaphlan', 'humann2', 'mOTUs2', 
                                'eggNOG', 'RDP'))) %>% 
  ggplot(aes(x=type, y=auroc, fill=fs.cutoff)) +
  geom_hline(yintercept = 0.5, colour='darkgrey') + 
  # geom_line(aes(group=id),col='lightgrey', size=.7) +
  geom_boxplot(outlier.alpha = 0.1, outlier.stroke = 0) + 
  theme_publication() + 
  # facet_grid(type~., space = 'free', scale='free') + 
  xlab('') + ylab('Best AUROC') + ylim(0.2, 1) + 
  # coord_flip(ylim = c(0,1), expand = FALSE)
  NULL
ggsave(g.ml.method, filename = 
         here("figures", "parameter_space", "ml_method.pdf"), 
       width = 60, height = 40, units = 'mm', useDingbats=FALSE)
# ggsave(g.filt.method, filename = 
         # here("figures", "parameter_space", "suppl_filt_method.pdf"), 
       # width = 100, height = 100, units = 'mm', useDingbats=FALSE)
ggsave(g.norm.method, filename = 
         here("figures", "parameter_space", "suppl_norm_method.pdf"), 
       width = 100, height = 60, units = 'mm', useDingbats=FALSE)
ggsave(g.fs.cutoff, filename =
         here("figures", "parameter_space", "suppl_fs_cutoff.pdf"),
       width = 100, height = 60, units = 'mm', useDingbats=FALSE)
