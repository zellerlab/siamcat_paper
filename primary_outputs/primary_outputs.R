# ##############################################################################
#
## Primary Outputs
##
## Running the pipeline for the data from Nielsen et al 2014
##  we produce the primary outputs of the package, namely the
##    Association plot
##    Evaluation plot
##    Interpretation plot
#
# ##############################################################################

library("tidyverse")
library("SIAMCAT")
library("curatedMetagenomicData")
library("ggembl")
library("ggpubr")
library("cowplot")
library("here")


# ##############################################################################
# load data

if (!file.exists(here('primary_outputs', 'models_final.RData'))){
  ### Metadata
  meta.nielsen.full <- combined_metadata %>%
    filter(dataset_name=='NielsenHB_2014')
  
  meta.nielsen <- meta.nielsen.full %>% 
    select(sampleID, subjectID, study_condition, disease_subtype,
           disease, age, country, number_reads, median_read_length, BMI) %>%
    mutate(visit=str_extract(sampleID, '_[0-9]+$')) %>%
    mutate(visit=str_remove(visit, '_')) %>%
    mutate(visit=as.numeric(visit)) %>%
    mutate(visit=case_when(is.na(visit)~0, TRUE~visit)) %>%
    group_by(subjectID) %>%
    filter(visit==min(visit)) %>%
    ungroup() %>%
    mutate(Sample_ID=sampleID) %>%
    mutate(Group=case_when(disease=='healthy'~'CTR',
                           TRUE~disease_subtype))
  
  meta.nielsen.uc <- meta.nielsen %>%
    filter(Group %in% c('CTR', 'UC')) %>%
    filter(country == 'ESP') %>%
    as.data.frame()
  rownames(meta.nielsen.uc) <- meta.nielsen.uc$sampleID
  meta.nielsen.uc$sampleID <- NULL
  
  x <- 'NielsenHB_2014.metaphlan_bugs_list.stool'
  feat <- curatedMetagenomicData(x=x, dryrun=FALSE)
  feat <- feat[[x]]@assayData$exprs
  feat <- feat[grep(x=rownames(feat), pattern='s__'),]
  feat <- feat[grep(x=rownames(feat),pattern='t__', invert = TRUE),]
  feat <- t(t(feat)/100)
  rownames(feat) <- str_extract(rownames(feat), 's__.*$')
  
  # ############################################################################
  # SIAMCAT analysis without the Danish control samples
  sc.obj <- siamcat(feat=feat, meta=meta.nielsen.uc, label='Group', case='UC')
  
  sc.obj <- filter.features(sc.obj, cutoff=1e-04,
                            filter.method = 'abundance')
  sc.obj <- filter.features(sc.obj, cutoff=0.05,
                            filter.method='prevalence',
                            feature.type = 'filtered')
  
  sc.obj <- check.associations(sc.obj, detect.lim = 1e-06,
                               alpha=0.1, max.show = 20,
                               plot.type = 'quantile.rect',
                               fn.plot = here('primary_outpus',
                                              'association_plot.pdf'))
  check.confounders(sc.obj, fn.plot = here('primary_outputs', 
                                           'confounders.pdf'))
  
  sc.obj <- normalize.features(sc.obj, norm.method = 'log.std',
                               norm.param = list(log.n0=1e-06, sd.min.q=0))
  sc.obj <- create.data.split(sc.obj, num.folds = 10, num.resample = 10)
  sc.obj <- train.model(sc.obj, method='lasso')
  sc.obj <- make.predictions(sc.obj)
  sc.obj <- evaluate.predictions(sc.obj)
  model.evaluation.plot(sc.obj, fn.plot = here('primary_outputs', 
                                               'eval_plot.pdf'))
  model.interpretation.plot(sc.obj, consens.thres = 0.8,
                            fn.plot = here('primary_outputs', 
                                           'interpretation.pdf'))
  
  # ############################################################################
  # SIAMCAT analysis with the Danish control samples included
  meta.nielsen.uc.dnk <- meta.nielsen %>%
    filter(Group %in% c('CTR', 'UC')) %>%
    # filter(country == 'ESP') %>% # this time we do not remove Danish samples
    as.data.frame()
  rownames(meta.nielsen.uc.dnk) <- meta.nielsen.uc.dnk$sampleID
  meta.nielsen.uc.dnk$sampleID <- NULL
  sc.obj.dnk <- siamcat(feat=feat, meta=meta.nielsen.uc.dnk,
                        label='Group', case='UC')
  sc.obj.dnk <- filter.features(sc.obj.dnk, cutoff=1e-04,
                                filter.method = 'abundance')
  sc.obj.dnk <- filter.features(sc.obj.dnk, cutoff=0.05,
                                filter.method='prevalence',
                                feature.type = 'filtered')
  check.confounders(sc.obj.dnk, fn.plot =  here('primary_outputs', 
                                                'confounders_dnk.pdf'))
  sc.obj.dnk <- check.associations(sc.obj.dnk, detect.lim = 1e-06,
                                   alpha=0.1, max.show = 20,
                                   plot.type = 'quantile.rect',
                                   fn.plot =  here('primary_outputs', 
                                                   'association_plot_dnk.pdf'))
  sc.obj.dnk <- normalize.features(sc.obj.dnk, norm.method = 'log.std',
                                   norm.param = list(log.n0=1e-06, sd.min.q=0))
  sc.obj.dnk <- create.data.split(sc.obj.dnk, num.folds = 10, num.resample = 10)
  sc.obj.dnk <- train.model(sc.obj.dnk, method='lasso')
  sc.obj.dnk <- make.predictions(sc.obj.dnk)
  sc.obj.dnk <- evaluate.predictions(sc.obj.dnk)
  model.evaluation.plot("Spanish samples only"=sc.obj,
                        "Danish and Spanish samples"=sc.obj.dnk,
                        fn.plot = here('primary_outputs', 
                                       'eval_plot_dnk.pdf'))
  
  # ############################################################################
  # SIAMCAT analysis to check if we can  distinguish the Spanish from the
  #    Danish controls?
  meta.nielsen.country <- meta.nielsen %>%
    filter(Group %in% c('CTR')) %>% # only control samples
    # filter(country == 'ESP') %>% # this time we do not remove Danish samples
    as.data.frame()
  rownames(meta.nielsen.country) <- meta.nielsen.country$sampleID
  meta.nielsen.country$sampleID <- NULL
  sc.obj.country <- siamcat(feat=feat, meta=meta.nielsen.country,
                            label='country', case='ESP')
  sc.obj.country <- filter.features(sc.obj.country, cutoff=1e-04,
                                    filter.method = 'abundance')
  sc.obj.country <- filter.features(sc.obj.country, cutoff=0.05,
                                    filter.method='prevalence',
                                    feature.type = 'filtered')
  sc.obj.country <- normalize.features(sc.obj.country, norm.method = 'log.std',
                                       norm.param = list(log.n0=1e-06,
                                                         sd.min.q=0))
  sc.obj.country <- create.data.split(sc.obj.country,
                                      num.folds = 10, num.resample = 10)
  sc.obj.country <- train.model(sc.obj.country, method='lasso')
  sc.obj.country <- make.predictions(sc.obj.country)
  sc.obj.country <- evaluate.predictions(sc.obj.country)
  model.evaluation.plot(sc.obj.country, fn.plot = here('primary_outputs',  
                                                       'eval_plot_country.pdf'))
  save(sc.obj, sc.obj.dnk, sc.obj.country,
       file = here('primary_outputs', 'models_final.RData'))
} else {
  load(here('primary_outputs', 'models_final.RData'))
}


# ##############################################################################
# plot everything in ggplot for more control

# ##############################################################################
# Figure 1

# association plot
temp <- associations(sc.obj)
temp$species <- rownames(temp)
temp <- temp %>%
  mutate(species_short=str_remove(species, '.*s__')) %>%
  filter(p.adj < 0.1) %>%
  arrange(fc) %>%
  mutate(species=factor(species, levels = species))

g.fc <- temp %>%
  ggplot(aes(x=species, y=fc, fill=fc > 0)) +
    geom_bar(stat='identity') +
    coord_flip() +
    theme_publication() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y=element_blank()) +
    ylab('gFC') + xlab('') +
    scale_fill_manual(values=rev(unique(temp$bcol)), guide=FALSE)
g.pval <- temp %>%
  ggplot(aes(x=species, y=-log10(p.adj))) +
    geom_bar(stat='identity') +
    coord_flip() +
    theme_publication() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y=element_blank()) +
    ylab('-') + xlab('')  +
    geom_hline(yintercept = -log10(0.1), colour='red')

df.label <- enframe(sc.obj@label$label, name='SampleID', value="Label")
mat.ab <- get.filt_feat.matrix(sc.obj)[as.character(temp$species),] %>%
  as.data.frame()
mat.ab$species <- rownames(mat.ab)

df.plot <- mat.ab %>%
  mutate(species=str_remove(species, '.*s__')) %>%
  gather(key='SampleID', value="abundance", -species) %>%
  mutate(abundance=log10(abundance+1e-06)) %>%
  mutate(species=factor(species, levels = temp$species_short))


g.ab <- df.plot %>%
  full_join(df.label) %>%
  mutate(Label=as.factor(Label)) %>%
  mutate(species=factor(species, levels = temp$species_short)) %>%
  ggplot(aes(x=species, y=abundance, fill=Label)) +
    geom_boxplot(outlier.shape = NA, size=0.1) +
    geom_jitter(aes(col=Label),
                position = position_jitterdodge(jitter.width = 0.08),
                size=0.6, stroke=0) +
    coord_flip() +
    theme_publication() +
    scale_fill_manual(values=unique(temp$bcol), guide=FALSE) +
    scale_colour_manual(values=str_remove(unique(temp$bcol), '85$'),
                        guide=FALSE) +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank()) +
    ylab('log10(rel. ab.)') + xlab('')

g.all <- plot_grid(g.ab, g.pval, g.fc, nrow=1, rel_widths = c(0.5, 0.25, 0.25))
ggsave(g.all, filename = here('figures', 'primary_output', 'assoc_plot.pdf'),
       useDingbats=FALSE, width = 100, height = 60, units = 'mm')

# confounder
g.bmi <- meta(sc.obj) %>%
  ggplot(aes(x=study_condition, y=BMI)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.08, stroke=0, size=1.3) +
    theme_publication() +
    xlab('') + ylab('BMI')
ggsave(g.bmi, filename = here('figures', 'primary_output', 'conf_bmi.pdf'),
       useDingbats=FALSE, width = 50, height = 60, units = 'mm')

# interpretation plot
feat.weights <- feature_weights(sc.obj)
feat.weights$species <- rownames(feat.weights)
rownames(feat.weights) <- NULL
feat.weight.selected <- feat.weights %>%
  mutate(species=str_remove(species, '.*s__')) %>%
  filter(percentage > 0.8) %>%
  arrange(median.rel.weight)


g.weights <- feat.weight.selected %>%
  mutate(species=factor(species, levels = rev(species))) %>%
  ggplot(aes(x=species, y=-median.rel.weight)) +
    geom_bar(stat='identity') +
    theme_publication() +
    coord_flip() +
    xlab('') + ylab('median rel. feat weight') +
    theme(axis.ticks.y=element_blank(),
          axis.text.y=element_blank()) +
    geom_text(aes(label=paste0(percentage*100, '%'), y=0.08))
ggsave(g.weights, filename = here('figures', 'primary_output', 'weights.pdf'),
       useDingbats=FALSE, width = 50, height = 60, units = 'mm')

# evaluation plot
model.evaluation.plot(sc.obj,
                      fn.plot = here('figures', 'primary_output', 'roc.pdf'))


# ##############################################################################
# Figure 2

# variance plot
feat <- get.filt_feat.matrix(sc.obj.dnk)
label <- label(sc.obj.dnk)$label
country <- as.factor(meta(sc.obj.dnk)$country)
names(country) <- rownames(meta(sc.obj.dnk))

var.label <- vapply(rownames(feat), FUN=function(x){
  x <- feat[x,]
  x <- rank(x)/length(x)
  ss.tot <- sum((x - mean(x))^2)/length(x)
  ss.o.i <- sum(vapply(unique(label), function(s){
    sum((x[label==s] - mean(x[label==s]))^2)
  }, FUN.VALUE = double(1)))/length(x)
  return(1-ss.o.i/ss.tot)
}, FUN.VALUE = double(1))
if (any(is.infinite(var.label))){
  var.label[is.infinite(var.label)] <- NA
}
var.batch <- vapply(rownames(feat), FUN=function(x){
  x <- feat[x,names(country)]
  x <- rank(x)/length(x)
  ss.tot <- sum((x - mean(x))^2)/length(x)
  ss.o.i <- sum(vapply(levels(country), function(s){
    sum((x[country==s] - mean(x[country==s]))^2)
  }, FUN.VALUE = double(1)))/length(x)
  return(1-ss.o.i/ss.tot)
}, FUN.VALUE = double(1))
if (any(is.infinite(var.batch))){
  var.batch[is.infinite(var.batch)] <- NA
}
df.plot <- tibble(label=var.label, batch=var.batch, species=names(var.label))
temp <- get.filt_feat.matrix(sc.obj.dnk)
df.plot$mean <- rowMeans(log10(temp+1e-05))
df.plot <- df.plot %>%
  mutate(species=str_remove(species, '^.*s__'))
g.conf <- df.plot %>%
  ggplot(aes(x=label, y=batch, size=mean)) +
  geom_point(stroke=0, alpha=0.2) +
  theme_publication()  +
  theme(aspect.ratio = 1) +
  xlim(0, 0.3) + ylim(0, 0.3) +
  scale_size_continuous(range = c(0.2, 4))
ggsave(g.conf, filename = here('figures', 'confounders', 'conf_var.pdf'),
       useDingbats=FALSE, width = 80, height = 80, units = 'mm')

# evaluation
model.evaluation.plot(sc.obj.dnk,
                      fn.plot = here('figures', 'confounders',
                                     'roc_all.pdf'))
model.evaluation.plot(sc.obj.country,
                      fn.plot = here('figures', 'confounders',
                                     'roc_country.pdf'))

# associations scatter plot
temp <- associations(sc.obj)
temp$species <- rownames(temp)
temp2 <- associations(sc.obj.dnk)
temp2$species <- rownames(temp2)
df.plot <- full_join(temp, temp2, by='species')

g.assoc <- df.plot %>%
  ggplot(aes(x=-log10(p.adj.x), y=-log10(p.adj.y))) +
    geom_abline(intercept = 0, slope = 1, col='darkgrey', linetype=3) +
    geom_point(alpha=0.2, stroke=0) +
    theme_publication() +
    xlab('Spanish samples only\n-log10(q val)') +
    ylab('Spanish and Danish samples\n-log10(q val)')
ggsave(g.assoc, filename = here('figures', 'confounders', 'assoc_esp_dnk.pdf'),
       width = 40, height = 55, useDingbats=FALSE, units = 'mm')

# Dorea plot
temp <- get.filt_feat.matrix(sc.obj.dnk)
label <- label(sc.obj.dnk)$label
country <- meta(sc.obj.dnk)[['country']]
x <- which(str_detect(rownames(temp), 'Dorea_formicigenerans'))


df.plot <- tibble(bug=log10(temp[x,names(label)] + 1e-05),
                  label=label, country=country) %>%
  mutate(x_value=paste0(country, '_', label))

g <- df.plot %>%
  ggplot(aes(x=x_value, y=bug)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.08, stroke=0, alpha=0.2) +
    theme_publication() +
    xlab('') +
    ylab("log10(Dorea_formicigenerans)") +
    stat_compare_means(comparisons = list(c('DNK_-1', 'ESP_-1'),
                                          c('DNK_-1', 'ESP_1'),
                                          c('ESP_-1', 'ESP_1')))
# right Wilcoxon test
wilcox.test(df.plot$bug~
            df.plot %>% 
              mutate(x_value=case_when(x_value=='ESP_-1'~'DNK_-1',
                                       TRUE~x_value)) %>% 
              pull(x_value))

ggsave(g, filename = here('figures', 'confounders', 'conf_dorea.pdf'),
       width = 38, height = 55, useDingbats=FALSE, units = 'mm')
