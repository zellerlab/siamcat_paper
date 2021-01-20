# ##############################################################################
#
## Reproduce metformin stuff with newer mOTU profiles
#
# ##############################################################################

library("tidyverse")
library("SIAMCAT")
library("ggembl")
library("ggpubr")
library("here")


fn.meta <- here("metformin", "meta_metf_t2d.tsv")
meta <- read.table(fn.meta, sep='\t', stringsAsFactors = FALSE,
                   check.names = FALSE, quote = '')
meta <- meta[meta$GROUP!='IGT',]
meta$METFORMIN[meta$METFORMIN!='Metf'] <- 'NoMetf'

# features
feat <- list()
for (i in c('Karlsson_2013', 'Qin_2012', 'metaHIT')){
  fn.feat <- here('data', 'features', 'motus', paste0(i, '_motus.tsv'))
  x <- read.table(fn.feat, sep='\t', check.names = FALSE,
                  stringsAsFactors = FALSE, quote = '',
                  comment.char = '',
                  header=TRUE) %>% 
    as.matrix()
  # adjust names for the CN cases
  if (i=='Qin_2012') {
    colnames(x) <- str_remove(colnames(x),'bgi-')
  }
  
  feat[[i]] <- x
}

feat <- do.call(cbind, feat)
feat <- feat[,colSums(feat) > 100]

feat <- prop.table(feat, 2)
print(dim(meta))
length(intersect(rownames(meta), colnames(feat)))
meta %>% 
  as_tibble(rownames = 'Sample_ID') %>% 
  filter(Sample_ID %in% colnames(feat)) %>% 
  group_by(COUNTRY) %>% 
  tally()

sc.obj <- siamcat(feat=feat, meta=meta, label='GROUP', case='T2D')

## Filtering
hist(log10(matrixStats::rowMaxs(feat)), 100)
sc.obj <- filter.features(sc.obj, cutoff=1e-03, 
                          filter.method = 'abundance', verbose=3)
hist(rowMeans(get.filt_feat.matrix(sc.obj) != 0), 100)
sc.obj <- filter.features(sc.obj, cutoff=0.05,
                          filter.method='prevalence', verbose=3, 
                          feature.type = 'filtered')
check.confounders(sc.obj, 
                  fn.plot = here('figures', 'metformin', 'confounder_plot.pdf'), 
                  verbose = 3)

feat.filt <- get.filt_feat.matrix(sc.obj)
label <- label(sc.obj)

var.label <- vapply(rownames(feat.filt), FUN=function(x){
  x <- feat.filt[x,]
  x <- rank(x)/length(x)
  ss.tot <- sum((x - mean(x))^2)/length(x)
  ss.o.i <- sum(vapply(unique(label$label), function(s){
    sum((x[label$label==s] - mean(x[label$label==s]))^2)
  }, FUN.VALUE = double(1)))/length(x)
  return(1-ss.o.i/ss.tot)
}, FUN.VALUE = double(1))
if (any(is.infinite(var.label))){
  var.label[is.infinite(var.label)] <- NA
}

meta.temp <- meta(sc.obj)

temp <- factor(meta.temp[['METFORMIN']])
names(temp) <- rownames(meta.temp)
temp <- temp[!is.na(temp)]
var.batch <- vapply(rownames(feat.filt), FUN=function(x){
  x <- feat.filt[x,names(temp)]
  x <- rank(x)/length(x)
  ss.tot <- sum((x - mean(x))^2)/length(x)
  ss.o.i <- sum(vapply(levels(temp), function(s){
    sum((x[temp==s] - mean(x[temp==s]))^2)
  }, FUN.VALUE = double(1)))/length(x)
  return(1-ss.o.i/ss.tot)
}, FUN.VALUE = double(1))
if (any(is.infinite(var.batch))){
  var.batch[is.infinite(var.batch)] <- NA
}
lim <- round(max(var.label, var.batch, na.rm=TRUE), digits = 2)
plot(var.label, var.batch, 
     xlab='Variance explained by label',
     xlim=c(0,lim), ylim=c(0,lim))

df.plot <- tibble(label=var.label,
                  batch=var.batch,
                  species=names(var.label))
df.plot$mean <- rowMeans(log10(feat.filt+1e-05))
df.plot$mean.size <- df.plot$mean+5
df.plot$mean.size <- df.plot$mean.size * 8/5 + 1

plot(var.label, var.batch, type='n', 
     xlab='Variance explained by label',
     ylab='Variance explained by metformin',
     xlim=c(0,lim), ylim=c(0,lim))
symbols(x=var.label, y=var.batch, circles=df.plot$mean.size, inches=1/3,
        bg=alpha("darkgrey", 0.4), fg='black', add=TRUE)
abline(0, 1, col='black', lty=3)


g1 <- df.plot %>% 
  ggplot(aes(x=label, y=batch, size=mean)) + 
    geom_abline(slope = 1, intercept = 0, col='darkgrey', linetype=3) +
    geom_point(alpha=0.4, stroke=0, col='#707372') + 
    theme_publication() + 
    scale_size_continuous(range = c(0.2, 4)) + 
    xlab('Variance explained by Group') + 
    ylab(paste0('Variance explained by Metformin')) + 
    xlim(0,round(max(var.label, var.batch), digits = 2)) + 
    ylim(0,round(max(var.label, var.batch), digits = 2)) +
    theme(plot.background = element_blank(), aspect.ratio = 1)

ggsave(g1, filename = here("figures", "metformin", 'metf_variance.pdf'),
       useDingbats=FALSE, width = 80, height = 80, units = 'mm')


x <- 'Proteobacteria sp. [ref_mOTU_v25_00095]'
label <- label(sc.obj)$label
my_comparisons <- list( c("ND control", "T2D metformin+"), 
                        c("ND control", "T2D metformin-"), 
                        c("T2D metformin+", "T2D metformin-") )
g <- meta(sc.obj) %>% 
  as_tibble() %>% 
  mutate(species=log10(feat.filt[x,names(label)] + 1e-05)) %>% 
  mutate(label=as.factor(label)) %>% 
  mutate(x_value=case_when(label==-1~'ND control',
                           label==1 & METFORMIN=='Metf'~'T2D metformin+',
                           TRUE~'T2D metformin-')) %>% 
  ggplot(aes(x=x_value, y=species, fill=x_value)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15), 
              aes(col=x_value), size=1) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values=c('#226a6575', '#e4943f75', '#51191d75'), 
                    guide=FALSE) + 
  scale_color_manual(values=c('#226a65', '#e4943f', '#51191d'), guide=FALSE) + 
  theme_publication() + 
  xlab('') + 
  ylab('Enterobacteriaceae sp.\nlog10(rel. ab.)') + 
  theme(plot.background = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1)) + 
  stat_compare_means(method='wilcox.test', comparisons = my_comparisons)
ggsave(g, filename = here("figures", 'metformin', "metf_ecoli.pdf"),
       useDingbats=FALSE, width = 50, height = 80, units = 'mm')

# ##############################################################################
# Metformin treatment classification

meta$GROUP2 <- paste0(meta$GROUP, '-', meta$METFORMIN)
comparisons <- list(c('T2D-NoMetf', 'CTR-NoMetf'), # T2D(metf-) vs ND controls
                    c('T2D-Metf', 'CTR-NoMetf'),   # T2D(metf+) vs ND controls
                    c('T2D-NoMetf', 'T2D-Metf'))   # T2D(metf-) vs T2D(metf+)

sc.list <- list()

temp <- read_tsv(here('data', 'motus_taxonomy.tsv'))
df.temp <- as.matrix(temp)
rownames(df.temp) <- df.temp[,2]
feat.temp <- feat
feat.temp <- feat.temp[-which(rownames(feat.temp) == '-1'),]
rownames(feat.temp) <- str_extract(rownames(feat.temp), 
                                   '(ref|meta)_mOTU_v25_[0-9]{5}')
sc.temp <- siamcat(feat=feat.temp)
tax_table(physeq(sc.temp)) <- df.temp
sc.temp <- summarize.features(sc.temp, level='genus')
feat.genus <- get.orig_feat.matrix(sc.temp)

for (x in comparisons){
  label <- create.label(label='GROUP2', case=x[1], control=x[2], meta=meta)
  sc.obj <- siamcat(feat=feat.genus, meta=meta, label=label)
  sc.obj <- filter.features(sc.obj, cutoff=1e-03, 
                            filter.method = 'abundance', verbose=3)
  sc.obj <- filter.features(sc.obj, cutoff=0.05,
                            filter.method='prevalence', verbose=3, 
                            feature.type = 'filtered')
  sc.obj <- normalize.features(sc.obj, norm.method = 'log.std', 
                               norm.param = list(log.n0=1e-05, sd.min.q=0))
  sc.obj <- create.data.split(sc.obj, num.folds = 5, num.resample = 10)
  sc.obj <- train.model(sc.obj, method='lasso')
  sc.obj <- make.predictions(sc.obj)
  sc.obj <- evaluate.predictions(sc.obj)
  sc.list[[paste(x, collapse = ':')]] <- sc.obj
}

model.evaluation.plot('T2D(metf-) vs CTR'=sc.list[[1]], 
                      'T2D(metf+) vs CTR'=sc.list[[2]], 
                      'T2D(metf-) vs T2D(metf+)'=sc.list[[3]],
                      fn.plot = here('figures', 'metformin', 'metf_aurocs.pdf'))
