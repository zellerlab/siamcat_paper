# ##############################################################################
#
## Script to compare supervised classification by ML to PERMANOVA-style
##    binary classification via within-group-distances
#
# ##############################################################################

library("tidyverse")
library("vegan")
library("here")
library("ggembl")
library("yaml")
library("pROC")
library("ggrepel")

within.vs.between.dists <- function(distance_input, label_input){
  
  distance_input <- as.matrix(distance_input)
  
  stopifnot(all(label_input$Sample_ID %in% rownames(distance_input)))
  
  diag(distance_input) <- NA
  distance_input[lower.tri(distance_input)] <- NA
  
  tbl.tmp <- tibble(Sample_ID_1=rownames(distance_input))
  all.data <- map(as.list(as.data.frame(distance_input)), .f = function(x){
    tbl.tmp %>% mutate(value=x)
  }) %>% bind_rows(.id = 'Sample_ID_2')
  
  all.data <- all.data %>% 
    full_join(label_input %>% 
                transmute(Sample_ID_1=Sample_ID, Group_1=Group),
              by='Sample_ID_1') %>% 
    full_join(label_input %>% 
                transmute(Sample_ID_2=Sample_ID, Group_2=Group),
              by='Sample_ID_2') %>% 
    filter(!is.na(value))
  
  res <- roc(controls=all.data %>% filter(Group_1==Group_2) %>% pull(value),
      cases=all.data %>% filter(Group_1!=Group_2) %>% pull(value),
      direction = '<')
  return(as.numeric(res$auc))
}


colours <- yaml.load_file(here("parameter_space", "data_info", "colours.yaml"))
type_colours <- unlist(colours$type)

# full metadata
meta.all <- read_tsv(here("parameter_space", "data_info", "meta_full.tsv"))

# all ML aurocs
auroc.all <- read_tsv(here("parameter_space", "files", "auroc_all.tsv")) %>% 
  mutate(identity_bray=NA_real_,log_euclidean=NA_real_,
         perm.identity_bray=NA_real_, perm.log_euclidean=NA_real_)

out.list <- list()
for (i in seq_len(nrow(auroc.all))){
  message(auroc.all$dataset.id[i], '-', auroc.all$case[i], '-', 
          auroc.all$type[i])
  log.n0 <- ifelse(auroc.all$type[i] %in% c('eggNOG', 'humann2'), 
                   1e-08, 1e-06)
  feat <- read.table(
    here("parameter_space", "files", "feat", 
         paste0('feat_rel_', auroc.all$dataset.id[i], '_',
                auroc.all$case[i], '_', auroc.all$type[i], 
                '.tsv')),
    sep='\t', stringsAsFactors = FALSE, check.names = FALSE,
    quote = '', comment.char = '', header = TRUE, row.names = 1) %>% 
    as.matrix()
  stopifnot(all(colnames(feat) %in% meta.all$Sample_ID))
  meta.temp <- meta.all %>% 
    filter(Study==auroc.all$dataset.id[i]) %>% 
    filter(Sample_ID%in%colnames(feat))
  nmds.list <- list()
  for (dist in c("identity_bray", "log_euclidean")){
    if (dist == "log_euclidean"){
      # log_Euclidean
      t1 <- vegdist(t(log10(feat + log.n0)), method='euclidean')
    } else if (dist == 'identity_bray'){
      # identity_bray
      t1 <- vegdist(t(feat), method = 'bray')
    } else {
      stop("Can't determine distance.")
    }
    auroc <- within.vs.between.dists(t1, meta.temp)
    tmp <- adonis(formula = t1 ~ Group, data = meta.temp, permutations = 1000)
    p.value <- tmp$aov.tab$`Pr(>F)`[1]
    auroc.all[[dist]][i] <- auroc
    auroc.all[[paste0('perm.', dist)]][i] <- p.value
    # nmds
    nmds <- metaMDS(t1, k=2, trace=0)
    df.plot <- as_tibble(nmds$points, rownames = 'Sample_ID') %>% 
      full_join(meta.temp, by='Sample_ID') %>% 
      mutate(Group=case_when(Group=='CTR'~'CTR', Group=='H'~'CTR', 
                             Group=='control'~'CTR', TRUE~'case')) %>% 
      mutate(method=dist)
    nmds.list[[dist]] <- df.plot
  }
  out.list[[i]] <- bind_rows(nmds.list)
}

names(out.list) <- auroc.all %>% 
  mutate(x=paste0(dataset.nice, '_', case.nice, '_', type)) %>% 
  pull(x)

df.plot <- bind_rows(out.list, .id='ID') %>% 
  separate(ID, into=c("dataset.nice", 'case.nice', 'type'), sep = '_')

highlight.set <- tibble(
  dataset.nice=c('Zhang 2015', 'Feng 2015', 'Qin 2014', 'He 2017',
                 'Schubert 2014', 'Lozupone 2013', 
                 'Turnbaugh 2009', 'Nielsen 2014'),
  type=c('mOTUs2', 'metaphlan', 'humann2', 'mOTUs2', 'RDP',
         'RDP', 'RDP', 'eggNOG'),
  case.nice=c('ART', 'CRC', 'LIV', 'CD', 'CDI', 'HIV', 'OB', 'CD')) %>% 
  left_join(auroc.all %>% transmute(dataset.nice=dataset.nice, type=type, 
                                    case.nice=case.nice, label.red=label.nice),
            by=c('dataset.nice', 'type', 'case.nice'))

# plot everything

for (d in c('identity_bray', 'log_euclidean')){
  g <- auroc.all %>% 
    mutate(n.all=n.ctr+n.case) %>% 
    mutate(border=!!sym(paste0('perm.', d)) < 0.005) %>% 
    full_join(highlight.set, by=c('dataset.nice', 'type', 'case.nice')) %>% 
    mutate(label.red=case_when(is.na(label.red)~'', TRUE~label.red)) %>% 
    ggplot(aes(x=!!sym(d), y=auroc, color = type, size = n.all)) +
      geom_abline(slope = 1, intercept = 0, colour='lightgrey', lty=3) + 
      geom_point(alpha = 0.5) +
      geom_point(data=. %>% filter(border), colour='black', shape=21, 
                 alpha=0.5, fill=NA) +
      scale_size(range = c(0.5, 3)) +
      scale_color_manual(values = type_colours) +
      guides(color=guide_legend(nrow=1,byrow=TRUE)) +
      theme_publication() +
      theme(legend.box = 'vertical', legend.position = 'bottom', 
            legend.direction = 'horizontal', 
            legend.margin=margin(t = -0.3, unit='cm'),
            legend.spacing.x = unit(-0.12, 'cm'), 
            aspect.ratio = 1) +
      scale_y_continuous(limits = c(0.47, 1), breaks = seq(0.5, 1, 0.1)) +
      scale_x_continuous(limits = c(0.47, 1), breaks = seq(0.5, 1, 0.1))  +
      xlab("Within-vs-between group distance AUROC") +
      ylab("Supervised classification AUROC") + 
      geom_text_repel(aes(label=label.red), color = 'black', size = 2.25,
                      nudge_x = 0.13, nudge_y = -0.1) +
      NULL
  
  ggsave(g, filename = here("figures", "parameter_space", 
                            paste0("scatter_", d, ".pdf")),
         height = 100, width = 85, units = 'mm', useDingbats=FALSE)  
  
  # plot individual ordinations
  g <- df.plot %>% 
    full_join(highlight.set, by=c('dataset.nice', 'type', 'case.nice')) %>% 
    filter(!is.na(label.red)) %>% 
    filter(method==d) %>% 
    ggplot(aes(x=MDS1, y=MDS2, col=Group)) +
      geom_point(shape=20, stroke=0) + 
      facet_wrap(~Study, nrow = 2, scales='free') + 
      theme_publication() +
      xlab("Axis 1") + ylab("Axis 2") + 
      theme(axis.text=element_blank(), 
            axis.ticks = element_blank()) + 
      scale_colour_manual(values=c('#D51E4650', '#AAA99D50'), guide=FALSE)
  ggsave(g, filename = here("figures", "parameter_space",
                            paste0("partial_", d, ".pdf")),
         height = 60, width = 85, units = 'mm', useDingbats=FALSE)
      
}
