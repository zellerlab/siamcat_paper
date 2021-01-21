# ##############################################################################
#
# Plot AUCs for the reproduction of the results in Duvallet et al
#
# ##############################################################################

# packages
library("tidyverse")
library("ggthemes")
library("cowplot")
library("ggembl")
library("here")

# disease colors from duvallet et al
disease_colors = c('CDI'= "#61AA60", 'EDD'= "#4b847a", 'IBD'= "#996CCE", 
                   'UC'= "#bc78c2", 'CD'= "#7a78c2", 'OB'= "#F0C948", 
                   'CRC'= "#F56484", 'ASD'= "#6992cf", 'T1D'= "#c98746",
                   'NASH'= "#4aac8b", 'LIV'= "#cc436f",  'CIRR'= "#cc436f",
                   'MHE'= "#cc436f", 'HIV'= "#B86958",  'PAR'= "#c07198",  
                   'ART'= "#d59847",  'RA'= "#d59847",  'PSA'= "#d59847")


# result from duvallet et al
results.df <- tibble(Study=c('Singh 2015', 'Schubert 2014', 'Vincent 2013',
                             'Goodrich 2014', 'Turnbaugh 2009', 
                             'Zupancic 2012', 'Ross 2015', 'Zhu 2013', 
                             'Baxter 2016', 'Zeller 2014', 'Wang 2012', 
                             'Chen 2012', 'Gevers 2014', 'Morgan 2012', 
                             'Papa 2012', 'Willing 2010', 
                             'Noguera-Julian 2016', 'Dinh 2015', 
                             'Lozupone 2013', 'Son 2015', 'Kang 2013', 
                             'Alkanani 2015', 'Mejia-Leon 2014', 
                             'Wong 2013', 'Zhu 2013', 'Scher 2013', 
                             'Zhang 2013', 'Scheperjans 2015'), 
                     Disease=c('EDD', 'CDI', 'CDI', 'OB', 'OB', 'OB', 'OB',
                               'OB', 'CRC', 'CRC', 'CRC', 'CRC', 'IBD', 
                               'IBD', 'IBD', 'IBD', 'HIV', 'HIV', 'HIV', 
                               'ASD', 'ASD', 'T1D', 'T1D', 'NASH', 'NASH', 
                               'ART', 'LIV', 'PAR'), 
                     `AUC-Duvallet` = c(0.96, 0.99, 0.91, 0.67, 0.84, 0.44, 
                                        0.49, 0.86, 0.77, 0.82, 0.9, 0.78, 
                                        0.71, 0.81, 0.84, 0.66, 0.67, 0.22, 
                                        0.92, 0.39, 0.76, 0.71, 0.77, 0.68, 
                                        0.93, 0.62, 0.8, 0.67))
# ##############################################################################
# minor modificaiton
results.df <- results.df %>% 
  mutate(x_value=paste(results.df$Study, results.df$Disease)) %>% 
  mutate(x_value=factor(x_value, levels=rev(x_value))) %>% 
  mutate(match=paste0(Disease, '_', str_extract(Study, pattern='^[A-Za-z-]+')))

# ##############################################################################
# add reproduction data

SIAMCAT.results <- read_tsv(here('reproduction', 
                                 'duvallet_reproduction', 
                                 'results.txt')) %>% 
  mutate(Disease=toupper(str_extract(data.tag, '^[a-z0-9]+'))) %>% 
  mutate(Study=str_remove(data.tag, '^[a-z0-9]+_')) %>% 
  mutate(Study=str_remove(Study, '_[A-Z0-9_]+$')) %>% 
  mutate(Study=str_to_title(Study)) %>% 
  mutate(Study=ifelse(Study=='Mejialeon', 'Mejia-Leon', Study),
         Study=ifelse(Study=='Noguerajulian', 'Noguera-Julian', Study)) %>% 
  mutate(match=paste0(Disease, '_', Study))

df.plot <- full_join(results.df, SIAMCAT.results, by='match')

# ##############################################################################
# plot
g <- df.plot %>% 
  ggplot(aes(x=x_value, y=`AUC-Duvallet`, col=Disease.x)) + 
    theme_embl(font.size = 12, panel.grid = 'major_x') + 
    theme(panel.grid.major.y = element_blank()) +
    geom_point(shape=18, size=3) +
    geom_segment(aes(x = as.numeric(x_value),
                 xend = as.numeric(x_value),
                 y=0.5, yend=`AUC-Duvallet`), linetype=3) + 
    geom_segment(aes(x = as.numeric(x_value)+0.25, 
                     xend = as.numeric(x_value)+0.25,
                     y=low, yend=high), size=0.75, alpha=0.75) +
    geom_point(aes(x=as.numeric(x_value)+0.25, y=auc, fill=Disease.x), 
              shape=23, size=3, col='black') + 
    coord_flip(ylim = c(0.5, 1), xlim=c(0, nrow(df.plot)+1), expand=FALSE) +
    scale_colour_manual(values=disease_colors, guide=FALSE) + 
    scale_fill_manual(values=disease_colors, guide=FALSE) + 
    ylab('AUROC') + 
    xlab('')



df.n <- read_delim(here('reproduction', 'duvallet_reproduction', 
                        './table1.dataset_info.md'), 
                   delim = '|') %>% 
  mutate(n.ctr=as.integer(trimws(` N_ctrl `))) %>% 
  mutate(n.dis=as.integer(trimws(` N_dis `))) %>% 
  mutate(data_set=trimws(`dataset_label `)) %>% 
  mutate(n.all=n.ctr+n.dis) %>%
  mutate(data_set=str_remove(data_set, ',')) %>% 
  filter(data_set %in% df.plot$x_value) %>% 
  mutate(disease=str_split_fixed(data_set, pattern = ' ', n=3)[,3])

g2 <-
  df.n %>% 
  mutate(data_set=factor(data_set, levels = levels(df.plot$x_value))) %>%
  ggplot(aes(x=data_set, y=n.all, col=disease)) +
    geom_bar(stat='identity', fill='white') + 
    theme_embl(font.size = 12, panel.grid = 'major_x') + 
    theme(panel.grid.major.y = element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) + 
    coord_flip(ylim=c(0, 150)) + 
    scale_color_manual(values=disease_colors, guide=FALSE) + 
    xlab('') + ylab('N samples') + 
    geom_text(aes(label=n.all, y=155), col='black')
    


g.all <- plot_grid(g, g2, rel_widths = c(0.75, .25))
ggsave(g.all, 
       filename = here('figures', 'reproduction', 'duvallet_reproduction.pdf'),
       width = 170, height = 160, useDingbats=FALSE, units = 'mm')
