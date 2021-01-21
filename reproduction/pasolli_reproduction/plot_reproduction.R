# ##############################################################################
#
# Plot AUCs for the reproduction of the results in Pasolli et al
#
# ##############################################################################

# packages
library("tidyverse")
library("ggembl")
library("cowplot")
library("here")


# disease colors from duvallet et al
disease_colors = c('CDI'= "#61AA60", 'EDD'= "#4b847a", 'IBD'= "#996CCE", 
                   'UC'= "#bc78c2", 'CD'= "#7a78c2", 'OB'= "#F0C948", 
                   'CRC'= "#F56484", 'ASD'= "#6992cf", 'T1D'= "#c98746",
                   'NASH'= "#4aac8b", 'LIV'= "#cc436f",  'CIRR'= "#cc436f",
                   'MHE'= "#cc436f", 'HIV'= "#B86958",  'PAR'= "#c07198",  
                   'ART'= "#d59847",  'RA'= "#d59847",  'PSA'= "#d59847")


# result from pasolli et al
results.df <- tibble(Study=c('Qin 2014 LIV', 'Zeller 2014 CRC', 
                             'Nielsen 2014 IBD','Chatelier 2013 OB', 
                             'Qin 2012 T2D', 'Karlsson 2013 T2D'), 
                     Disease=c('LIV', 'CRC', 'IBD', 'OB', 'T1D', 'T1D'), 
                     `AUC-Pasolli` = c(0.945, 0.873, 0.890, 
                                       0.655, 0.744, 0.762),
                     n.samples=c(232, 121, 110, 278, 344, 96))

# ##############################################################################
# add reproduction data
SIAMCAT.results <- tibble(Study=c('Qin 2014 LIV', 'Zeller 2014 CRC', 
                                  'Nielsen 2014 IBD', 'Chatelier 2013 OB', 
                                  'Qin 2012 T2D', 'Karlsson 2013 T2D'), 
                          auc = c(0.9426108, 0.8692922, 0.8908235, 0.6749794, 
                                  0.7429344, 0.7648091),
                          low=c(0.9101852, 0.8003736, 0.8266867, 
                                0.6070811, 0.6914618, 0.6685375),
                          high=c(0.9750364, 0.9382108, 0.9549603, 
                                 0.7428778, 0.794407, 0.8610808))

df.plot <- full_join(results.df, SIAMCAT.results, by='Study')

# ##############################################################################
# plot
g <- df.plot %>% 
  mutate(Study=factor(Study, levels=Study)) %>% 
  ggplot(aes(x=Study, y=`AUC-Pasolli`, col=Disease)) + 
    geom_point(shape=18, size=4) +
    geom_segment(aes(x = as.numeric(Study),
                     xend = as.numeric(Study),
                     y=0.5, yend=`AUC-Pasolli`), linetype=3) + 
    geom_segment(aes(x = as.numeric(Study)+0.25, 
                     xend = as.numeric(Study)+0.25,
                     y=low, yend=high), size=0.75, alpha=0.75) +
    geom_point(aes(x=as.numeric(Study)+0.25, y=auc, fill=Disease), 
               shape=23, size=4, col='black') + 
    coord_flip(ylim = c(0.5, 1), xlim=c(0, nrow(df.plot)+1), expand=FALSE) +
    scale_colour_manual(values=disease_colors, guide=FALSE) + 
    scale_fill_manual(values=disease_colors, guide=FALSE) + 
    ylab('AUC') + 
    xlab('') + 
    theme_embl(font.size = 12, panel.grid = 'major_y') + 
    theme(panel.grid.major.y = element_blank())

g2 <- df.plot %>% 
  mutate(Study=factor(Study, levels=Study)) %>% 
  ggplot(aes(x=Study, y=n.samples, col=Disease)) +
  geom_bar(stat='identity', fill='white') + 
  theme_embl(font.size = 12, panel.grid = 'major_x') + 
  theme(panel.grid.major.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + 
  coord_flip() + 
  scale_color_manual(values=disease_colors, guide=FALSE) + 
  xlab('') + ylab('N samples')

g.all <- plot_grid(g, g2, rel_widths = c(0.75, .25))
ggsave(g.all, 
       filename = here('figures', 'reproduction', 'pasolli_reproduction.pdf'),
       width = 170, height = 50, useDingbats=FALSE, units = 'mm')
