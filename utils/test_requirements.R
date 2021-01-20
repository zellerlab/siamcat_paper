packages <- c('tidyverse', 'SIAMCAT', 'pROC', 'ggembl', 'progress', 'here',
              'readxl', 'feather', 'curatedMetagenomicData', 'matrixStats',
              'ggpubr', 'vegan', 'reshape2', 'dplyr', 'ggplot2', 'flexclust',
              'GGally', 'ggrepel', 'gridExtra', 'cowplot', 'yaml', 'optparse',
              'car', 'ggthemes', 'filesstring')

warn <- c()
for (p in packages){
  x <- suppressMessages(suppressWarnings(require(p, character.only = TRUE)))
  if (!x) warn <- c(warn, p)
}

if (length(warn) > 1){
  message('The packages\n\n\t',
          paste0(warn, collapse = ', '),
          '\n\nare not installed!')
} else if (length(warn) == 1){
  message('The package\n\n\t',
          paste0(warn, collapse = ', '),
          '\n\nis not installed!')
} else {
  message('You are good to go! Yay')
}
