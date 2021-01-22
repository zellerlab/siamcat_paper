# ##############################################################################
#
## What variance is explained by the different hyperparameter options
##    in the ML workflow?
#
# ##############################################################################

library("tidyverse")
library("here")
library("ggembl")
library("car")


# load final_results
final_results <- read_tsv(here("parameter_space", "job_info",
                               "full_results.tsv"))

all.tasks <- read_tsv(here("parameter_space", "data_info", "all_tasks.tsv"))

out_list <- list()
for (i in seq_len(nrow(all.tasks))){

  # Remove rows where auroc is NA
  tmp <- final_results %>%
    filter(type==all.tasks$type[i],
           dataset.id==all.tasks$dataset.id[i],
           case==all.tasks$case[i]) %>%
    filter(!is.na(auroc))

  # Make sure all independent variables are correctly encoded and the have the
  # correct data type
  tmp <- tmp %>%
    mutate(filt.method=as.factor(filt.method)) %>%
    mutate(norm.method=as.factor(norm.method)) %>%
    mutate(ml.method=as.factor(ml.method)) %>%
    mutate(fs.method=case_when(is.na(fs.method)~'None', TRUE~fs.method)) %>%
    mutate(fs.method=as.factor(fs.method)) %>%
    mutate(fs.cutoff=case_when(is.na(fs.cutoff)~0, TRUE~fs.cutoff)) %>%
    mutate(fs.cutoff=as.factor(fs.cutoff))

  # There is still one pair of singular variables:
  #   fs.cutoff0 <--> fs.methodNone

  mylm <- lm(formula = auroc ~ filt.method + norm.method +
               ml.method + fs.method + fs.cutoff, data = tmp)
  # 1. Calcualte the sequential sum of squares as one way to calculate the
  #   proportion of variance explained (type I, standard for anova R).
  #   "How much MORE variance can be explained by adding B to a model that
  #   already contains A?"
  anova_typeI <- anova(mylm)
  # 2. Calcualte the sequential sum of squares as one way to calculate the
  #   proportion of variance explained (type III, standard for anova in SPSS)
  #   Given a model that contains every variable except A, how much MORE
  #   variance can be explained by a model containing all those variables
  #   and A?" For further info regaring the difference between those two
  #   ways to determine the proportion of variance explained, see this
  #   excellent stack overflow post:
  #   https://stats.stackexchange.com/questions/20452/how-to-interpret-type-i-type-ii-and-type-iii-anova-and-manova
  anova_typeIII <- Anova(mylm, type = "III", singular.ok = TRUE)
  # 3. Train a simple linear regression to predictAUC from a single variable
  #   and calculate the explained sum of squares. Train one model per variable.
  #   Compare the individual explained sum of squares to the explained sum of
  #   squares for the model containing all variables. Will not add up, since
  #   variables are not all orthogonal to each other.
  vars <- c("filt.method", "norm.method", "ml.method",
            "fs.method", "fs.cutoff")
  anova_single_feature <- lapply(vars, function(x){
    mylm <- lm(as.formula(
            paste("auroc",
            x,
            sep = " ~ ")), data = tmp)
    anova_sf <- anova(mylm)
    return(anova_sf[, 2, drop = FALSE])
  })

  names(anova_single_feature) <- vars

  out_list[[i]] <- list(
    "anova_typeI" = anova_typeI,
    "anova_typeIII" = anova_typeIII,
    "anova_single_feature" = anova_single_feature)
  names(out_list)[i] <- all.tasks %>% select(1:3) %>%
    slice(i) %>% paste(collapse = '-')
}



# get results out of the list
all.res <- list()
all.res[['anova_typeI']] <- map(
  map(out_list, "anova_typeI"), .f = function(x){
    tmp <- (x$`Sum Sq`/sum(x$`Sum Sq`))*100
    names(tmp) <- rownames(x)
    tmp %>% enframe()
    }) %>%
  bind_rows(.id='dataset')

all.res[['anova_typeIII']] <- map(
  map(out_list, "anova_typeIII"), .f = function(x){
    tmp <- (x$`Sum Sq`/sum(x$`Sum Sq`))*100
    names(tmp) <- rownames(x)
    tmp %>% enframe()
  }) %>%
  bind_rows(.id='dataset')

all.res[['anova_single_feature']] <- map(
  map(out_list, "anova_single_feature"), .f = function(x){
    tmp <- map(x, .f = function(y){((y/sum(y))*100)[1,]}) %>% bind_rows()
    tmp <- tibble(name=colnames(tmp), value=t(tmp)[,1])
  }) %>%
  bind_rows(.id='dataset')

all.res <- bind_rows(all.res, .id='method') %>%
  filter(name!='(Intercept)')

# plot
refactor <- all.res %>%
  group_by(name) %>%
  summarise(x=median(value)) %>%
  arrange(x)
all.res <- all.res %>%
  mutate(name=factor(name, levels = refactor$name))

g <- all.res %>%
  ggplot(aes(x = name, y = value)) +
    geom_boxplot(outlier.shape = NA, col='#707372') +
    geom_jitter(colour='#70737275', width = 0.1, stroke=0) +
    theme_publication() +
    ylab("Proportion of variance") + xlab("") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_grid(~method)

ggsave(g, filename = here("figures", "parameter_space",
                          "variance_explained.pdf"),
       height = 70, width = 140, units = 'mm', useDingbats=FALSE)
