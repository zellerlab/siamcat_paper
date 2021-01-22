# ##############################################################################
# function to compute cross-prediction

f.cross.prediction <- function(sc.obj.train, d, c, ctr_type=NULL){

  res.list <- list()
  main.roc <- eval_data(sc.obj.train)$roc
  auroc <- as.numeric(main.roc$auc)
  threshold <- main.roc$thresholds[which(main.roc$specificities > 0.9)[1]]
  case.pred <- enframe(rowMeans(pred_matrix(sc.obj.train)),
                       name='Sample_ID', value='prediction') %>%
    full_join(enframe(sc.obj.train@label$label,
                      name='Sample_ID', value='Group'),
              by='Sample_ID') %>%
    mutate(Group=case_when(Group==-1~'controls', TRUE~'case')) %>%
    filter(Group=='case')
  cross.tasks <- motu.tasks %>%
    mutate(x=paste0(dataset.id, '_', case)) %>%
    filter(x!=paste0(d, '_', c))
  if (!is.null(ctr_type)){
    if (is.character(ctr_type)){
      if (length(ctr_type) > 1){
        stopifnot(all(ctr_type %in% cross.tasks$dataset.id))
        cross.tasks <- cross.tasks %>%
          filter(!dataset.id %in% ctr_type)
      } else if (length(ctr_type) == 1 & ctr_type=='other_ctr'){
        cross.tasks <- cross.tasks %>%
          filter(!dataset.id %in% c('Yachida_2019', 'Qin_2012'))
      }
    } else {
      stop("Wrong input for `ctr_type`!")
    }
  }
  for (j in cross.tasks$x){
    load(here('parameter_space', 'sc',
              paste0('sc_', j, '_mOTUs2.RData')))
    sc.obj <- make.predictions(sc.obj.train, sc.obj, verbose = 0)
    sc.obj <- evaluate.predictions(sc.obj, verbose = 0)
    cross.pred <- enframe(rowMeans(pred_matrix(sc.obj)),
                          name='Sample_ID', value='prediction') %>%
      full_join(enframe(sc.obj@label$label,
                        name='Sample_ID', value='Group'),
                by='Sample_ID') %>%
      mutate(Group=case_when(Group==-1~names(sc.obj@label$info[1]),
                             Group==1~names(sc.obj@label$info[2])))
    cp.auroc.values <- vapply(unique(cross.pred$Group), FUN = function(x){
      temp <- bind_rows(cross.pred %>% filter(Group==x), case.pred)
      res <- auc(response=temp$Group, predictor=temp$prediction,
                 direction='<', levels=c(x, 'case'))
      return(as.double(res))
    }, FUN.VALUE = double(1)) %>%
      enframe(name='Group', value='cp.auroc')

    df.temp <- cross.pred %>%
      group_by(Group) %>%
      summarise(n.group=n(), n.pos=sum(prediction>threshold),
                .groups='drop') %>%
      mutate(frac=n.pos/n.group) %>%
      full_join(cp.auroc.values, by='Group') %>%
      mutate(dataset.id=cross.tasks %>% filter(x==j) %>% pull(dataset.id),
             dataset.test=j,
             dataset.train=paste0(d, '_',c),
             dataset.id.train=d,
             case.train=c,
             auroc=auroc, 
             ext.auroc=as.numeric(sc.obj@eval_data$auroc))
    res.list[[j]] <- df.temp
  }
  res.list <- bind_rows(res.list)
  # add self-predictions
  pred <- enframe(rowMeans(pred_matrix(sc.obj.train)),
                  name='Sample_ID', value='prediction') %>%
    full_join(enframe(sc.obj.train@label$label,
                      name='Sample_ID', value='Group'),
              by='Sample_ID') %>%
    mutate(Group=case_when(Group==-1~names(sc.obj.train@label$info[1]),
                           Group==1~names(sc.obj.train@label$info[2])))
  df.temp <- pred %>%
    group_by(Group) %>%
    summarise(n.group=n(), n.pos=sum(prediction>threshold),
              .groups='drop') %>%
    mutate(frac=n.pos/n.group) %>%
    mutate(cp.auroc=auroc,
           dataset.id=d,
           dataset.test=paste0(d, '_',c),
           dataset.train=paste0(d, '_',c),
           dataset.id.train=d,
           case.train=c,
           auroc=auroc,
           ext.auroc=NA_real_)
  res.list <- bind_rows(res.list, df.temp)
  return(res.list)
}
