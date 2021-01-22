# ##############################################################################
#
## Function to train control-augmented models
#
# ##############################################################################

library("progress")

train.model.ctr <- function(siamcat,
                            correct.list,
                            method = c("lasso",
                                       "enet",
                                       "enet-0.5"),
                            stratify = TRUE,
                            modsel.crit = list("auc"),
                            min.nonzero.coeff = 1,
                            param.set = NULL,
                            perform.fs = FALSE,
                            param.fs = list(thres.fs = 100, method.fs = "AUC",
                                            direction="absolute"),
                            feature.type='normalized',
                            n=2,
                            verbose = 1) {
  
  if (verbose > 1)
    message("+ starting train.model")
  s.time <- proc.time()[3]
  
  # check and get features
  if (!feature.type %in% c('original', 'filtered', 'normalized')){
    stop("Unrecognised feature type, exiting...\n")
  }
  
  feat.ctr <- do.call(cbind, correct.list)
  meta.ctr <- map(seq_along(correct.list),
                  .f = function(x){tibble(
                    Sample_ID=colnames(correct.list[[x]]),
                    Study=paste0('Study_', x))}) %>%
    bind_rows()
  
  if (feature.type == 'original'){
    feat <- get.orig_feat.matrix(siamcat)
    feat.ctr <- feat.ctr[rownames(feat),]
  } else if (feature.type == 'filtered'){
    if (is.null(filt_feat(siamcat, verbose=0))){
      stop('Features have not yet been filtered, exiting...\n')
    }
    feat <- get.filt_feat.matrix(siamcat)
    feat.ctr <- feat.ctr[rownames(feat),]
  } else if (feature.type == 'normalized'){
    if (is.null(norm_feat(siamcat, verbose=0))){
      stop('Features have not yet been normalized, exiting...\n')
    }
    feat <- get.norm_feat.matrix(siamcat)
    # frozen normalization
    temp <- norm_params(siamcat, verbose=0)
    feat.ctr <- feat.ctr[temp$retained.feat,]
    feat.ctr.log <- log10(feat.ctr + temp$log.n0)
    feat.ctr.norm <- (feat.ctr.log - temp$feat.mean)/temp$feat.adj.sd
    feat.ctr <- feat.ctr.norm
  }
  rownames(feat.ctr) <- make.names(rownames(feat.ctr))
  
  # make sure the names fit
  rownames(feat) <- make.names(rownames(feat))
  
  
  # checks
  label <- label(siamcat)
  if (label$type == "TEST"){
    stop('Model can not be trained to SIAMCAT object with a TEST label.',
         ' Exiting...')
  }
  if (is.null(data_split(siamcat, verbose=0))){
    stop("SIAMCAT object needs a data split for model training! Exiting...")
  }
  data.split <- data_split(siamcat)
  
  # check modsel.crit
  if (!all(modsel.crit %in% c("auc", "f1", "acc", "pr", "auprc"))) {
    warning("Unkown model selection criterion... Defaulting to AU-ROC!\n")
    measure <- list(mlr::auc)
  } else {
    measure <- list()
  }
  if (verbose > 2)
    message("+++ preparing selection measures")
  for (m in modsel.crit) {
    if (m == "auc") {
      measure[[length(measure) + 1]] <- mlr::auc
    } else if (m == "acc") {
      measure[[length(measure) + 1]] <- mlr::acc
    } else if (m == "f1") {
      measure[[length(measure) + 1]] <- mlr::f1
    } else if (m == "pr" || m == "auprc") {
      auprc <- makeMeasure(
        id = "auprc",
        minimize = FALSE,
        best = 1,
        worst = 0,
        properties = c("classif", "req.pred",
                       "req.truth", "req.prob"),
        name = "Area under the Precision
        Recall Curve",
        fun = function(task,
                       model,
                       pred,
                       feats,
                       extra.args) {
          measureAUPRC(
            getPredictionProbabilities(pred),
            pred$data$truth,
            pred$task.desc$negative,
            pred$task.desc$positive
          )
        }
      )
      measure[[length(measure) + 1]] <- auprc
    }
  }
  
  # Create matrix with hyper parameters.
  hyperpar.list <- list()
  
  # Create List to save models.
  models.list <- list()
  num.runs <- data.split$num.folds * data.split$num.resample
  bar <- 0
  if (verbose > 1)
    message(paste("+ training", method, "models on", num.runs,
                  "training sets"))
  
  if (verbose > 1 & perform.fs){
    message('+ Performing feature selection ',
            'with following parameters:')
    for (i in seq_along(param.fs)) {
      message(paste0('    ', names(param.fs)[i], ' = ',
                     ifelse(is.null(param.fs[[i]]), 'NULL', param.fs[[i]])))
    }
  }
  
  if (verbose == 1 || verbose == 2)
    pb <- progress_bar$new(total = num.runs)
  
  for (fold in seq_len(data.split$num.folds)) {
    if (verbose > 2)
      message(paste("+++ training on cv fold:", fold))
    
    for (resampling in seq_len(data.split$num.resample)) {
      if (verbose > 2)
        message(paste("++++ repetition:", resampling))
      
      fold.name <-
        paste0("cv_fold",
               as.character(fold),
               "_rep",
               as.character(resampling))
      fold.exm.idx <- match(data.split$
                              training.folds[[resampling]][[fold]],
                            names(label$label))
      
      ### subselect examples for training
      label.fac <-
        factor(label$label,
               levels = sort(label$info))
      train.label <- label.fac[fold.exm.idx]
      data <-
        as.data.frame(t(feat)[fold.exm.idx,])
      
      # add controls
      n.ctr <- sum(train.label==-1)
      n.ctr.add <- n.ctr*n
      
      n.ctr.add.indiv <- n.ctr.add/length(correct.list)
      ctr.add <- c()
      meta.ctr.temp <- meta.ctr
      
      while(TRUE){
        if (any(table(meta.ctr.temp$Study) < n.ctr.add.indiv)){
          x <- sort(table(meta.ctr.temp$Study))
          add.2 <- meta.ctr.temp %>%
            filter(Study==names(x)[1]) %>%
            pull(Sample_ID)
          ctr.add <- c(ctr.add, add.2)
          meta.ctr.temp <- meta.ctr.temp %>%
            filter(Study!=names(x)[1])
          n.ctr.add.indiv <- (n.ctr.add-length(add.2))/
            length(unique(meta.ctr.temp$Study))
          if (nrow(meta.ctr.temp) ==0){
            break()
          }
        } else {
          ctr.add <- c(ctr.add, meta.ctr.temp %>%
                         group_by(Study) %>%
                         sample_n(size=n.ctr.add.indiv) %>%
                         pull(Sample_ID))
          break()
        }
      }
      
      
      data.ctr.add <- feat.ctr[,ctr.add]
      ctr.label <- rep(-1, ncol(data.ctr.add))
      names(ctr.label) <- colnames(data.ctr.add)
      
      data <- rbind(data, t(data.ctr.add))
      train.label <- c(train.label, factor(ctr.label, levels=c(-1, 1))) - 2
      train.label[train.label==0] <- 1
      train.label <- factor(train.label, levels = c(-1, 1))
      
      stopifnot(nrow(data) == length(train.label))
      stopifnot(all(rownames(data) == names(train.label)))
      
      #feature selection
      if (perform.fs) {
        
        stopifnot(all(c('method.fs', 'thres.fs', 'direction') %in%
                        names(param.fs)))
        
        # test method.fs
        if (!param.fs$method.fs %in% c('Wilcoxon', 'AUC', 'gFC')) {
          stop('Unrecognised feature selection method...\n')
        }
        
        # assert the threshold
        if (param.fs$method.fs == 'Wilcoxon') {
          stopifnot(param.fs$thres.fs < 1 && param.fs$thres.fs > 0)
        } else {
          stopifnot(param.fs$thres.fs > 10)
        }
        stopifnot(param.fs$thres.fs < ncol(data))
        
        if (param.fs$method.fs == 'Wilcoxon') {
          assoc <- vapply(data,
                          FUN=function(x, label){
                            d <- data.frame(x=x, y=label);
                            t <- wilcox.test(x~y, data=d)
                            return(t$p.val)
                          },
                          FUN.VALUE=double(1),
                          label=train.label)
          data <- data[,which(assoc < param.fs$thres.fs)]
        } else if (param.fs$method.fs == 'AUC') {
          assoc <- vapply(data,
                          FUN=SIAMCAT:::get.single.feat.AUC,
                          FUN.VALUE = double(1),
                          label=train.label,
                          pos=max(label$info),
                          neg=min(label$info))
          if (param.fs$direction == 'absolute'){
            assoc[assoc < 0.5] <- 1 - assoc[assoc < 0.5]
          } else if (param.fs$direction == 'negative'){
            assoc <- 1 - assoc
          }
          asso <- assoc[assoc > 0.5]
          data <- data[,names(which(
            rank(-assoc) <= param.fs$thres.fs))]
        } else if (param.fs$method.fs == 'gFC') {
          assoc <- vapply(data,
                          FUN=get.quantile.FC,
                          FUN.VALUE = double(1),
                          label=train.label,
                          pos=max(label$info),
                          neg=min(label$info))
          if (param.fs$direction == 'absolute'){
            assoc <- abs(assoc)
          } else if (param.fs$direction == 'negative'){
            assoc <- -assoc
          }
          assoc <- assoc[assoc > 0]
          data <- data[,names(which(
            rank(-assoc) <= param.fs$thres.fs))]
          
        }
        
        stopifnot(ncol(data) > 0)
        if (verbose > 2) {
          message(paste0('++ retaining ', ncol(data),
                         ' features after feature selection with ',
                         param.fs$method.fs, '-threshold ',
                         param.fs$thres.fs))
        }
      }
      
      data$label <- train.label
      
      ### internal cross-validation for model selection
      if (method=='enet-0.5'){
        model <-
          SIAMCAT:::train.plm(
            data = data,
            method = 'enet',
            measure = measure,
            min.nonzero.coeff = min.nonzero.coeff,
            param.set = list(alpha=0.5),
            neg.lab = min(label$info)
          )
      } else {
        model <-
          SIAMCAT:::train.plm(
            data = data,
            method = method,
            measure = measure,
            min.nonzero.coeff = min.nonzero.coeff,
            param.set = param.set,
            neg.lab = min(label$info)
          )
      }
      
      bar <- bar + 1
      
      if (!all(model$feat.weights == 0)) {
        models.list[[bar]] <- model
      } else {
        warning("Model without any features selected!\n")
      }
      
      if (verbose == 1 || verbose == 2)
        pb$tick()
    }
  }
  
  model_list(siamcat) <- list(
    models = models.list,
    model.type = method,
    feature.type = feature.type)
  
  e.time <- proc.time()[3]
  
  if (verbose > 1)
    message(paste(
      "+ finished train.model in",
      formatC(e.time - s.time,
              digits = 3),
      "s"
    ))
  if (verbose == 1)
    message(paste("Trained", method, "models successfully."))
  
  return(siamcat)
}
