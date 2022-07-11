# helper functions for training WID-HPV using 10-fold cross validation
el_classifier <- function(beta_tr, type_tr,
                          beta_val, type_val,
                          alpha = 0.0,
                          lambda = "(-4,4)",
                          slope = "default"){
  
  require(Hmisc)
  require(rms)
  require(pROC)
  require(glmnet)
  
  lambda <- if (lambda == "(-4,4)"){
    exp(seq(-4,4,by=0.1))
  } else {
    if (lambda == "unrestricted"){
      NULL
    } else {
      x <- strsplit(lambda, split = "[(|,|)+]")
      lambda.min <- x[[1]][2]
      lambda.max <- x[[1]][3]
      exp(seq(lambda.min,lambda.max,by=0.1))
    }
  }
  
  
  # create training and validation datasets
  dat_tr <- list(x = beta_tr,
                 y = type_tr)
  
  dat_val <- list(x = beta_val,
                  y = type_val)
  
  # perform k-fold cross validation
  NFOLD <- 10
  foldid <- rep(1:NFOLD,ceiling(ncol(dat_tr$x)/NFOLD))[1:ncol(dat_tr$x)]
  
  fit.cv <- cv.glmnet(x = t(dat_tr$x), 
                      y = as.factor(as.vector(dat_tr$y)),
                      lambda = lambda,
                      type.measure="auc",
                      alpha=alpha,
                      family="binomial",
                      nfold=NFOLD,
                      foldid=foldid)
  
  
  # compute training auc
  predictor <- predict(fit.cv,
                       newx=t(dat_tr$x),
                       s="lambda.min")
  
  tr_predictor <- predict(fit.cv,
                          newx=t(dat_tr$x),
                          s="lambda.min",
                          type='response')
  
  tr_roc <- roc(dat_tr$y, as.numeric(predictor),quiet = TRUE)
  
  # compute validation auc
  predictor <- predict(fit.cv,
                       newx=t(dat_val$x),
                       s="lambda.min")
  
  val_predictor <- predict(fit.cv,
                           newx=t(dat_val$x),
                           s="lambda.min",
                           type='response')
  
  val_roc <- roc(dat_val$y, as.numeric(predictor), quiet = TRUE)
  
  # calibration slope
  x <- if(slope == "default"){
    ifelse(type_val=="Control",1,0)
  } else {
    ifelse(type_val=="Control",0,1)
  }
  r <- val.prob(val_predictor[,1], x)
  
  # Warning message if lambda.min too close to minimum/maximum
  if (fit.cv$lambda.min < (min(lambda)+min(lambda)*0.1) | fit.cv$lambda.min > (max(lambda) - max(lambda)*0.1)) warning("Optimum lambda is close to limits of lambda (range may be too small)")
  
  # return results
  return(list(tr_roc=tr_roc,
              val_roc=val_roc,
              fit.cv=fit.cv,
              oob_auc=mean(fit.cv$cvm),
              tr_predictor=tr_predictor,
              val_predictor=val_predictor,
              r=r))
  
}

train_linear_classifier <- function(beta_tr,
                                    type_tr,
                                    beta_val,
                                    type_val,
                                    alpha,
                                    lambda = "(-4,4)",
                                    n_seq,
                                    slope = "default"){
  
  if(max(n_seq) > nrow(beta_tr)){
    stop('n larger than beta number of rows')
  }
  
  val_auc <- rep(NA, length(n_seq))  # validation AUC
  tr_auc <- rep(NA, length(n_seq))   # training AUC
  oob_auc <- rep(NA, length(n_seq))  # out-of-bag validation AUC estimates
  lambda <- lambda # Customisable lambda
  cal_intercept <- rep(NA, length(n_seq)) # calibration intercept
  cal_slope <- rep(NA, length(n_seq)) # calibration slope
  val_predictor <- rep(NA, length(n_seq))
  
  counter <- 1
  for (n in n_seq){
    cat('\nRunning n =', n,'\n')
    res <- el_classifier(beta_tr[1:n,], type_tr,
                         beta_val[1:n,], type_val,
                         alpha=alpha,
                         lambda = lambda,
                         slope = slope)
    val_auc[counter] <- as.numeric(res$val_roc$auc)
    tr_auc[counter] <- as.numeric(res$tr_roc$auc)
    oob_auc[counter] <- as.numeric(res$oob_auc)
    cal_intercept[counter] <- as.numeric(res$r["Intercept"])
    cal_slope[counter] <- as.numeric(res$r["Slope"])
    val_predictor[counter] <- as.numeric(res$val_predictor)
    counter <- counter + 1
    cat('\nAUC =', round(res$val_roc$auc, digits=3),'\nIntercept = ', round(res$r["Intercept"], 3), '\nSlope = ', round(res$r["Slope"],3), '\n')
  }
  
  return(list(n_seq=n_seq,
              tr_auc=tr_auc,
              val_auc=val_auc,
              val_predictor=val_predictor,
              oob_auc=oob_auc,
              cal_intercept=cal_intercept,
              cal_slope=cal_slope))
}

merge_epi_imm <- function(beta_epi, beta_imm){
  
  if(!identical(colnames(beta_epi),colnames(beta_imm))){
    stop('epithelial and immune subsets have different sample names')
  } 
  if(nrow(beta_epi) != nrow(beta_imm)){
    stop('epithelial and immune subsets have different numbers of rows')
  }  
  
  # no of CpGs
  n <- nrow(beta_epi)
  
  # all CpG names
  cg_names <- rep(NA, 2*n)
  cg_names[seq(1,2*n,by=2)] <- rownames(beta_epi)
  cg_names[seq(2,2*n,by=2)] <- rownames(beta_imm)
  
  # Take the unique CpG names
  ucg <- unique(cg_names)
  
  # beta_epi_imm will contain all unique CpGs
  beta_epi_imm <- data.frame(matrix(NA, nrow=length(ucg), ncol=ncol(beta_epi)))
  rownames(beta_epi_imm) <- ucg
  colnames(beta_epi_imm) <- colnames(beta_epi)
  
  ind <- match(ucg, rownames(beta_epi))
  beta_epi_imm[!is.na(ind),] <- beta_epi[ind[!is.na(ind)],]
  
  ind <- match(ucg, rownames(beta_imm))
  beta_epi_imm[!is.na(ind),] <- beta_imm[ind[!is.na(ind)],]
  
  return(beta_epi_imm)
  
}