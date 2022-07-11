
# Description: train linear elastic net classifers as a function of the number of input variables.
# 
# Author: James E. Barrett
# Contact: regmjeb@ucl.ac.uk
# Date: 21 October 2019
# Edit: 10 November 2020: customisable lambda fuction 
# Edit: 24 March 2020: adding in calibration intercept/slope for evaluation
# Edit Charlotte: 21 February 2022: customizable k for k-fold for cross-validation
# using family=binomial (logistic regression)

# Important: the oob_auc is the mean of the cvm across different n-inputs,
# with cvm the mean cross-validated error for different lambdas in each training round
# for small datasets, glmnet autmatically changes the type.measure = "auc" to "deviance".
# In this case, the mean(cvm) does not equal the oob_auc and will be greater than 1


train_linear_classifier <- function(beta_tr, # training set for n-fold cross-validation (features/variables= beta @ CpG)
                                    type_tr, # pheno training set (binomial)
                                    beta_val, # holdout test data set
                                    type_val, # pheno test set
                                    alpha, # learning rate parameter, amount of change in coeficients on each update
                                    lambda = "(-4,4)", # input customizable lambda function, power to which e should be raised to calculate lambda,
                                    # with lambda the regularization rate, multiplied with regularization term for tuning
                                    # increase/decrease complexity and risk for over-/under-fitting
                                    k = 10, # number of folds in the training set used for crossvalidation, default=10
                                    n_seq, # number of input variables to train on
                                    slope = "default"){
  
  if(max(n_seq) > nrow(beta_tr)){
    stop('n larger than beta number of rows')
  }
  if (k < 3) {
    stop('k should be greater or equal than 3 for k-fold cross validation')
  }
  
  
  val_auc <- rep(NA, length(n_seq))  # init vector for validation AUC
  tr_auc <- rep(NA, length(n_seq))   # init vector for training AUC
  oob_auc <- rep(NA, length(n_seq))  # init vector for out-of-bag validation AUC estimates
  lambda <- lambda # Customisable lambda
  cal_intercept <- rep(NA, length(n_seq)) # init vector for calibration intercepts
  cal_slope <- rep(NA, length(n_seq)) # calibration slope
  val_predictor <- rep(NA, length(n_seq))
  
  #store evaluation scores obtained for each number of input CpGs
  counter <- 1
  for (n in n_seq){
    cat('\nRunning n =', n,'\n')
    res <- el_classifier(beta_tr[1:n,], type_tr,
                         beta_val[1:n,], type_val,
                         alpha=alpha,
                         lambda = lambda,
                         slope = slope,
                         k=k)
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


el_classifier <- function(beta_tr, type_tr,
                          beta_val, type_val,
                          alpha = 0.0,
                          lambda = "(-4,4)",
                          slope = "default",
                          k = 10){
  
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
  
  # perform k-fold cross validation with glmnet
  NFOLD <- k
  foldid <- rep(1:NFOLD,ceiling(ncol(dat_tr$x)/NFOLD))[1:ncol(dat_tr$x)]
  
  # actual model parameters for the FULL dataset
  fit.cv <- cv.glmnet(x = t(dat_tr$x), 
                      y = as.factor(as.vector(dat_tr$y)),
                      lambda = lambda,
                      type.measure="auc", # loss used for cross-validation
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
  
  # # -------- ROC curve ------- #
  # y <- rep(0,length(dat_val$y))
  # y[dat_val$y=='Endometrial'] <- 1
  # 
  # x <- as.numeric(predictor)
  # 
  # y <- y[order(x, decreasing = TRUE)]
  # TPR=cumsum(y)/sum(y)
  # FPR=cumsum(!y)/sum(!y)
  # 
  # plot(c(0,FPR),c(0,TPR),type='l',
  #      xlim=c(0,1),ylim=c(0,1),
  #      xlab='False positives',
  #      ylab='True postivies',
  #      main='PC1 ROC curve')
  # points(c(-1,2),c(-1,2),lty='dashed',type='l')
  # grid()
  # 
  
  # print results
  #cat('Training AUC = ',tr_roc$auc,'\n','Validation AUC = ',val_roc$auc,'\n',sep='')
  
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
