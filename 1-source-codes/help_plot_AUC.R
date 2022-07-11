## helper functions for plotting AUC for model selection

combine_sets <- function(resid, x, type = "default"){
  if (type == "default"){
    m <- matrix(nrow = length(resid$n_seq), ncol = 5)
    colnames(m) <- c("n", "AUC_val", "AUC_oob", "AUC_tr", "Type")
    
    m[,1] <- resid$n_seq
    m[,2] <- resid$val_auc
    m[,3] <- resid$oob_auc
    m[,4] <- resid$tr_auc
    m[,5] <- rep(x, length(resid$n_seq))
    
    return(m)
  } else {
    m <- matrix(nrow = length(resid$n_seq), ncol = 6)
    colnames(m) <- c("n", "AUC_val", "AUC_oob", "AUC_tr", "Type", "Slope")
    
    m[,1] <- resid$n_seq
    m[,2] <- resid$val_auc
    m[,3] <- resid$oob_auc
    m[,4] <- resid$tr_auc
    m[,5] <- rep(x, length(resid$n_seq))
    m[,6] <- resid$cal_slope
    return (m)
  }
}

plot_performance_compare <- function(data, col, sets = 6){
  colour.pal.d8 <- c("#EA7580","#F2949A","#F6B3A1","#D8C99E","#3AB9AC","#109DB7","#0C6EA5","#172869")
  colour.pal.c <- colorRampPalette(colour.pal.d8)
  
  data %>%
    ggplot() +
    geom_line(aes(x = n,
                  y = col,
                  colour = Type)) +
    ylim(c(0.2, 1)) +
    xlab("Number of input CpGs") +
    ylab("AUC") +
    theme_minimal() +
    scale_colour_manual(values = colour.pal.c(sets),
                        aesthetics = "colour")+
    theme(legend.position = "top",
          legend.title = element_blank())
}


plot_performance_compare_slope <- function(data, col, sets = 6){
  colour.pal.d8 <- c("#EA7580","#F2949A","#F6B3A1","#D8C99E","#3AB9AC","#109DB7","#0C6EA5","#172869")
  colour.pal.c <- colorRampPalette(colour.pal.d8)
  
  data %>%
    ggplot() +
    geom_line(aes(x = n,
                  y = col,
                  colour = Type)) +
    ylim(c(-1, 1)) +
    xlab("Number of input CpGs") +
    ylab("Slope") +
    theme_minimal() +
    scale_colour_manual(values = colour.pal.c(sets),
                        aesthetics = "colour")+
    theme(legend.position = "top",
          legend.title = element_blank())
}

consensus_roc <- function(type, prob){
  
  #thres.seq <- seq(0,1,by=0.00001)
  thres.seq <- sort(as.numeric(prob))
  
  TPR <- rep(NA, length(thres.seq))
  FPR <- rep(NA, length(thres.seq))
  
  counter <- 1
  for (thres in thres.seq){
    
    prop <- apply(prob>thres, MARGIN = 2, FUN = 'sum')/nrow(prob)
    
    TPR[counter] <- sum(prop[type!='Control']<0.5)/sum(type!='Control')
    FPR[counter] <- sum(prop[type=='Control']>=0.5)/sum(type=='Control')
    
    counter <- counter + 1
  }
  
  x <- 1-FPR[match(unique(FPR),FPR)]
  y <- TPR[match(unique(FPR),FPR)]
  
  auc <- sum(y[2:(length(y))]*diff(x))
  
  # # for comparison
  # plot(roc(type,as.numeric(prob)))
  
  return(list(TPR=TPR,
              FPR=FPR,
              auc=auc))
}