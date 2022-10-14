plot_roc <- function(type, index, title = NULL,
                     direction = "<",
                     colour = "gray20",
                     style = "default",
                     size = 5){

  require(pROC)
    roc <- roc(type, index, direction = direction)
    ci <- ci.auc(roc, conf.level = 0.95)
    
    anno <- paste('AUC = ', round(roc$auc[[1]], 2),'\n(95% CI: ', round(ci[1], 2),'-',round(ci[3],2),')',sep='')
    
    
    if(style == "lancet"){
      anno <- gsub("[.]", "·", anno)
    }

    
  ggplot(mapping=aes(x=1-roc$specificities,
                     y=(roc$sensitivities))) +
    geom_path(size = 0.75,
              colour = colour) +
    geom_line(aes(x=c(0,1),y=c(0,1)),color='gray',
              size = 0.5) +
    xlab("1 - Specificity") +
    ylab("Sensitivity") +
    annotate("label", 
             x = 0.70,
             y = 0.2,
             label = anno,
             size = size,
             fill = "#3C548899",
             colour = "black") +
    ggtitle(title) +
    theme(panel.background = element_blank(),
          title = element_text(size = 7),
          aspect.ratio = 1)
  
}
