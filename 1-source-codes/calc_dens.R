## calculate density distribution beta
calc_dens <- function(beta){
  
  beta <- as.data.frame(beta) # make sure beta is data frame format
  
  #calculate density distribution for first sample
  d <- density(beta[,1], na.rm=FALSE, bw=0.02,from = -0.05, to = 1.05)
  x <- d$x
  y <- d$y
  dens_distr <- data.frame(beta=x,density=y, EPIC_ID=colnames(beta)[1])
  
  #calculate density distribution for remaining samples and append
  for (i in 2:length(colnames(beta))){
    d <- density(beta[,i], na.rm=FALSE, bw=0.02,from = -0.05, to = 1.05)
    x <- d$x
    y <- d$y
    df <- data.frame(beta=x,density=y, EPIC_ID=colnames(beta)[i])
    dens_distr <- rbind(dens_distr,df)
  }
  
  return(dens_distr)
  
}