
# fit linear model to calculate delta-beta value in training subset of the discovery set (P1)

library(tidyverse)
library(dplyr)

dir.create("2-output") #make dir for output

beta <- readRDS("./0-output/beta_P1.Rds") # read in beta-matrix
pheno <- readRDS("./1-output/pheno_P1_ic.Rds") # read in phenotypic information

pheno <- pheno %>% filter(split=="training set") %>% droplevels #take training subset pheno
beta <- beta[,match(row.names(pheno),colnames(beta))] # take training subset beta
identical(colnames(beta), rownames(pheno)) #sanity-check

# Set factor levels and data types correctly
type_tmp <- rep('Control', nrow(pheno))
type_tmp[pheno$HPV=='POS'] <- 'POS'

pheno$type <- factor(type_tmp)
pheno$ic <- as.numeric(pheno$ic)
pheno$age <- as.numeric(pheno$age)

pheno <- pheno %>% select(type,age,ic) # keep only relevant variables

# fit linear models
cat('Estimating delta-betas...\n')
db <- matrix(NA, nrow=nrow(beta),ncol=2+ncol(pheno))
rownames(db) <- rownames(beta)
colnames(db) <- c('db_epithelial','db_immune',colnames(pheno))

ic.control <- pheno$ic[which(pheno$type=='Control')]
ic.case <- pheno$ic[which(pheno$type!='Control')]

pB <- txtProgressBar(min=1,max=nrow(beta), width =50L, style = 3)
for (i in 1:nrow(beta)){
  setTxtProgressBar(pB, i)
  
  ldat <- data.frame(beta=as.numeric(beta[i,]),
                     pheno)
  
  lfit <- lm(beta ~ ., data=ldat)
  db[i,3:ncol(db)] <- summary(lfit)$coefficients[2:(ncol(pheno)+1),4]
  
  beta.control <- as.numeric(beta[i,which(pheno$type=='Control')])
  beta.case <- as.numeric(beta[i,which(pheno$type!='Control')])
  
  dat.control <- data.frame(beta=beta.control, ic=ic.control)
  dat.case <- data.frame(beta=beta.case, ic=ic.case)
  
  fit.control <- lm(beta ~ ic, data=dat.control)
  fit.case <- lm(beta ~ ic, data=dat.case)
  
  delta_beta_epithelial <- round(fit.case$coefficients[1]-fit.control$coefficients[1],digits=3)
  delta_beta_immune <- round(-fit.control$coefficients[1] - fit.control$coefficients[2]
                             +fit.case$coefficients[1] + fit.case$coefficients[2],digits=3)
  db[i,1] <- delta_beta_epithelial
  db[i,2] <- delta_beta_immune
  
}
close(pB)
cat('done\n\n')

saveRDS(db, file="./2-output/delta-beta.Rds")

