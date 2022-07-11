
####
# calculate delta-beta consider both early and late apoptosis as "cases"
# note: three epithelial cell lines, one lymphocyte-like (immune cell) cell line
####

library(dplyr)
library(tidyverse)

dir.create("./10-output") #create directory for output

# load beta and pheno for all samples
beta <- readRDS("./0-output/beta_mims.Rds")
pheno <- readRDS("../../0-data/pheno_mims.Rds")

# subset beta and arrange ~ pheno
beta <- beta[match(rownames(pheno),colnames(beta))]
identical(rownames(pheno), colnames(beta))

#create factor type
pheno <- pheno %>%
  mutate(type = case_when(
    BH3.mims_type == "control" ~ "Control",
    BH3.mims_type %in% c("early apoptosis", "late apoptosis") ~ "case"
  ))

# calculate delta-beta's
db_mean_early <- c() #init
db_mean_late <- c() #init
db_mean <- c() #init
for (i in 1:nrow(beta)){
  db_mean_early[i] <- mean(as.numeric(beta[i,which(pheno$BH3.mims_type=='early apoptosis')])) - mean(as.numeric(beta[i,which(pheno$BH3.mims_type=='control')]))
  db_mean_late[i] <- mean(as.numeric(beta[i,which(pheno$BH3.mims_type=='late apoptosis')])) - mean(as.numeric(beta[i,which(pheno$BH3.mims_type=='control')]))
  db_mean[i] <- mean(as.numeric(beta[i,which(pheno$BH3.mims_type!='control')])) - mean(as.numeric(beta[i,which(pheno$BH3.mims_type=='control')]))
}

db <- data.frame(db_mean_early=db_mean_early, db_mean_late=db_mean_late, db_mean=db_mean)
rownames(db) <- rownames(beta)

saveRDS(db, file="./10-output/delta-beta.Rds")