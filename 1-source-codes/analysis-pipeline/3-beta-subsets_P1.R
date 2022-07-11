
# Create subsets for training WID-HPV
# Note that training and hold out test sets were already defined for another study and kept the same here

library(tidyverse)
library(dplyr)

dir.create("3-output") #make dir for output
path_to_output <- "./3-output"

# create training and hold out validation subsets and save in data folder
beta <- readRDS("./0-output/beta_P1.Rds") #read-in beta-matrix
pheno <- readRDS("../../0-data/pheno_P1.Rds") #read-in pheno

pheno_tr <- pheno %>% filter(split=="training set") %>% droplevels #take training subset pheno
beta_tr <- beta[,match(row.names(pheno_tr),colnames(beta))] # take training subset beta
identical(colnames(beta_tr), rownames(pheno_tr)) #sanity-check
saveRDS(file=beta_tr,"./3-output/beta_CO_tr.Rds") # save training set

pheno_val <- pheno %>% filter(split=="hold out test set") %>% droplevels
beta_val <- beta[,match(row.names(pheno_val),colnames(beta))]
identical(colnames(beta_val), rownames(pheno_val))
saveRDS(file=beta_val,"./3-output/beta_CO_val.Rds") # save hold out test set

# read in delta-beta matrix
db <- readRDS("./2-output/delta-beta.Rds")

# save top 30000 CpGs epithelial component
ord <- order(abs(db[,1]),decreasing = TRUE)[1:30000]
beta_tr_epi <- beta_tr[ord,]
beta_val_epi <- beta_val[ord,]
saveRDS(beta_tr_epi, file=paste0(path_to_output,"/beta_tr_epi.Rds"))
saveRDS(beta_val_epi, file=paste0(path_to_output,"/beta_val_epi.Rds"))

# save top 30000 CpGs immune cell component
ord <- order(abs(db[,2]),decreasing = TRUE)[1:30000]
beta_tr_imm <- beta_tr[ord,]
beta_val_imm <- beta_val[ord,]
saveRDS(beta_tr_imm, file=paste0(path_to_output,"/beta_tr_imm.Rds"))
saveRDS(beta_val_imm, file=paste0(path_to_output,"/beta_val_imm.Rds"))

# save top 30000 CpGs with lowest p-value linear fit
ord <- order(db[,3],decreasing = FALSE)[1:30000]
beta_tr_pv <- beta_tr[ord,]
beta_val_pv <- beta_val[ord,]
saveRDS(beta_tr_pv, file=paste0(path_to_output,"/beta_tr_pv.Rds"))
saveRDS(beta_val_pv, file=paste0(path_to_output,"/beta_val_pv.Rds"))