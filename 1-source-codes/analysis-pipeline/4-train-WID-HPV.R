
# Train WID-HPV and save parameters for choosing best model

library(glmnet)
library(pROC)
library(rms)
library(dplyr)
library(tidyverse)

source("../train_linear_index-10-fold.R")


# Load subsetted beta matrices
beta_tr_epi <- readRDS(file = './3-output/beta_tr_epi.Rds')
beta_val_epi <- readRDS(file = './3-output/beta_val_epi.Rds')
beta_tr_imm <- readRDS(file = './3-output/beta_tr_imm.Rds')
beta_val_imm <- readRDS(file = './3-output/beta_val_imm.Rds')
beta_tr_pv <- readRDS(file = './3-output/beta_tr_pv.Rds')
beta_val_pv <- readRDS(file = './3-output/beta_val_pv.Rds')

# load pheno and select HPVneg/pos samples
pheno <- readRDS("../../0-data/pheno_P1.Rds") #read-in pheno
pheno <- pheno %>% filter(type %in% c("Control_NEG", "Control_POS"))

type_tmp <- rep('Control', nrow(pheno))
type_tmp[pheno$HPV=='POS'] <- 'POS'
pheno$type <- factor(type_tmp)

ind <- match(colnames(beta_tr_epi), rownames(pheno))
type_tr <- as.factor(as.character(pheno$type[ind]))
ic_tr <- as.numeric(pheno$ic[ind])

ind <- match(colnames(beta_val_epi), rownames(pheno))
type_val <- as.factor(as.character(pheno$type[ind]))
ic_val <- as.numeric(pheno$ic[ind])

n.seq <-c(2,10, 20, 30, 40, 50, 60, 70, 80, 90, 100,500,seq(1000,10000,by=1000))

# Train and validate classifiers

dir.create("4-output") #make dir for output

# ------- Ridge -------- #
res_ridge_epi <- train_linear_classifier(beta_tr = beta_tr_epi,
                                         type_tr = type_tr,
                                         beta_val = beta_val_epi,
                                         type_val = type_val,
                                         alpha = 0.0,
                                         n_seq = n.seq,
                                         lambda = "unrestricted",
                                         slope = "inv")

saveRDS(res_ridge_epi, file='./4-output/res_ridge_epi.Rds')



# ------- Lasso -------- #
res_lasso_epi <- train_linear_classifier(beta_tr = beta_tr_epi,
                                         type_tr = type_tr,
                                         beta_val = beta_val_epi,
                                         type_val = type_val,
                                         alpha = 1.0,
                                         n_seq = n.seq,
                                         lambda = "unrestricted",
                                         slope = "inv")

saveRDS(res_lasso_epi, file='./4-output/res_lasso_epi.Rds')


# ------- Ridge -------- #
res_ridge_imm <- train_linear_classifier(beta_tr = beta_tr_imm,
                                         type_tr = type_tr,
                                         beta_val = beta_val_imm,
                                         type_val = type_val,
                                         alpha = 0.0,
                                         n_seq = n.seq,
                                         lambda = "unrestricted",
                                         slope = "inv")

saveRDS(res_ridge_imm, file='./4-output/res_ridge_imm.Rds')



# ------- Lasso -------- #
res_lasso_imm <- train_linear_classifier(beta_tr = beta_tr_imm,
                                         type_tr = type_tr,
                                         beta_val = beta_val_imm,
                                         type_val = type_val,
                                         alpha = 1.0,
                                         n_seq = n.seq,
                                         lambda = "unrestricted",
                                         slope = "inv")

saveRDS(res_lasso_imm, file='./4-output/res_lasso_imm.Rds')


#==============================================================================#
# Merge epithelial and immune ranked subsets

beta_tr <- merge_epi_imm(beta_tr_epi, beta_tr_imm)
beta_val <- merge_epi_imm(beta_val_epi, beta_val_imm)

# ------- Ridge -------- #
res_ridge_epi_imm <- train_linear_classifier(beta_tr = beta_tr,
                                             type_tr = type_tr,
                                             beta_val = beta_val,
                                             type_val = type_val,
                                             alpha = 0.0,
                                             n_seq = n.seq,
                                             lambda = "unrestricted",
                                             slope = "inv")

saveRDS(res_ridge_epi_imm, file='./4-output/res_ridge_epi_imm.Rds')

# ------- Lasso -------- #
res_lasso_epi_imm <- train_linear_classifier(beta_tr = beta_tr,
                                             type_tr = type_tr,
                                             beta_val = beta_val,
                                             type_val = type_val,
                                             alpha = 1.0,
                                             n_seq = n.seq,
                                             lambda = "unrestricted",
                                             slope = "inv")

saveRDS(res_lasso_epi_imm, file='./4-output/res_lasso_epi_imm.Rds')

# only significant DMPs
db <- readRDS("./2-output/delta-beta.Rds")
db <- db[p.adjust(db[,3])<0.05,]

n <- as.numeric(nrow(beta_tr))

# ------- Ridge -------- #
res_ridge_sig <- train_linear_classifier(beta_tr = beta_tr,
                                             type_tr = type_tr,
                                             beta_val = beta_val,
                                             type_val = type_val,
                                             alpha = 0.0,
                                             n_seq = n,
                                             lambda = "unrestricted",
                                             slope = "inv")

saveRDS(res_ridge_sig, file='./4-output/res_ridge_sig.Rds')

# ------- Lasso -------- #
res_lasso_sig <- train_linear_classifier(beta_tr = beta_tr,
                                             type_tr = type_tr,
                                             beta_val = beta_val,
                                             type_val = type_val,
                                             alpha = 1.0,
                                             n_seq = n,
                                             lambda = "unrestricted",
                                             slope = "inv")
saveRDS(res_lasso_sig, file='./4-output/res_lasso_sig.Rds')

# F-statistic
ind <- match(colnames(beta_tr), rownames(pheno))
pheno_tr <- pheno[ind,]
identical(colnames(beta_tr), rownames(pheno_tr))

f <- numeric(nrow(beta_tr))
n <- 1:nrow(beta_tr)

f <- sapply(n, function(i){
  x <- lm(as.numeric(beta_tr[i,]) ~ as.factor(pheno_tr$type))
  f[i] <- summary(x)$f[1]
})

names(f) <- rownames(beta_tr)

f <- f[order(f, decreasing = TRUE)][1:30000]

ind <- match(names(f), rownames(beta_tr))
beta_tr <- beta_tr[ind,]

ind <- match(names(f), rownames(beta_val))
beta_val <- beta_val[ind,]

ind <- match(colnames(beta_val), rownames(pheno))
pheno_val <- pheno[ind,]
identical(colnames(beta_val), rownames(pheno_val))

pheno_tr$type <- droplevels(as.factor(pheno_tr$type))
pheno_val$type <- droplevels(as.factor(pheno_val$type))

n.seq <-c(2,10, 20, 30, 40, 50, 60, 70, 80, 90, 100,500,seq(1000,10000,by=1000))

# Train index
res_ridge <- train_linear_classifier(beta_tr = beta_tr,
                                     type_tr = pheno_tr$type,
                                     beta_val = beta_val,
                                     type_val = pheno_val$type,
                                     alpha = 0,
                                     n_seq = n.seq,
                                     lambda = "unrestricted",
                                     slope = "inv")

save(res_ridge, file = "./4-output/res_ridge_f.Rdata")

res_lasso <- train_linear_classifier(beta_tr = beta_tr,
                                     type_tr = pheno_tr$type,
                                     beta_val = beta_val,
                                     type_val = pheno_val$type,
                                     alpha = 1,
                                     n_seq = n.seq,
                                     lambda = "unrestricted",
                                     slope = "inv")

save(res_lasso, file = "./4-output/res_lasso_f.Rdata")

res_el <- train_linear_classifier(beta_tr = beta_tr,
                                  type_tr = pheno_tr$type,
                                  beta_val = beta_val,
                                  type_val = pheno_val$type,
                                  alpha = 0.5,
                                  n_seq = n.seq,
                                  lambda = "unrestricted",
                                  slope = "inv")

save(res_el, file = "./4-output/res_el.Rdata")



