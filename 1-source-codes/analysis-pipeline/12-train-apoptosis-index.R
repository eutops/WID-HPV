#### use four different approaches to train an apoptosis index using 3-fold cross-validation

## Approach 1 top 10000 CpGs with largest absolute mean delta-beta
## Approach 2 top 10000 CpGs ~ Fstatsitic, sort absolute mean delta-beta
## Approach 3 input significant CpGs (p.adj) for type but not cell-type, sort absolute mean delta-beta
## Approach 4 input only significant DMPs (limna-Champ, p-value <0.01), sort absolute mean delta-beta

# split data for training: one split 1/4 test (hold out for final validation), 3/4 training, but ensure only one control sample in the test set

library(dplyr)
library(tidyverse)
library(fs)

dir.create("./12-output") #create directory for output
source("../train_linear_index-k-fold.R") #functions for k-fold cross validation
n.seq <-c(2,10, 20, 30, 40, 50, 60, 70, 80, 90, 100,500,seq(1000,10000,by=1000))

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

# define train-test split
pheno$type <- as.factor(pheno$type)
pheno_tr <- pheno %>% filter(!rownames(pheno) %in% c("206119360111_R02C01","206119360111_R06C01","206119360122_R04C01")) %>% droplevels
pheno_val <- pheno %>% filter(rownames(pheno) %in% c("206119360111_R02C01","206119360111_R06C01","206119360122_R04C01")) %>% droplevels

####----------- Approach 1 top 10000 CpGs with largest absolute mean delta-beta -----------####

setwd("./12-output")
dir.create("./1-output")

# load previous delta-beta calculations, complete beta-matrix and pheno
deltabeta <- readRDS("../10-output/delta-beta.Rds")

# take top 30,000 delta-beta CpGs, but only 10000 will be used for training the apoptosis index
db_sorted <- deltabeta %>% arrange(desc(abs(db_mean))) %>% slice(1:30000) 

# subset beta for training and validation
beta_tr <- beta %>% 
  filter(rownames(beta) %in% rownames(db_sorted)) %>%
  select(rownames(pheno_tr)) %>% droplevels
identical(rownames(pheno_tr), colnames(beta_tr))

beta_val <- beta %>% select(rownames(pheno_val)) %>% droplevels
identical(rownames(pheno_val), colnames(beta_val))

# Train classifier - Elastic net regularization with 3-fold cross-validation
res_el <- train_linear_classifier(beta_tr =beta_tr,
                                  type_tr = pheno_tr$type,
                                  beta_val = beta_val,
                                  type_val = pheno_val$type,
                                  alpha = 0.5,
                                  n_seq = n.seq,
                                  k=3,
                                  lambda = "unrestricted")

# due to the small amount of samples, glmnet set type.measure to "deviance"
# there oob_auc is not an auc but the mean of the cvm (mean cross-validated error = deviance for logistic regression)

saveRDS(res_el, file = "./1-output/res_el.Rds")

####----------- Approach 2 top 10000 CpGs ~ Fstatistic -----------####

dir.create("./2-output/")

# calculate F-statistic and corresponding p-values
f <- c() #init
p <- c() #init

for (i in 1:nrow(beta)){
  x <- lm(as.numeric(beta[i,]) ~ as.factor(pheno$type))
  f <- c(f,summary(x)$f[1])
  p <- c(p,summary(x)$coefficients[,"Pr(>|t|)"][[2]])
}
df <- data.frame(row.names=rownames(beta), Fstatistic = f, pValue = p)
save(df, file = "./2-output/f.Rdata")

# Subset significant Cpgs
df_sig <- df %>% filter(pValue < 0.05) %>% droplevels #15359 CpGs left
cpgs <- rownames(df_sig)
df_sig <- df_sig %>% rownames_to_column(var="CpG")

# sort full set of significant CpGs according to db_mean and subset beta for training and validation
db_subset <- deltabeta %>% 
  filter(rownames(deltabeta) %in% cpgs) %>% 
  select(db_mean) %>%
  rownames_to_column(var="CpG")
df_sig <- left_join(df_sig,db_subset) %>% column_to_rownames(var="CpG")
df_sig <- df_sig %>% arrange(desc(abs(db_mean)))

beta_tr <- beta %>% filter(rownames(beta) %in% c(rownames(df_sig))) %>%
  select(rownames(pheno_tr)) %>% droplevels
identical(rownames(pheno_tr), colnames(beta_tr))

beta_val <- beta %>% filter(rownames(beta) %in% c(rownames(df_sig))) %>%
  select(rownames(pheno_val)) %>% droplevels
identical(rownames(pheno_val), colnames(beta_val))

# Train classifier - Elastic net regularization with 3-fold cross-validation
res_el <- train_linear_classifier(beta_tr,
                                  type_tr = pheno_tr$type,
                                  beta_val = beta_val,
                                  type_val = pheno_val$type,
                                  alpha = 0.5,
                                  n_seq = n.seq,
                                  k=3,
                                  lambda = "unrestricted")

saveRDS(res_el, file = "./2-output/res_el.Rds")


####----------- Approach 3 input significant CpGs (p.adj) for type but not cell-type -----------####

dir.create("./3-output/")

# rearrange pheno
pheno <- pheno %>% mutate(line=case_when(
  Sample_type %in% c("Cal51","HT1080","HeLa") ~ "epithelial",
  Sample_type == "Nalm6" ~ "immune"
))
pheno$line <- as.factor(pheno$line)

# calculate relevant p-values and adjust
p_type <- c() #init
p_line <- c()

for (i in 1:nrow(beta)){
  x <- lm(as.numeric(beta[i,]) ~ pheno$type + pheno$line)
  p_type <- c(p_type,summary(x)$coefficients[,"Pr(>|t|)"][[2]])
  p_line <- c(p_line,summary(x)$coefficients[,"Pr(>|t|)"][[3]])
}

df <- data.frame(row.names=rownames(beta), pValue_type = p_type, pValue_line = p_line)
save(df, file = "./3-output/p.Rdata")

# Subset significant Cpgs
df_sig <- df %>% filter(pValue_type < 0.05 & pValue_line >= 0.05) %>% droplevels #13572 CpGs left
cpgs <- rownames(df_sig)
df_sig <- df_sig %>% rownames_to_column(var="CpG")

# sort full set of significant CpGs according to db_mean and subset beta for training and validation
db_subset <- deltabeta %>% 
  filter(rownames(deltabeta) %in% cpgs) %>% 
  select(db_mean) %>%
  rownames_to_column(var="CpG")
df_sig <- left_join(df_sig,db_subset) %>% column_to_rownames(var="CpG")
df_sig <- df_sig %>% arrange(desc(abs(db_mean)))

beta_tr <- beta %>% filter(rownames(beta) %in% c(rownames(df_sig))) %>%
  select(rownames(pheno_tr)) %>% droplevels
identical(rownames(pheno_tr), colnames(beta_tr))

beta_val <- beta %>% filter(rownames(beta) %in% c(rownames(df_sig))) %>%
  select(rownames(pheno_val)) %>% droplevels
identical(rownames(pheno_val), colnames(beta_val))

# Train classifier - Elastic net regularization with 3-fold cross-validation
res_el <- train_linear_classifier(beta_tr,
                                  type_tr = pheno_tr$type,
                                  beta_val = beta_val,
                                  type_val = pheno_val$type,
                                  alpha = 0.5,
                                  n_seq = n.seq,
                                  k=3,
                                  lambda = "unrestricted")

saveRDS(res_el, file = "./3-output/res_el.Rds")

####----------- Approach 4 input only significant DMPs (limna-Champ, p-value <0.01) -----------####

dir.create("./4-output/")

# subset significant DMPs and subset beta for training and validation
DMP <- readRDS("../11-output/DMP.Rds")
DMP <- DMP$Control_to_case
DMP <- DMP %>% arrange(desc(abs(deltaBeta)))

beta_tr <- beta %>% filter(rownames(beta) %in% c(row.names(DMP))) %>%
  select(rownames(pheno_tr)) %>% droplevels
identical(rownames(pheno_tr), colnames(beta_tr))

beta_val <- beta %>% filter(rownames(beta) %in% c(row.names(DMP))) %>%
  select(rownames(pheno_val)) %>% droplevels
identical(rownames(pheno_val), colnames(beta_val))

# Train classifier - Elastic net regularization with 3-fold cross-validation

n.seq <-c(2,10, 20, 30, 40, 50, 60, 70, 80, 90, 100,500,seq(1000,6308,by=1000))
res_el <- train_linear_classifier(beta_tr,
                                  type_tr = pheno_tr$type,
                                  beta_val = beta_val,
                                  type_val = pheno_val$type,
                                  alpha = 0.5,
                                  n_seq = n.seq,
                                  k=3,
                                  lambda = "unrestricted")

saveRDS(res_el, file = "./4-output/res_el.Rds")


