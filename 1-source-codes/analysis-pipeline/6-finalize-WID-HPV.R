
# Finalize WID-HPV by running best performing model/subset features, and scaling to hold out test

library(glmnet)
library(pROC)
library(rms)
library(dplyr)
library(tidyverse)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

source("../train_linear_index-10-fold.R")

#------------------------------------
# 1. Prepare data
#------------------------------------

beta_tr_epi <- readRDS(file = './3-output/beta_tr_epi.Rds')
beta_val_epi <- readRDS(file = './3-output/beta_val_epi.Rds')

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

#------------------------------------
# 2. Train
#------------------------------------
n <- 500
beta_tr <- beta_tr_epi[1:n,]
beta_val <- beta_val_epi[1:n,]

dir.create("6-output") #make dir for output

# Lasso
res <- el_classifier(beta_tr = beta_tr,
                     type_tr = type_tr, 
                     beta_val = beta_val,
                     type_val = type_val,
                     alpha = 1.0,
                     lambda = "unrestricted",
                     slope = "inv")
save(res, file = "./6-output/res.Rdata")

#------------------------------------
# 3. Save coefs
#------------------------------------

index_coef <- coef(res$fit.cv, s = "lambda.min")
names <- rownames(index_coef)
index_coef <- as.numeric(index_coef)
names(index_coef) <- names
index_coef <- index_coef[index_coef!=0]
saveRDS(index_coef, file='./6-output/index_coef.Rds')

pheno_val <- pheno[match(colnames(beta_val), rownames(pheno)),]
pheno_val$res <- res$val_predictor
save(pheno_val, file = "./6-output/pheno_val.Rdata")

#------------------------------------
# 4. Scaling index in external validation sett
# & saving output
#------------------------------------

index <- index_coef
# Compute index
ind <- na.omit(match(names(index), rownames(beta_val)))
b <- beta_val[ind,]

ind <- na.omit(match(rownames(b), names(index)))
w <- index[ind] 

if(!identical(names(w), rownames(b))){
  stop('***names mismatch***')
}

scale <- colSums(b*w, na.rm = TRUE)

pheno_val$scale <- scale
saveRDS(scale, file='./6-output/scale.Rds')

#------------------------------------
# 5. annotate index CpGs
#------------------------------------

## load in EPIC annotation
data("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
anno = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno <- as.data.frame(anno@listData) # holds all information

# list is nicely sorted, only get info you need
df <- data.frame(Name=anno$Name, 
                 Chromosome=anno$chr,
                 Gene=anno$UCSC_RefGene_Name,
                 Infinium.Probe.Type=anno$Type, 
                 CpG.context=anno$Relation_to_Island,
                 Gene.context=anno$UCSC_RefGene_Group)

#keep only first element
df$Gene <- gsub("\\;.*","",df$Gene)
df$Gene.context <- gsub("\\;.*","",df$Gene.context)

# assuming no annotation means non-coding
df$Gene.context[df$Gene.context==""] <-"IGR" 

# Keep only CpGs index
index_CpGs <- index_coef[-1]
index_CpGs <- sort(index_CpGs, decreasing=TRUE) #sort according to weight index
index_CpGs <- names(index_CpGs)
df <- df %>% filter(Name %in% index_CpGs) %>% droplevels()

## Add info on hyper/hypomethylation (epithelial component in training set)
db <- readRDS("./2-output/delta-beta.Rds") # load in delta-beta's
db <- as.data.frame(db)
db <- db %>% 
  filter(row.names(db) %in% index_CpGs)%>% # retain only index-Cpgs
  mutate(methylation.epithelial = case_when(db_epithelial <0 ~ "hypomethylated",db_epithelial >0 ~ "hypermethylated")) %>%
  droplevels
db <- db %>% select(db_epithelial,methylation.epithelial) %>% rownames_to_column(var="Name")

## combine annotation and delta-beta info
df <- full_join(df, db)
df <- df[match(index_CpGs, df$Name),] #sort Cpgs according to weight index

## save
saveRDS(df, file="./6-output/annotation_WID-HPV_index_CpGs.Rds")











