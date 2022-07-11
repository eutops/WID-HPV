
# Finalize apoptosis index by running best subset features, and scaling to hold out test
# Model: Elastic net (3-fold CV, loss=deviance) after taking "f-subset" - n = 3000

library(dplyr)
library(tidyverse)
library(ggplot2)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

dir.create("./14-output") #create directory for output
source("../train_linear_index-k-fold.R") #functions for k-fold cross validation

#------------------------------------
# 1. Prepare data
#------------------------------------

# load pheno and beta
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

#select only signifcant CpGs, combine with delta-beta and sort ~ abs(delta-beta)
load("./12-output/2-output/f.Rdata")
df_sig <- df %>% filter(pValue < 0.05) %>% droplevels #15359 significantCpGs left
cpgs <- rownames(df_sig)
df_sig <- df_sig %>% rownames_to_column(var="CpG")
deltabeta <- readRDS("./10-output/delta-beta.Rds")
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

#------------------------------------
# 2. Train
#------------------------------------

n <- 3000
beta_tr <- beta_tr[1:n,]
beta_val <- beta_val[1:n,]

res <- el_classifier(beta_tr = beta_tr,
                     type_tr = pheno_tr$type, 
                     beta_val = beta_val,
                     type_val = pheno_val$type,
                     alpha = 0.5,
                     k=3,
                     lambda = "unrestricted")
save(res, file = "./14-output/res.Rdata")

#------------------------------------
# 3. Save coefs
#------------------------------------

#take coeficients for minimal lambda(most simple model, minimizing risk for overfit)
index_coef <- coef(res$fit.cv, s = "lambda.min")
names <- rownames(index_coef)
index_coef <- as.numeric(index_coef)
names(index_coef) <- names
index_coef <- index_coef[index_coef!=0]
saveRDS(index_coef, file='./14-output/index_coef.Rds')


#------------------------------------
# 4. Scale index according to the hold out test data set
# & saving output
#------------------------------------

#### code to scale in test dataset
pheno_val <- pheno[match(colnames(beta_val), rownames(pheno)),]
pheno_val$res <- res$val_predictor
save(pheno_val, file = "./14-output/pheno_val.Rdata")


index <- index_coef
## Compute index
ind <- na.omit(match(names(index), rownames(beta_val)))
b <- beta_val[ind,]

ind <- na.omit(match(rownames(b), names(index)))
w <- index[ind] 

if(!identical(names(w), rownames(b))){
  stop('***names mismatch***')
}

scale <- colSums(b*w, na.rm = TRUE)

pheno_val$scale <- scale

pheno_val %>%
  ggplot(aes(x = type,
             y = scale)) +
  geom_boxplot()

saveRDS(scale, file='./14-output/scale.Rds')

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

## Add info on hyper/hypomethylation (mean overall in training set)
db <- readRDS("./10-output/delta-beta.Rds") # load in delta-beta's
db <- db %>% 
  filter(row.names(db) %in% index_CpGs)%>% # retain only index-Cpgs
  mutate(methylation.mean = case_when(db_mean <0 ~ "hypomethylated",db_mean >0 ~ "hypermethylated")) %>%
  droplevels
db <- db %>% select(methylation.mean) %>% rownames_to_column(var="Name")

## combine annotation and delta-beta info
df <- full_join(df, db)
df <- df[match(index_CpGs, df$Name),] #sort Cpgs according to weight index

## save
saveRDS(df, file="./14-output/annotation_apoptosis_index_CpGs.Rds")










