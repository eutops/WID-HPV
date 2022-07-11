#calculate apoptosis index for BH3-mims, external OC set and CIN datasets

library(dplyr)
library(tidyverse)
library(GEOquery)

dir.create("./15-output") #create directory for output

# core function
index_apop <- function(beta){
  # Load coefficients
  scale <- readRDS('./14-output/scale.Rds')
  index <- readRDS('./14-output/index_coef.Rds')
  
  # Compute index
  ind <- na.omit(match(names(index), rownames(beta)))
  b <- beta[ind,]
  
  ind <- na.omit(match(rownames(b), names(index)))
  w <- index[ind] 
  
  if(!identical(names(w), rownames(b))){
    stop('***names mismatch***')
  }
  
  apop_index <- colSums(b*w, na.rm = TRUE)
  apop_index <- ((apop_index - mean(scale)) / sd(scale))*(-1) # reverse so that control=negative,case=positive
  return(apop_index)
}

####----------- BH3-mims data -----------####

# load pheno and beta, arrange
beta <- readRDS("./0-output/beta_mims.Rds")
pheno <- readRDS("../../0-data/pheno_mims.Rds")
beta <- beta[match(rownames(pheno),colnames(beta))]
identical(rownames(pheno), colnames(beta))

#Calculate apoptosis index
apop_index <- index_apop(beta)

##combine with info pheno
pheno$apop_index <- apop_index[match(rownames(pheno), names(apop_index))]

#save 
saveRDS(pheno, file="./15-output/pheno_mims_apo_index.Rds")


####----------- platinum+azacitine treatment in HGSOC cell lines (external validation apoptosis index) -----------####
## HGSOC cells treated with DNMTi and carboplatin (n=8): 
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168225

gse <- getGEO("GSE168225")
pheno <- pData(gse$GSE168225_series_matrix.txt.gz)
beta <- exprs(gse$GSE168225_series_matrix.txt.gz)

#sum(is.na(beta))

### because there were NA's in beta-matrix, download and preprocess raw data instead
getGEOSuppFiles("GSE168225")
untar("GSE168225/GSE168225_RAW.tar", exdir = "GSE168225/idat")
head(list.files("GSE168225/idat", pattern = "idat"))
idatFiles <- list.files("GSE168225/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)

# standard preprocessing done, load final clean beta matrix
beta <- readRDS(file="../../0-data/beta_HGSOC.Rds")

#rename beta-colnames so that fits with pheno
colnames(beta) = gsub(pattern = "_.*", replacement = "", x = colnames(beta))

#plot beta-distribution
source("../calc_dens.R")
dens <- calc_dens(beta)
dens %>%
  ggplot(aes(x = beta,
             y = density, color=EPIC_ID)) + 
  geom_line()+
  xlab("Beta") +
  ylab("Density") +
  theme(legend.title=element_blank(), legend.text = element_text(size=8)) +
  guides(color=guide_legend(ncol=2)) 

# two samples show a very bad beta-distribution, these will be removed:
beta <- beta %>% select(!c("GSM5134231","GSM5134237"))
pheno <- pheno %>% filter(!rownames(pheno) %in% c("GSM5134231","GSM5134237"))
identical(rownames(pheno),colnames(beta))

#label pheno
pheno <- pheno %>%
  mutate(type=case_when(
    `treatment:ch1` == "Untreated" ~ "control",
    `treatment:ch1` == "Treated with 5microM azacitidine and 20microM carboplatin" ~ "azacitidine+carboplatin"
  ))

#Calculate apoptosis index
apop_index <- index_apop(beta)

##combine with info pheno
pheno$apop_index <- apop_index[match(rownames(pheno), names(apop_index))]

#save 
saveRDS(pheno, file="./15-output/pheno_GSE168225.Rds")



####----------- CIN data -----------####

##Calculate apoptosis index - P1
# load pheno and beta, arrange
beta <- readRDS("./0-output/beta_P1.Rds")
pheno <- readRDS("./7-output/pheno_P1_ic_WID_HPV.Rds")
beta <- beta[match(rownames(pheno),colnames(beta))]
identical(rownames(pheno), colnames(beta))

#Calculate apoptosis index
apop_index <- index_apop(beta)

##combine with info pheno
pheno$apop_index <- apop_index[match(rownames(pheno), names(apop_index))]

#save 
saveRDS(pheno, file="./15-output/pheno_P1_ic_WID_HPV_apo_index.Rds")
