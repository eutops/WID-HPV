### use Champ to determine significant DMPS overall (group early and late apoptosis)

library(ChAMP)
library(dplyr)
library(tidyverse)

dir.create("./11-output") #create directory for output

####Prepare pheno and beta
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

phenotype <- c(pheno$type)

## DMP for control versus BH3mims
DMP <- champ.DMP(beta = beta,
                  pheno = phenotype,
                  compare.group = NULL,
                  adjPVal = 1, #extremely small sample size, all adjPVal are ~1 so for adjPVal=0.05, no results!
                  adjust.method = "BH", 
                  arraytype = "EPIC"
                 ) 
# retain the 6308 DMPs with a P.Value <0.05
DMP <- DMP$Control_to_case %>% filter(P.Value <0.05)
DMP <- list(Control_to_case = DMP)
saveRDS(DMP,file="./11-output/DMP.Rds")