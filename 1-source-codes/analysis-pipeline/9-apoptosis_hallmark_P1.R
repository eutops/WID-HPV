
# Analysis methylation status promotor region apoptosis hall mark genes

library(ggplot2)
library(msigdbr)
library(dplyr)
library(stringr)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

dir.create("9-output") #make dir for output

# Get apoptosis genes
apoptosis <- msigdbr(species = "Homo sapiens",
                     category = "H")
apoptosis <- apoptosis %>%
  dplyr::filter(gs_id == "M5902") # apoptosis hallmark genes

apoptosis$human_gene_symbol # sanity checking the gene names

# Find apoptosis CpGs in TSS200/1500 and island on array
data <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
data <- data %>%
  as.data.frame() %>%
  mutate(genes = str_split(UCSC_RefGene_Name, ";", simplify = TRUE)[,1]) %>%
  mutate(reg = str_split(UCSC_RefGene_Group, ";", simplify = TRUE)[,1]) %>%
  filter(genes %in% apoptosis$human_gene_symbol & reg %in% c("TSS200", "TSS1500") & Relation_to_Island=="Island")
length(unique(data$genes)) # 103 out of 161 genes with promoter CpGs found

# Load beta for CIN data
beta <- readRDS("./0-output/beta_P1.Rds") #read-in beta-matrix
intersect <- intersect(rownames(data), rownames(beta))
beta <- beta[match(intersect, rownames(beta)),]
data <- data[match(intersect, rownames(data)),]
identical(rownames(beta), rownames(data))

# arrange pheno
pheno <- readRDS("../../0-data/pheno_P1.Rds") #read-in pheno
pheno <- pheno %>%
  dplyr::slice(match(colnames(beta), rownames(pheno)))
identical(rownames(pheno), colnames(beta))

# Add apoptosis genes
pheno$sum_apoptosis <- apply(beta, 2, sum)
pheno$mean_apoptosis <- apply(beta, 2, mean)
pheno$median_apoptosis <- apply(beta, 2, median)

# p53 methylation
data <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
data <- data %>%
  as.data.frame() %>%
  mutate(genes = str_split(UCSC_RefGene_Name, ";", simplify = TRUE)[,1]) %>%
  mutate(reg = str_split(UCSC_RefGene_Group, ";", simplify = TRUE)[,1]) %>%
  filter(genes == "TP53" & reg %in% c("TSS200", "TSS1500") & Relation_to_Island=="Island")

pheno$mean_p53 <- apply(beta, 2, mean)
pheno$median_p53 <- apply(beta, 2, median)
pheno$sum_p53 <- apply(beta, 2, sum)

save(pheno, file = "./9-output/dat_apoptosis_hallmark.Rdata")