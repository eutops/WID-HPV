
### gene ontology testing for DMPs and index CpGs, background all CpGs on EPIC (that passed QC)

library(missMethyl)
library(tidyverse)
library(dplyr)
library(org.Hs.eg.db)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(GOxploreR)

dir.create("8-output") #make dir for output

######---------------- annotate DMPs -----------

## load in EPIC annotation
data("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
anno = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno <- as.data.frame(anno@listData) # holds all information

# select and save significant DMPs
db <- readRDS("./2-output/delta-beta.Rds")
db <- db[p.adjust(db[,3])<0.05,]
cpg_8386_annot <- anno %>% filter(Name %in% row.names(db)) #note that five probes are not annotated. These will be removed for the gene ontology testing
saveRDS("./8-output/cpg_8386_annot.Rds")
sig.cpg = as.character(c(cpg_8386_annot$Name)) 

# background CpGs
beta <- readRDS("./0-output/beta_P1.Rds") #read-in beta-matrix
all.cpg = row.names(beta)

# take subset beta_merged, only DMPs in HPV-ve, HPV+ve in TRAINING! set
beta_tr <- beta %>% filter(row.names(beta) %in% row.names(db)) 
pheno <- readRDS("../../0-data/pheno_P1.Rds") #read-in pheno
pheno_tr <- pheno %>% filter(split == "training set")
beta_tr <- beta_tr %>% dplyr::select(c(rownames(pheno_tr)))

identical(row.names(db),row.names(beta_tr))

#calculate mean db
db_mean <- c() #init
for (i in 1:nrow(beta_tr)){
  db_mean[i] <- mean(as.numeric(beta_tr[i,which(pheno_tr$type!='Control_NEG')])) - mean(as.numeric(beta_tr[i,which(pheno_trt$type=='Control_NEG')]))
}

db$db_mean <- db_mean # add mean delta-beta overall
saveRDS(db, file="./8-output/deltabeta_train_dmp.Rds")

######---------------- gene ontology testing all DMPs ----------------######

#### GO
DMPs_vs_all_result_GO <- gometh(
  sig.cpg,
  all.cpg = NULL,
  collection = c("GO"),
  array.type = c("EPIC"),
  plot.bias = FALSE,
  prior.prob = TRUE,
  anno = NULL,
  equiv.cpg = TRUE,
  fract.counts = TRUE,
  genomic.features = c("ALL"),
  sig.genes = TRUE # output significantly differentially methylated genes that overlap
)

saveRDS(DMPs_vs_all_result_GO, file="./8-output/DMPs_vs_all_result_GO.Rds")

## count number of DMPs per significant gene
result <- DMPs_vs_all_result_GO %>% filter(FDR<0.05) %>% droplevels #filter out only significant GO categories after FDRadj

# make vector with all significant genes
genes <- c() #init
for (i in 1:length(result$SigGenesInSet)){
  out <- unlist(strsplit(result$SigGenesInSet[i],","))
  genes <- c(genes,out)
}
genes <- unique(genes) #2812 genes


#annotate DMPs~ annotation file used by missMethyl
ann <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
flatAnn <- missMethyl:::.getFlatAnnotation(anno=ann)
sig.DMPs <- flatAnn %>% filter(cpg %in% row.names(db))

saveRDS(sig.DMPs, file="./8-output/DMPs_annot_missMethyl.Rds")

#get count
gene_count_DMP <- data.frame(unclass(table(sig.DMPs$alias)))
colnames(gene_count_DMP) <- c("DMP_count")
gene_count_DMP <- gene_count_DMP %>% filter(row.names(gene_count_DMP) %in% genes)   #filter only genes in GO enriched terms

saveRDS(gene_count_DMP, file="./8-output/DMP_count_siggenes_GO.Rds")


######---------------- gene ontology testing 27 WID-HPV index CpGs ----------------######

# load WID-HPV index CpGs
cpg_index_annot <- readRDS("./6-output/annotation_WID-HPV_index_CpGs.Rds")
sig.cpg = as.character(c(cpg_index_annot$Name)) 

#### GO index vs all
index_vs_all_result_GO <- gometh(
  sig.cpg,
  all.cpg = NULL,
  collection = c("GO"),
  array.type = c("EPIC"),
  plot.bias = FALSE,
  prior.prob = TRUE,
  anno = NULL,
  equiv.cpg = TRUE,
  fract.counts = TRUE,
  genomic.features = c("ALL"),
  sig.genes = TRUE # output significantly differentially methylated genes that overlap
)

#check <- index_vs_all_result_GO %>% filter(FDR<0.05) %>% droplevels #this yields nothing
check <- index_vs_all_result_GO %>% filter(P.DE<0.05) %>% droplevels #232 terms
saveRDS(index_vs_all_result_GO, file="./8-output/WIDHPV_vs_all_result_GO.Rds")

## Filter results for summary
table <- index_vs_all_result_GO %>% 
  filter(P.DE<0.05) %>% 
  rownames_to_column(var="ID") %>%
  droplevels %>% #232 terms
  dplyr::rename(Category=ONTOLOGY,Term=TERM,Genes=SigGenesInSet, pval=P.DE, adj_pval=FDR)
#table$GenehitsPerTerm <- sapply(table$Genes, function(x) length(unlist(strsplit(x, ",")))) #count nr of Genehits in each term, given already in DE

#Biological process only
table_BP <- table %>% filter(Category=="BP") %>% droplevels()

# only retain terms with more than 1 CpG hit
table_BP_filtered <- table_BP %>% filter(DE >= 2)
goterm <- table_BP_filtered %>% pull(ID)
#distRankingGO(goterm = goterm, domain = "BP", plot = TRUE)
#visRsubDAGBP(goterm = goterm, organism = "Human")

# for each gene combination, get highest and lowest level on DAG
combGene <- unique(table_BP_filtered$Genes) # all gene combinations (min.2 genes for BP)
levels <- GOTermBPOnLevel(goterm = goterm) %>% dplyr::rename(ID=Term) #get levels from DAG
table_BP_filtered <- full_join(levels,table_BP_filtered) #combine levels with table

#get for first element in combGene min/max df
df <- table_BP_filtered %>% filter(Genes == combGene[1])
dfMin <- df[which.min(df$Level),] #lowest level
dfMax <- df[which.max(df$Level),] #lowest level
dfMinMax <- bind_rows(dfMin,dfMax)

for (i in 2:length(combGene)){ #loop over remaining ones in list
  df <- table_BP_filtered %>% filter(Genes == combGene[i])
  dfMin <- df[which.min(df$Level),] #lowest level
  dfMax <- df[which.max(df$Level),] #lowest level
  dfMinMax <- bind_rows(dfMinMax,dfMin,dfMax)
}

# combine this info with prioritized version of the BP_filtered table
# removing redundant intermediate terms, retain only lowest level ends
# note: high-level terms are not retained
priority <- prioritizedGOTerms(lst = goterm, organism = "human", sp = FALSE, domain = "BP") 
priorityOut <- table_BP_filtered %>% filter(ID %in% c(priority$HF)) # 8 terms left
result <- bind_rows(dfMinMax, priorityOut) %>% distinct()

# simplify
result <- result %>%
  dplyr::select(ID,Level,Term,N,pval,Genes) %>%
  dplyr::arrange(Genes,Level)
saveRDS(result, "./8-output/summaryGOindex.Rds")


