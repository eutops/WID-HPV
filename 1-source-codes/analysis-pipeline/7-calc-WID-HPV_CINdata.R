
# calculate the WID-HPV for the full cytology dataset and save along with original phenotypic information and estimated ic content

dir.create("7-output") #make dir for output

# core function
index_HPV <- function(beta){
  # Load coefficients linear index
  scale <- readRDS('./6-output/scale.Rds')
  index <- readRDS('./6-output/index_coef.Rds')
  
  # Compute index
  ind <- na.omit(match(names(index), rownames(beta)))
  b <- beta[ind,]
  
  ind <- na.omit(match(rownames(b), names(index)))
  w <- index[ind] 
  
  if(!identical(names(w), rownames(b))){
    stop('***names mismatch***')
  }
  
  index_HPV <- colSums(b*w, na.rm = TRUE)
  index_HPV <- (index_HPV - mean(scale)) / sd(scale)
  return(index_HPV)
}

# load data
pheno <- readRDS("../../0-data/pheno_P1.Rds") #read-in pheno
beta <- readRDS("./data/beta_P1.Rds") #read-in beta-matrix

# calculate WID-HPV and save final pheno
HPV_index <- index_HPV(beta_P1)
pheno$WID_HPV <- HPV_index[match(row.names(pheno), names(HPV_index))]
saveRDS(pheno,file="./7-output/pheno_P1_ic_WID_HPV.Rds")
