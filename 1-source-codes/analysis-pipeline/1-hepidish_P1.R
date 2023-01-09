
#calculate IC content and IC subtypes for full cytology dataset (P1)

library(EpiDISH)
library(tidyverse)
library(dplyr)

beta <- readRDS("./0-output/beta_P1.Rds") #read-in beta-matrix
pheno <- readRDS("../../0-data/pheno_P1.Rds") #read-in pheno
pheno <- pheno[match(colnames(beta), row.names(pheno)),] #match order samples
identical(row.names(pheno), colnames(beta)) #sanity-check

dir.create("1-output") #make dir for output

#run epidish using reference for epithelial cells
out.l <- epidish(beta.m = beta, ref.m = centEpiFibIC.m, method = "RPC") 
save(out.l, file='./1-output/epidish_output.Rdata') #save output

pheno$ic <- out.l$estF[,3] # add estimated ic content to pheno
identical(rownames(out.l$estF), row.names(pheno)) #sanity-check
pheno <- saveRDS(pheno, file="./1-output/pheno_P1_ic.Rds") #save pheno

#run hepidish for deconvolution immune subset
frac.m <- hepidish(beta.m = beta, ref1.m = centEpiFibIC.m, ref2.m = centBloodSub.m, 
                   h.CT.idx = 3, method = 'RPC', maxit = 500) 
save(frac.m, file='./1-output/hepidish_output.Rdata') #save output


