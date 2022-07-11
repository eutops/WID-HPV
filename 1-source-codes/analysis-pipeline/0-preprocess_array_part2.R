# Description: pipeline to load and process raw data from the Illumina EPIC array
# PART2

####---------------------------------------------------------------------
# For each Experiment: 
# 1) Combine detection p-values from the different BeadChips (plates), determine which probes failed
# 2) Combine corrected and normalized beta-densities from the different BeadChips, after removing failed/flagged probes and failed samples
# FINAL OUTPUT: beta_merged.Rds

# To continue the pipeline after preprocessing, for the cin data save the cleaned beta matrix as follows:
# beta_P1 <- beta_merged
# saveRDS(beta_P1, file="./0-output/beta_P1.Rds")

# For the apoptsosis data save as:
# beta_mims <- beta_merged
# saveRDS(beta_mims, file="./0-output/beta_mims.Rds")

####---------------------------------------------------------------------

####### Section to be completed for each experiment

path_to_output <- "<YOUR-PATH-2>/example_output/" 
#when restarting an R-session after running PART1 of the pipeline, 
#the preferred output folder should be specified again

# based on the beta-distributions, you may want to manually remove some additional samples.
# In that case, make a list (.csv file) and specify path to list below (including file name)
path_to_samples_to_remove <- NULL

path_to_flagged <- "<YOUR-PATH-3>/flagged_probes/" 
# enter path to flagged-probes folder
# example contains a list of 2932 non CpG probes, 82108 probes flagged by Zhou et al. (10.1093/nar/gkw967), 6102 SNP probes, 537 probes on Y-chromosome

####### Section to be completed -end

library(impute)

path_to_plates <- file.path(path_to_output, "plates_detP_beta")
path_to_log <- file.path(path_to_output, "log")

DETECTION_P_THRESHOLD <- 0.01   # maximum detection p-value
FAILED_SAMPLE_THESHOLD <- 0.1   # maximum proportion of failed samples per probe


#combine all flagged probes, to be removed 
setwd(path_to_flagged)
flagged = list.files(path = path_to_flagged, pattern="*.csv") #all files with flagged probes
list2env(
  lapply(setNames(flagged, make.names(gsub("*.csv$", "", flagged))), 
         read.csv, header=FALSE), envir = .GlobalEnv) #read in files with flagged probes
flagged <- gsub("*.csv$", "", flagged) #all objects with flagged probes

flagged_names <- vector()
for(cpgs in flagged){
  dat <- get(cpgs)
  dat <- as.vector(dat[,1])
  flagged_names <- c(dat, flagged_names)
} 

# List all samples (manually identified after PART1 of the pre-processing pipeline) to be removed 
if(!is.null(path_to_samples_to_remove)){
    dat <- read.csv(path_to_samples_to_remove,header = FALSE)
    bad_sample_list <- as.vector(dat[,1])
} else {
    bad_sample_list <- NULL
  }

##===========================================================##
# 1) Combine detection p-values from the different BeadChips (plates)

# Load and combine the beta and detP matrices
input_dir_list <- dir(path_to_plates)
plates <- NULL
for (chip in input_dir_list){
  plates <- c(plates,strsplit(chip,split='_')[[1]][1])
}
plates <- unique(plates) #list of Beadchips for which pre-processing pipeline PART1 was completed

load(paste(path_to_plates,'/',plates[1],'_detP.Rdata',sep='')) # load detP matrix for first chip
detP_merged <- detP
if(length(plates)>1){ 
  for(p in plates[2:length(plates)]){ #loop over files for remaining chips
    load(paste(path_to_plates,'/',p,'_detP.Rdata',sep='')) #load next detP matrix
    detP_merged <- merge(detP_merged, detP,
                         by.x='row.names', by.y='row.names', sort=FALSE) #merge next detP matrix
    row.names(detP_merged) <- detP_merged$Row.names # reassign the row names
    detP_merged$Row.names <- NULL # remove this column created by the merge function
  }
} # combine all detP matrices


failed_probes <- rownames(detP_merged)[rowSums(detP_merged>DETECTION_P_THRESHOLD)>(ncol(detP_merged)*FAILED_SAMPLE_THESHOLD)] # Determine which probes failed
rm_names <- unique(c(flagged_names,failed_probes)) # combine all probes to be removed

rm(detP, chip, p); gc() # clean-up global environment

##===========================================================##
# 2) Combine corrected and normalized beta-densities with detP matrixes from the different BeadChips

na_count <- 0
beta_merged <- NULL

for(p in plates){
  
  load(paste(path_to_plates,'/',p,'_beta_ssNOOB_filtered_norm.Rdata',sep=''))
  load(paste(path_to_plates,'/',p,'_detP.Rdata',sep=''))
  
  # remove probes and any bad samples   
  ind.row <- match(rownames(beta_ssNOOB_filtered_norm),rm_names)
  ind.col <- match(colnames(beta_ssNOOB_filtered_norm),bad_sample_list)
  
  beta_plate_filtered <- beta_ssNOOB_filtered_norm[is.na(ind.row), is.na(ind.col)]
  rm(beta_ssNOOB_filtered_norm);invisible(gc())
  
  # match detection p-values to beta matrix
  ind.row <- match(rownames(beta_plate_filtered),rownames(detP))
  ind.col <- match(colnames(beta_plate_filtered),colnames(detP))
  detP_plate <- detP[ind.row, ind.col]
  
  # replace failed data points with NA
  beta_plate_filtered[detP_plate > DETECTION_P_THRESHOLD] <- NA
  na_count <- na_count + sum(is.na(beta_plate_filtered))
  rm(detP);rm(detP_plate);invisible(gc())
  
  # Imputation
  r <- rowSums(is.na(beta_plate_filtered))
  if(sum(r > 0.8*ncol(beta_plate_filtered))>0){
    beta_plate_filtered <- beta_plate_filtered[(r < 0.8*ncol(beta_plate_filtered))>0, ]
    cat('!! Probe failure in more than 80% of samples for',
        sum(r > 0.8*ncol(beta_plate_filtered)),
        'probes on plate',p,'(removed from beta matrix)\n')
  }
  
  out <- capture.output(beta_plate_filtered_imputed <- impute.knn(beta_plate_filtered,k=10,rowmax=0.8)$data)
  rm(beta_plate_filtered);invisible(gc()) # clean-up global environment
  
  
  if(is.null(beta_merged)){
    
    beta_merged <- beta_plate_filtered_imputed
    rm(beta_plate_filtered_imputed);invisible(gc()) # clean-up global environment
    
  } else {
    
    beta_merged <- merge(beta_merged, beta_plate_filtered_imputed,
                         by.x='row.names', by.y='row.names', sort=FALSE)
    rm(beta_plate_filtered_imputed);invisible(gc()) # clean-up global environment
    row.names(beta_merged) <- beta_merged$Row.names # reasign the row names
    beta_merged$Row.names <- NULL # remove this column created by the merge function
  }
  
  cat(paste('Merged plate ',p,'\n',sep=''))
}

cat(paste(na_count,
          ' (',
          round(na_count/(ncol(beta_merged) * nrow(beta_merged)),digits=4),
          '%) ',
          'data points failed detection p-value test\n',sep=''))


cat('\nMerged beta matrix has',
    nrow(beta_merged),
    'CpGs and',
    ncol(beta_merged),
    'samples\n\n')

info <- paste(na_count,
                      ' (',
                      round(na_count/(ncol(beta_merged) * nrow(beta_merged)),digits=4),
                      '%) ',
                      'data points failed detection p-value test. ','Merged beta matrix has ',
                      nrow(beta_merged),
                      'CpGs and ',
                      ncol(beta_merged),
                      ' samples.',sep='')

##===========================================================##
# Save output

save(failed_probes, file=paste(path_to_log,'/failed_probes.Rdata',sep='')) #log failed probes removed
write.table(info, file=file.path(path_to_log, "beta_merged_info.csv"), row.names = FALSE, col.names = FALSE)
save(beta_merged, file=paste(path_to_output,'/beta_merged.Rdata',sep=''))





