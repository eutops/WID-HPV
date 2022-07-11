# Description: pipeline to load and process raw data from the Illumina EPIC array
# PART1

####---------------------------------------------------------------------
# For each Beadchip (plate) separately: 
# 1) Load idat files & extract p-values 
# 2) Filter low-quality samples
# 3) Do background intensity correction and a dye bias correction
# 4) Do BMIQ normalization
# 5) Save beta density as plotly plot -> these should be manually inspected, 
# although rare, some samples may have to be removed still
####---------------------------------------------------------------------

####### Section to be completed for each experiment
path_to_idat <- "<YOUR-PATH-1>/example_data/raw_data/" 
                   # specify path to folder containing all array data. 
                   # should contain sub-folders, holding for each BeadChip the individual idat files
                   # names of the sub-folders should be identical to the name of each BeadChip,
                   # listed in the vector "plates"

plates = c("205117390092",
           "205117390097",
           "205117390146",
           "205117390160") # enter names of BeadChips to be processed

path_to_output <- "<YOUR-PATH-2>/example_output/" #enter path to a separate folder where you want output stored

####### Section to be completed -end

library(minfi)
library(ChAMP)
library(plotly)
library(ewastools)

dir.create(file.path(path_to_output, "log"))
path_to_log <- file.path(path_to_output, "log")
dir.create(file.path(path_to_output, "plates_detP_beta"))
path_to_plates <- file.path(path_to_output, "plates_detP_beta")


# Define global thresholds
INTENSITY_THRESHOLD <- 9.5     # minimum median intensity required
DETECTION_P_THRESHOLD <- 0.01  # maximum detection p-value
FAILED_PROBE_THESHOLD <- 0.1   # maximum proportion of failed probes per sample


for (chip in plates){ #start loop over each Beadchip
  
# Initialise a list to log various parameters
log_data <- list(plate_name = chip, # name of plate
                 n_samples=NA,                # no of samples on plate
                 n_probes=NA,                 # no of probes in matrix
                 n_samples_removed=NA,        # no of samples that fail QC
                 snp_outlier_metric=NA,       # from ewastools::snp_outliers
                 rm_sample_list=NA)                     

##===========================================================##
# 1) Load data and extract p-values
RGset <- read.metharray.exp(base = file.path(path_to_idat, chip),
                            verbose = TRUE,
                            force = TRUE,
                            recursive = TRUE) # Load data and create RGset object
detP <- minfi::detectionP(RGset, type = "m+u") # Extract detection p-values from RGset

##===========================================================##
# 2) Filter low quality samples

Mset <- preprocessRaw(RGset) #Qc median intensity extraction
qc <- getQC(Mset) # Extract quality control information

beta_snp <- getSnpBeta(RGset)
genotypes <- call_genotypes(beta_snp, learn = FALSE, maxiter = 50) # Extract SNP outlier metric using ewastools
log_data$snp_outlier_metric <- snp_outliers(genotypes)

# Filter any samples with median (un)methylated intensity less than threshold
low_intensity_samples <- rownames(qc)[qc$mMed<INTENSITY_THRESHOLD | qc$uMed<INTENSITY_THRESHOLD]

# Filter samples with too many failed probes
failed_samples <- colnames(detP)[colSums(detP>DETECTION_P_THRESHOLD) > (nrow(detP) * FAILED_PROBE_THESHOLD)]
samples_to_remove <- unique(c(low_intensity_samples, failed_samples)) # can add bad_sample_list here

log_data$n_samples_removed <- length(samples_to_remove)

rm_ind <- match(samples_to_remove, colnames(RGset))
if(length(rm_ind)>0){
  RGset_filtered <- RGset[,-rm_ind]
  log_data$rm_sample_list <- samples_to_remove #log samples removed
} else {
  RGset_filtered <- RGset
}

##===========================================================##
# 3) background intensity correction and a dye bias correction

ssNOOB_filtered <- preprocessNoob(RGset_filtered, dyeCorr=TRUE, verbose=TRUE, dyeMethod='single') # Background intensity correction and dye bias correction
beta_ssNOOB_filtered <- getBeta(ssNOOB_filtered) # Extract corresponding beta values
log_data$n_samples <- ncol(beta_ssNOOB_filtered)
log_data$n_probes <- nrow(beta_ssNOOB_filtered)

##===========================================================##
# 4) BMIQ normalization (probe bias correction)
beta_ssNOOB_filtered_norm <- champ.norm(beta=beta_ssNOOB_filtered,arraytype="EPIC",cores=4) #arraytype="EPIC" or "450K"

##===========================================================##
# 5) Save beta density as plotly plot
d <- density(beta_ssNOOB_filtered_norm[,1], na.rm=TRUE, bw=0.02,
             from = -0.05, to = 1.05)
p <- plot_ly(x=d$x,y=d$y,
             mode = 'line',
             type='scatter',
             text=colnames(beta_ssNOOB_filtered_norm)[1])

for(j in 2:ncol(beta_ssNOOB_filtered_norm)){
  d <- density(beta_ssNOOB_filtered_norm[,j], na.rm=TRUE, bw=0.02,
               from = -0.05, to = 1.05)
  p <- p %>% add_trace(x=d$x,y=d$y,
                       mode = 'line',
                       type='scatter',
                       text=colnames(beta_ssNOOB_filtered_norm)[j])
}

p <- p %>%   layout(title = chip,
                    yaxis = list(title = 'Beta'),
                    xaxis = list(title = 'Density'),
                    showlegend=FALSE)%>%
  config(displayModeBar = FALSE)


##===========================================================##
# Save outputs

save(qc, file=paste(path_to_log,'/',
                    chip,'_qc.Rdata',sep=''))

save(beta_ssNOOB_filtered_norm, file=paste(path_to_plates,'/',
                                           chip,'_beta_ssNOOB_filtered_norm.Rdata',sep=''))

save(detP, file=paste(path_to_plates,'/',
                      chip,'_detP.Rdata',sep=''))

save(p, file=paste(path_to_log,'/',
                   chip,'_beta_density_plot.Rdata',sep=''))

save(log_data, file=paste(path_to_log,'/',
                         chip,'_log_data.Rdata',sep=''))

} # end for-loop
