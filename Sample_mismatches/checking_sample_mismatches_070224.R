# Code for checking for inconsistencies in data between Battenberg and DPClust files in PCAWG data
# 02.05.23/07.02.24
# Maria Jakobsdottir <maria.jakobsdottir@manchester.ac.uk>

library(tidyr)
library(dplyr)
# Samples
# files
metadata <- read.delim("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/release_may2016.v1.4.tsv", header = T, sep = "\t")
annotation <-  read.delim("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/evolution_and_heterogeneity/icgc_sample_annotations_summary_table.txt", header = T, sep = "\t")

# List of samples
samples <- subset(metadata, dcc_project_code %in% c("BRCA-UK","BRCA-EU","BRCA-US",
                                                    "OV-AU","OV-US",
                                                    "PRAD-UK","PRAD-CA","PRAD-US"))
samples$tumor_wgs_aliquot_id <- gsub(",.*","",samples$tumor_wgs_aliquot_id)
dim(samples)
samples <- merge(samples, annotation[,c("tumour_aliquot_id","wgd_status")], by.x = "tumor_wgs_aliquot_id", by.y = "tumour_aliquot_id")
dim(samples)

mismatched_samples <- as.data.frame(matrix(nrow = 0, ncol = 4))
colnames(mismatched_samples) <- c("sample_number","icgc_donor_id","tumor_wgs_aliquot_id","status")

for(i in 1:nrow(samples)){
  if(all(c(file.exists(paste0("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/Battenberg/rerun005_final_subclones/",
                              samples$tumor_wgs_aliquot_id[i],"_subclones.txt")),
           file.exists(paste0("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/DPClust/",
                              samples$tumor_wgs_aliquot_id[i],"/",
                              samples$tumor_wgs_aliquot_id[i],"_allDirichletProcessInfo.txt")),
           file.exists(paste0("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/final_consensus_12oct_passonly/snv_mnv/",
                              samples$tumor_wgs_aliquot_id[i],".consensus.20160830.somatic.snv_mnv.vcf.gz"))
  ))){
    
    # Read in files
    bb <- read.delim(paste0("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/Battenberg/rerun005_final_subclones/",
                            samples$tumor_wgs_aliquot_id[i],"_subclones.txt"),
                     sep = "\t", header = T)
    mult <- read.delim(paste0("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/DPClust/",
                              samples$tumor_wgs_aliquot_id[i],"/",
                              samples$tumor_wgs_aliquot_id[i],"_allDirichletProcessInfo.txt"))
    muts <- read.delim(paste0("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/final_consensus_12oct_passonly/snv_mnv/",
                              samples$tumor_wgs_aliquot_id[i],".consensus.20160830.somatic.snv_mnv.vcf.gz"),
                       header = FALSE, comment.char = "#")
    colnames(muts) <- c("CHROM",  "POS",     "ID",      "REF",     "ALT",     "QUAL",    "FILTER",  "INFO")
    
    
    muts_all <- muts[,c("CHROM","POS","POS","REF","ALT")]
    colnames(muts_all) <- c("chr","start","end","ref","alt")
    
    for(j in 1:nrow(bb)){
      
      tmp_bb <- bb[j,]
      
      tmp_mult <- subset(mult,
                         chr == tmp_bb[1,c("chr")] & 
                           start <= tmp_bb[1,c("endpos")] & 
                           end >= tmp_bb[1,c("startpos")])
      if((nrow(tmp_mult) > 0)  & !all(is.na(tmp_mult$nMaj1))){
        sub_tmp_mult <- tmp_mult[!is.na(tmp_mult$nMaj1),]
        
        if(!all(apply(sub_tmp_mult[,c("nMaj1","nMin1","frac1")], 1, function(x) x == tmp_bb[1,c("nMaj1_A","nMin1_A","frac1_A")]))){
          print(paste(j))
          
          sample_data <- as.data.frame(matrix(nrow = 0, ncol = 4))
          colnames(sample_data) <- c("sample_number","icgc_donor_id","tumor_wgs_aliquot_id","status")
          sample_data[1,1:4] <- c(i,samples$icgc_donor_id[i],samples$tumor_wgs_aliquot_id[i],"mismatch")
          mismatched_samples <- rbind(mismatched_samples,sample_data)
        }
      }
    }
  }
}

mismatched_samples <- distinct(mismatched_samples)
# 71 suspected mismatches
nrow(mismatched_samples)/nrow(samples) # 14.3% suspected mismatches?

write.table(mismatched_samples, "/mnt/bmh01-rds/UoOxford_David_W/b05055gj/suspected_PCAWG_BRCA_OV_PRAD_suspect_mismatches_unchecked_070224.txt", sep = "\t", row.names = FALSE, quote = FALSE)
