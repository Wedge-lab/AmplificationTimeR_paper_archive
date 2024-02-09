# code to make a list of all chromosome 8q segments identified in PCAWG data
# written by Maria
# Run in interactive session on 16.11.23
# using R 4.2.2 on csf3

library(tidyr)
# files
metadata <- read.delim("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/release_may2016.v1.4.tsv", header = T, sep = "\t")
annotation <-  read.delim("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/evolution_and_heterogeneity/icgc_sample_annotations_summary_table.txt", header = T, sep = "\t")

# List of samples
samples <- subset(metadata, dcc_project_code %in% c("BRCA-UK","BRCA-EU","BRCA-US",
                                                    "OV-AU","OV-US",
                                                    "PRAD-UK","PRAD-CA","PRAD-US"))
samples <- separate_rows(samples, tumor_wgs_aliquot_id, sep = ",")
dim(samples)
samples <- merge(samples, annotation[,c("tumour_aliquot_id","wgd_status")], by.x = "tumor_wgs_aliquot_id", by.y = "tumour_aliquot_id")
dim(samples)
scna.files <- paste0("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/Battenberg/rerun005_final_subclones/",
       samples$tumor_wgs_aliquot_id,"_subclones.txt")
scna.files <- as.data.frame(scna.files)
colnames(scna.files) <- c("scna.filename")
scna.files <- as.data.frame(scna.files[file.exists(scna.files$scna.filename),])
colnames(scna.files) <- c("scna.filename")

c8_segs <- as.data.frame(matrix(nrow = 0, ncol = 14))
colnames(c8_segs) <- c("samplename","chr","startpos","endpos","n1A","n2A","Tumour_tot_CN",
                       "nMaj1_A","nMin1_A","frac1_A",
                       "nMaj2_A","nMin2_A","frac2_A","length")

for(i in 1:nrow(scna.files)){
  
  sample.cn <-  read.table(scna.files$scna.filename[i], header = T, sep = "\t")
  sample.cn$length <- (sample.cn$end - sample.cn$start)+1
  sample.cn$samplename <- scna.files$scna.filename[i]
  
  sample.cn$n1A <- sample.cn$nMin1_A + sample.cn$nMaj1_A
  sample.cn$n2A <- sample.cn$nMin2_A + sample.cn$nMaj2_A
  
  sample.cn$n1_tot_cn <- (sample.cn$n1A*sample.cn$frac1_A)
  sample.cn$n2_tot_cn <- (sample.cn$n2A*sample.cn$frac2_A)
  sample.cn$tot_cn <- rowSums(sample.cn[,c("n1_tot_cn","n2_tot_cn")], na.rm = TRUE)
  
  tum_tot_cn <- sum(sample.cn$tot_cn*sample.cn$length)/(sum(sample.cn$length))
  
  sample.cn$Tumour_tot_CN <- tum_tot_cn
  
  sample.c8 <- subset(sample.cn, chr == 8 &
                        startpos <= 146364021 & # end of chromosome from UCSC hg19 genome browser
                        endpos >= 46047109) # centromere from UCSC hg19 genome browser
  
  c8_segs <- rbind(c8_segs, sample.c8[,c("samplename","chr","startpos","endpos","n1A","n2A","Tumour_tot_CN",
                                         "nMaj1_A","nMin1_A","frac1_A",
                                         "nMaj2_A","nMin2_A","frac2_A","length")])
  
  print(paste(i,"of",nrow(scna.files)))
}

write.table(c8_segs, paste0("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/20231116_AmplificationTimeR_review_response/BRCA_OV_PRAD_all_8q_segments_",Sys.Date(),".txt"), sep = "\t", row.names = FALSE)
