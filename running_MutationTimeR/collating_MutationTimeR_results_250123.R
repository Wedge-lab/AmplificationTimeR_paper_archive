# script to pick out timed MYC segments from MutationTimeR output
# For AmplificationTimeR paper
# Run on results where the entire sample was run, rather than just MYC segment
# Have to run code in this way because it has to estimate ploidy and WGD etc

# 25.01.23

library(tidyr)

mut_time_files <- list.files("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/MutationTimeR_BRCA_OV_PRAD/", pattern = "12.txt")
mut_time_files <- as.data.frame(mut_time_files)
mut_time_files$filename <- mut_time_files$mut_time_files

mut_time_files <- separate(mut_time_files, col = "mut_time_files",
                           into = c("do_ID","wgs_aliquot_id","tool","muts"),
                           sep = "_")
mut_time_files <- mut_time_files[-1449,]


output_table <- as.data.frame(matrix(ncol = 21, nrow = 0))

colnames(output_table) <- c("seqnames","start","end","width","strand","major_cn","minor_cn",
                            "clonal_frequency","type","time","time.lo",
                            "time.up","time.2nd","time.2nd.lo","time.2nd.up","time.star","n.snv_mnv",
                            "do_ID","wgs_aliquot_id","tool","muts")

# ROI
amp_chr <- 8 # MYC chromosome
amp_start <- 128748315 # MYC start hg19
amp_stop <-  128753680 # MYC stop hg19

for(i in 1:nrow(mut_time_files)){
  mt_data <- read.delim(paste0("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/MutationTimeR_BRCA_OV_PRAD/",
                               mut_time_files$filename[i]), header = TRUE, sep = "\t")
  mt_data$do_ID <- mut_time_files$do_ID[i]
  mt_data$wgs_aliquot_id <- mut_time_files$wgs_aliquot_id[i]
  mt_data$tool <- mut_time_files$tool[i]
  mt_data$muts <- mut_time_files$muts[i]
  
  mt_data <- subset(mt_data, seqnames == amp_chr &
                      start <= amp_stop &
                      end >= amp_start)
  
  # select highest copy number subclone if there is more than 1
  if(nrow(mt_data) == 0){
    
  }else if(nrow(mt_data) == 1){
    output_table <- rbind(output_table, mt_data[,c("seqnames","start","end","width","strand","major_cn","minor_cn",
                                                 "clonal_frequency","type","time","time.lo",
                                                 "time.up","time.2nd","time.2nd.lo","time.2nd.up","time.star","n.snv_mnv",
                                                 "do_ID","wgs_aliquot_id","tool","muts")])
  }else if(nrow(mt_data) == 2){
    mt_data$sum <- mt_data$major_cn + mt_data$minor_cn
    mt_data <- mt_data[which.max(mt_data$sum),]
    
    output_table <- rbind(output_table, mt_data[,c("seqnames","start","end","width","strand","major_cn","minor_cn",
                                                   "clonal_frequency","type","time","time.lo",
                                                   "time.up","time.2nd","time.2nd.lo","time.2nd.up","time.star","n.snv_mnv",
                                                   "do_ID","wgs_aliquot_id","tool","muts")])
  }
  
  
  
}

write.table(output_table, "/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/MutationTimeR_BRCA_OV_PRAD/summary_MutationTimeR_MYC_2023-01-25.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
