# Code to calculate gain timing using AmplificationTimeR
# For AmplificationTimeR manuscript

# Maria Jakobsdottir
# 23.11.23

# Run on PCAWG BRCA-UK, OV-AU, and PRAD-UK
# Use all mutations, C>T at CpG, and SBS1 and SBS5

### Libraries
# devtools::install_github("Wedge-lab/AmplificationTimeR", force = TRUE)
library(AmplificationTimeR)
library(tidyr)
library(BSgenome.Hsapiens.UCSC.hg19)

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

# prepare output table
output_tab <- as.data.frame(matrix(nrow = nrow(samples), ncol = 48))
colnames(output_tab) <- c("sample","region","highest_copy_number","event_order","num_mutations_used","clonality_status","flags",
                          "t_1","t_1_mean_bootstrap","t_1_lower_ci","t_1_upper_ci",
                          "t_2","t_2_mean_bootstrap","t_2_lower_ci","t_2_upper_ci",
                          "t_3","t_3_mean_bootstrap","t_3_lower_ci","t_3_upper_ci",
                          "t_4","t_4_mean_bootstrap","t_4_lower_ci","t_4_upper_ci",
                          "t_5","t_5_mean_bootstrap","t_5_lower_ci","t_5_upper_ci",
                          "t_6","t_6_mean_bootstrap","t_6_lower_ci","t_6_upper_ci",
                          "t_7","t_7_mean_bootstrap","t_7_lower_ci","t_7_upper_ci",
                          "t_8","t_8_mean_bootstrap","t_8_lower_ci","t_8_upper_ci",
                          "t_9","t_9_mean_bootstrap","t_9_lower_ci","t_9_upper_ci",
                          "t_10","t_10_mean_bootstrap","t_10_lower_ci","t_10_upper_ci","tumor_wgs_aliquot_id")
output_tab$tumor_wgs_aliquot_id <- samples$tumor_wgs_aliquot_id

output_tab_all_muts <- output_tab
output_tab_ctcpg <- output_tab
output_tab_sbs15 <- output_tab

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
    
    if(samples$wgd_status[i] == "no_wgd"){
      wgd_status <- FALSE
    }else if(samples$wgd_status[i] == "wgd"){
      wgd_status <- TRUE
    }
    
    # Paramters of amplification being investigated
    amp_chr <- 8 # MYC chromosome
    amp_start <- 128748315 # MYC start hg19
    amp_stop <-  128753680 # MYC stop hg19
    amp_wgd <- wgd_status
    amp_name <- samples$sample[i]
    
    tryCatch({time_all_muts <- time_amplification(cn_data = bb,
                                                  multiplicity_data = mult,
                                                  mutation_data = muts_all,
                                                  muts_type = "SBS1 and SBS5",
                                                  sample_id = samples$icgc_donor_id[i],
                                                  amplification_chrom = amp_chr,
                                                  amplification_start = amp_start,
                                                  amplification_stop = amp_stop,
                                                  is_WGD = amp_wgd,
                                                  genome = "hg19")
    
    output_tab_all_muts[i,1:46] <- time_all_muts}, error=function(e){})
    
    if(file.exists(paste0("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/single_mutation_signature_attribution/",
                          samples$icgc_donor_id[i],"_",
                          samples$tumor_wgs_aliquot_id[i],"_individual_mutation_signature_attributions_2022-11-21.txt"))){
      
      muts_sbs <- read.delim(paste0("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/single_mutation_signature_attribution/",
                                    samples$icgc_donor_id[i],"_",
                                    samples$tumor_wgs_aliquot_id[i],"_individual_mutation_signature_attributions_2022-11-21.txt"),
                             header = TRUE, sep = "\t")
      muts_sbs$chrom <- gsub("chr","",muts_sbs$chrom)
      
      # subset for clocklike mutations
      muts_sbs15 <- subset(muts_sbs, Signature %in% c("SBS1","SBS5"))
      muts_sbs15 <- muts_all[paste0(muts_all$chr,":",muts_all$start) %in% paste0(muts_sbs15$chrom,":",muts_sbs15$start),]
      
      tryCatch({time_sbs15 <- time_amplification(cn_data = bb,
                                                 multiplicity_data = mult,
                                                 mutation_data = muts_sbs15,
                                                 muts_type = "SBS1 and SBS5",
                                                 sample_id = samples$icgc_donor_id[i],
                                                 amplification_chrom = amp_chr,
                                                 amplification_start = amp_start,
                                                 amplification_stop = amp_stop,
                                                 is_WGD = amp_wgd,
                                                 genome = "hg19")
      
      output_tab_sbs15[i,1:46] <- time_sbs15}, error=function(e){})
    }
    
    
    
    tryCatch({time_ctcpg <- time_amplification(cn_data = bb,
                                               multiplicity_data = mult,
                                               mutation_data = muts_all,
                                               muts_type = "All",
                                               sample_id = samples$icgc_donor_id[i],
                                               amplification_chrom = amp_chr,
                                               amplification_start = amp_start,
                                               amplification_stop = amp_stop,
                                               is_WGD = amp_wgd,
                                               genome = "hg19")
    
    output_tab_ctcpg[i,1:46] <- time_ctcpg}, error=function(e){})
    
    
  }
  print(i)
}



write.table(output_tab_all_muts, paste("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/20230103_AmplificationTimeR_BRCA_OV_PRAD/","PCAWG_MYC_BRCA-UK_OV-AU_PRAD-UK","_AmplificationTimeR_all_muts_",Sys.Date(),".txt", sep = ""),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(output_tab_sbs15, paste("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/20230103_AmplificationTimeR_BRCA_OV_PRAD/","PCAWG_MYC_BRCA-UK_OV-AU_PRAD-UK","_AmplificationTimeR_SBS1_5_muts_",Sys.Date(),".txt", sep = ""),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(output_tab_ctcpg, paste("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/20230103_AmplificationTimeR_BRCA_OV_PRAD/","PCAWG_MYC_BRCA-UK_OV-AU_PRAD-UK","_AmplificationTimeR_C_T_CpG_",Sys.Date(),".txt", sep = ""),
            sep = "\t", row.names = FALSE, quote = FALSE)

writeLines(capture.output(sessionInfo()), paste0("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/20230103_AmplificationTimeR_BRCA_OV_PRAD/AmplificationTimeR_run_sessionInfo_",Sys.Date(),".txt"))
