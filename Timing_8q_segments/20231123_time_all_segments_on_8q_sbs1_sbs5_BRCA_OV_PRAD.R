# Run AmplificationTimeR on all chromosome 8 segments
# for loop as there aren't that many samples

# install and load package
# devtools::install_github("Wedge-lab/AmplificationTimeR@suggested_updates", force = TRUE)
library(AmplificationTimeR)
library(BSgenome.Hsapiens.UCSC.hg19)

# Annotation
annotation <-  read.delim("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/evolution_and_heterogeneity/icgc_sample_annotations_summary_table.txt", header = T, sep = "\t")

# read in c8 segments
c8_segments <- read.delim("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/20231116_AmplificationTimeR_review_response/BRCA_OV_PRAD_all_8q_segments_2023-11-16.txt", sep = "\t", header = TRUE)

# prepare output table
output_table <- as.data.frame(matrix(nrow = 0, ncol = 47))
colnames(output_table) <- c("sample","region","highest_copy_number","event_order",
                            "num_mutations_used","clonality_status","flags",
                            "t_1","t_1_mean_bootstrap","t_1_lower_ci","t_1_upper_ci",
                            "t_2","t_2_mean_bootstrap","t_2_lower_ci","t_2_upper_ci",
                            "t_3","t_3_mean_bootstrap","t_3_lower_ci","t_3_upper_ci",
                            "t_4","t_4_mean_bootstrap","t_4_lower_ci","t_4_upper_ci",
                            "t_5","t_5_mean_bootstrap","t_5_lower_ci","t_5_upper_ci",
                            "t_6","t_6_mean_bootstrap","t_6_lower_ci","t_6_upper_ci",
                            "t_7","t_7_mean_bootstrap","t_7_lower_ci","t_7_upper_ci",
                            "t_8","t_8_mean_bootstrap","t_8_lower_ci","t_8_upper_ci",
                            "t_9","t_9_mean_bootstrap","t_9_lower_ci","t_9_upper_ci",
                            "t_10","t_10_mean_bootstrap","t_10_lower_ci","t_10_upper_ci")

# run amplificationtimer on each segment and add output to output table
for(i in 1:nrow(c8_segments)){#
  
  amp_chr <- c8_segments$chr[i]
  amp_start <- c8_segments$startpos[i]
  amp_stop <- c8_segments$endpos[i]
  
  amp_name <- gsub("_subclones.txt","",c8_segments$samplename[i])
  amp_name <- gsub("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/Battenberg/rerun005_final_subclones/","",amp_name)
  
  annotation_row <- which(annotation$tumour_aliquot_id == amp_name)
  
  if(annotation$wgd_status[annotation_row] == "no_wgd"){
    amp_wgd <- FALSE
  }else if(annotation$wgd_status[annotation_row] == "wgd"){
    amp_wgd <- TRUE
  }
  
  
  amp_cn <- read.delim(c8_segments$samplename[i], header = TRUE, sep = "\t")
  
  if(all(file.exists(paste0("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/DPClust/",amp_name,"/",
                            amp_name,"_allDirichletProcessInfo.txt")),
         file.exists(paste0("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/single_mutation_signature_attribution/",
                            annotation$icgc_donor_id[annotation_row],"_",
                            amp_name,"_individual_mutation_signature_attributions_2022-11-21.txt")),
         file.exists(paste0("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/final_consensus_12oct_passonly/snv_mnv/",
                            amp_name,".consensus.20160830.somatic.snv_mnv.vcf.gz")))){
    
    amp_mult <- read.delim(paste0("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/DPClust/",amp_name,"/",
                                  amp_name,"_allDirichletProcessInfo.txt"), header = TRUE, sep = "\t")
    
    amp_muts <- read.delim(paste0("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/final_consensus_12oct_passonly/snv_mnv/",
                                  amp_name,".consensus.20160830.somatic.snv_mnv.vcf.gz"), header = FALSE, sep = "\t", comment.char = "#")
    
    colnames(amp_muts) <- c("CHROM",  "POS",     "ID",      "REF",     "ALT",     "QUAL",    "FILTER",  "INFO")
    
    
    muts_all <- amp_muts[,c("CHROM","POS","POS","REF","ALT")]
    colnames(muts_all) <- c("chr","start","end","ref","alt")
    
    muts_sbs <- read.delim(paste0("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/single_mutation_signature_attribution/",
                                  annotation$icgc_donor_id[annotation_row],"_",
                                  amp_name,"_individual_mutation_signature_attributions_2022-11-21.txt"),
                           header = TRUE, sep = "\t")
    muts_sbs$chrom <- gsub("chr","",muts_sbs$chrom)
    
    # subset for clocklike mutations
    muts_sbs15 <- subset(muts_sbs, Signature %in% c("SBS1","SBS5"))
    muts_sbs15 <- muts_all[paste0(muts_all$chr,":",muts_all$start) %in% paste0(muts_sbs15$chrom,":",muts_sbs15$start),]
    # muts_sbs15$ref <- substr(muts_sbs15$Mutation.Type, 1,1)
    # muts_sbs15$alt <- substr(muts_sbs15$Mutation.Type, 3,3)
    # muts_sbs15$chr <- muts_sbs15$chrom
    # muts_sbs15 <- muts_sbs15[,c("chr","start","end","ref","alt")]
    
    tryCatch({
      amp_output <- time_amplification(cn_data = amp_cn,
                                       multiplicity_data = amp_mult,
                                       mutation_data = muts_sbs15,
                                       muts_type = "SBS1 and SBS5",
                                       amplification_chrom = amp_chr,
                                       amplification_start = amp_start,
                                       amplification_stop = amp_stop,
                                       sample_id = amp_name,
                                       is_WGD = amp_wgd,
                                       genome = "hg19")
      
      output_table <- rbind(output_table,amp_output)
    }, error=function(e){})
  }
}
# head(output_table)

# save original results
write.table(output_table, paste0("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/20231116_AmplificationTimeR_review_response/BRCA_OV_PRAD_AmplificationTimeR_C8_results_sbs1_sbs5_",Sys.Date(),".txt"), sep = "\t", row.names = FALSE, quote = FALSE)

###
### remove mismatched samples and only save one sample per patient ###
###

'%!in%' <- function(x,y)!('%in%'(x,y))
# files
metadata <- read.delim("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/release_may2016.v1.4.tsv", header = T, sep = "\t")
annotation <-  read.delim("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/evolution_and_heterogeneity/icgc_sample_annotations_summary_table.txt", header = T, sep = "\t")

# List of samples analysed
samples <- subset(metadata, dcc_project_code %in% c("BRCA-UK","BRCA-EU","BRCA-US",
                                                    "OV-AU","OV-US",
                                                    "PRAD-UK","PRAD-CA","PRAD-US"))
dim(samples)
# remove later aliquot IDs (after comma)
samples$tumor_wgs_aliquot_id <- gsub(",.*","",samples$tumor_wgs_aliquot_id)
dim(samples)
# merge with annotation to get extra info
samples <- merge(samples, annotation[,c("tumour_aliquot_id","wgd_status")], by.x = "tumor_wgs_aliquot_id", by.y = "tumour_aliquot_id")
dim(samples)

# mismatched samples
mismatched_samples <- read.delim("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/suspected_PCAWG_BRCA_OV_PRAD_mismatches_020523.txt", sep = "\t", header = TRUE)

# remove mismatched samples
sbs_muts_c8 <- output_table
sub_sbs_muts_c8 <- subset(sbs_muts_c8, sample %!in%  mismatched_samples$tumor_wgs_aliquot_id)

# merge with extra data and keep only one sample per patient
sub_sbs_muts_c8 <- merge(sub_sbs_muts_c8, samples, by.x = "sample", by.y = "tumor_wgs_aliquot_id")

# save data with mismatched samples removed and only one sample per patient
write.table(sub_sbs_muts_c8,
            paste0("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/20231116_AmplificationTimeR_review_response/BRCA_OV_PRAD_AmplificationTimeR_C8_results_sbs1_sbs5_",Sys.Date(),"_removed_mismatch_single_aliquot.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# save session info
writeLines(capture.output(sessionInfo()), paste0("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/20231116_AmplificationTimeR_review_response/BRCA_OV_PRAD_AmplificationTimeR_C8_results_sbs1_sbs5_sessionInfo",Sys.Date(),".txt"))

