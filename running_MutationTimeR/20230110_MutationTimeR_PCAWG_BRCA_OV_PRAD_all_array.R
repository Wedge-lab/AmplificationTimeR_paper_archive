# Code to calculate gain timing using MutationTimeR
# For AmplificationTimeR manuscript
# Run as a batch/array job, because otherwise it takes forever

# Maria Jakobsdottir
# 29.11.22

# Run on PCAWG BRCA-UK, OV-AU, and PRAD-UK
# Use all mutations, C>T at CpG, and SBS1 and SBS5

### Libraries
library(MutationTimeR)
library(Repitools)
library(tidyr)

# sample name tumor_wgs_aliquot_id
args = commandArgs(trailingOnly=TRUE)
file_sample <- args[1]

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

selected_sample <- subset(samples, tumor_wgs_aliquot_id == file_sample)

if(all(file.exists(paste0("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/Battenberg/rerun005_final_subclones/",
                          selected_sample$tumor_wgs_aliquot_id,"_subclones.txt")),
       file.exists(paste0("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/final_consensus_12oct_passonly/snv_mnv/",
                          selected_sample$tumor_wgs_aliquot_id,".consensus.20160830.somatic.snv_mnv.vcf.gz")),
       file.exists(paste0("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/Battenberg/rerun005_final_subclones/",
                          selected_sample$tumor_wgs_aliquot_id,"_rho_and_psi.txt"))
       )){
  
  # Read in data
  bb <- read.delim(paste0("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/Battenberg/rerun005_final_subclones/",
                          selected_sample$tumor_wgs_aliquot_id,"_subclones.txt"),
                   sep = "\t", header = T)
  vcf <- readVcf(paste0("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/final_consensus_12oct_passonly/snv_mnv/",
                        selected_sample$tumor_wgs_aliquot_id,".consensus.20160830.somatic.snv_mnv.vcf.gz"))

  purity <- read.delim(paste0("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/Battenberg/rerun005_final_subclones/",
                              selected_sample$tumor_wgs_aliquot_id,"_rho_and_psi.txt"),
                       sep = "\t", header = T)
  tab <- TabixFile(paste0("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/final_consensus_12oct_passonly/snv_mnv/",
                          selected_sample$tumor_wgs_aliquot_id,".consensus.20160830.somatic.snv_mnv.vcf.gz"))
  
  
  ### Processing Battenberg input
  # separate out subclonal segments
  bb_subclonal <- subset(bb, frac1_A < 1)
  
  
  
  # create iranges object for all 1_A solutions
  bb_iranges <- IRanges(start = bb$startpos, end = bb$endpos)
  bb_rle <- Rle(bb$chr)
  bb_ranges <- GRanges(bb_rle, bb_iranges, major_cn=bb$nMaj1_A , minor_cn=bb$nMin1_A, clonal_frequency=(purity$rho[1] * bb$frac1_A)) # Copy number segments, needs columns  major_cn, minor_cn and clonal_frequency of each segment
  
  # create iranges object for all 2_A solutions (subclones)
  bb_iranges_sub <- IRanges(start = bb_subclonal$startpos, end = bb_subclonal$endpos)
  bb_rle_sub <- Rle(bb_subclonal$chr)
  bb_ranges_sub <- GRanges(bb_rle_sub, bb_iranges_sub, major_cn=bb_subclonal$nMaj2_A , minor_cn=bb_subclonal$nMin2_A, clonal_frequency=(purity$rho[1] * bb_subclonal$frac2_A))
  
  # combine iranges objects
  bb_ranges_all <- append(bb_ranges, bb_ranges_sub)
  
  
  
  ### Timing segments
  ## All mutations
  bb_ranges_all_all <- bb_ranges_all
  mt <- mutationTime(vcf, bb_ranges_all_all)
  
  mcols(bb_ranges_all_all) <- cbind(mcols(bb_ranges_all_all),mt$T)
  
  results_df_all <- as.data.frame(bb_ranges_all_all)
  results_df_all <- results_df_all[,-9]
  
  write.table(results_df_all, paste("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/MutationTimeR_BRCA_OV_PRAD/",
                                    selected_sample$icgc_donor_id,"_",
                                    selected_sample$tumor_wgs_aliquot_id,"_MutationTimeR_all_muts_full_",
                                    Sys.Date(),".txt", sep = ""),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  ####
  # Mutation subtypes
  if(file.exists(paste0("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/single_mutation_signature_attribution/",
                        selected_sample$icgc_donor_id,"_",
                        selected_sample$tumor_wgs_aliquot_id,"_individual_mutation_signature_attributions_2022-11-21.txt"))){
    muts <- read.delim(paste0("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/single_mutation_signature_attribution/",
                              selected_sample$icgc_donor_id,"_",
                              selected_sample$tumor_wgs_aliquot_id,"_individual_mutation_signature_attributions_2022-11-21.txt"),
                       header = TRUE, sep = "\t")
    muts$chrom <- gsub("chr","",muts$chrom)
    
    ### Processing vcf
    sbs15_muts <- subset(muts, Signature %in% c("SBS1","SBS5"))
    
    # make range object for sbs1/5 mutations
    sbs15_ranges <- GRanges(seqnames=sbs15_muts$chrom, 
                            ranges=IRanges(start=sbs15_muts$start,
                                           end=sbs15_muts$end))
    # read vcf for only sbs1/5 mutations
    vcf_sbs15 <- readVcf(tab, "hg19", param=sbs15_ranges)
    
    ctcpg_muts <- subset(muts, Mutation.Type == "C>T" & Mutation.Subtype %in% c("ACG","TCG","GCG","CCG"))
    
    # make range object for sbs1/5 mutations
    ctcpg_ranges <- GRanges(seqnames=ctcpg_muts$chrom, 
                            ranges=IRanges(start=ctcpg_muts$start,
                                           end=ctcpg_muts$end))
    # read vcf for only sbs1/5 mutations
    vcf_ctcpg <- readVcf(tab, "hg19", param=ctcpg_ranges)
    
    ## SBS1 and SBS5 mutations
    bb_ranges_all_sbs15 <- bb_ranges_all
    mt_sbs15 <- mutationTime(vcf_sbs15, bb_ranges_all_sbs15)
    
    mcols(bb_ranges_all_sbs15) <- cbind(mcols(bb_ranges_all_sbs15),mt_sbs15$T)
    
    results_df_sbs15 <- as.data.frame(bb_ranges_all_sbs15)
    results_df_sbs15 <- results_df_sbs15[,-9]
    
    write.table(results_df_sbs15, paste("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/MutationTimeR_BRCA_OV_PRAD/",
                                        selected_sample$icgc_donor_id,"_",
                                        selected_sample$tumor_wgs_aliquot_id,"_MutationTimeR_SBS1_5_full_",
                                        Sys.Date(),".txt", sep = ""),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    ## C>T at CpG mutations
    bb_ranges_all_ctcpg <- bb_ranges_all
    mt_ctcpg <- mutationTime(vcf_ctcpg, bb_ranges_all_ctcpg)
    
    mcols(bb_ranges_all_ctcpg) <- cbind(mcols(bb_ranges_all_ctcpg),mt_ctcpg$T)
    
    results_df_ctcpg <- as.data.frame(bb_ranges_all_ctcpg)
    results_df_ctcpg <- results_df_ctcpg[,-9]
    
    write.table(results_df_ctcpg, paste("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/MutationTimeR_BRCA_OV_PRAD/",
                                        selected_sample$icgc_donor_id,"_",
                                        selected_sample$tumor_wgs_aliquot_id,"_MutationTimeR_CT_at_CpG_full_",
                                        Sys.Date(),".txt", sep = ""),
                sep = "\t", row.names = FALSE, quote = FALSE)
  }
  }
  
  


writeLines(capture.output(sessionInfo()), paste0("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/MutationTimeR_BRCA_OV_PRAD/sessionInfo_",Sys.Date(),".txt"))

