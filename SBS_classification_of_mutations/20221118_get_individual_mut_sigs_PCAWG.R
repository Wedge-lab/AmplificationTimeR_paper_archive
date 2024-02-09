#!/usr/bin/env Rscript

# R script for working out mutational signature of individual mutations 
# Designed for samples in PCAWG, can be re-worked to make it more compatible with other data
# !!! doesn't currently filter pass mutations because they were already filtered !!!
# Maria Jakobsdottir <maria.jakobsdottir@manchester.ac.uk>
# 18.11.22

##### Libraries #####
library(SomaticSignatures)
library(BSgenome.Hsapiens.UCSC.hg19)
##### Functions #####

## Mini functions ##
# inside function to split mutation with a ">"
split_mutation <- function(mutation){
  split_mut <- paste0(strsplit(mutation, split = character(0))[[1]][1],
                      ">",
                      strsplit(mutation, split = character(0))[[1]][2])
  return(split_mut)
}

# inside function to combine alteration with context
context_mutation <- function(alteration, context){
  context_mut <- paste0(strsplit(context, split = character(0))[[1]][1],
                        strsplit(alteration, split = character(0))[[1]][1],
                        strsplit(context, split = character(0))[[1]][3])
  return(context_mut)
}

# Function to get the highest mutation probability
highest_mut_prob <- function(mutation, prob_matrix){
  
  match_row <- prob_matrix[which((prob_matrix$Mutation.Type == mutation[1]) & (prob_matrix$Mutation.Subtype == mutation[2])),]
  # column 5 is the SBS1 column, anything before it is text
  
  highest_colname <- colnames(prob_matrix)[(which.max(match_row[,5:ncol(match_row)])+4)]
  highest_value <- match_row[,highest_colname]
  
  return(c(highest_colname,highest_value))
  
}

## Main function ##

#' @param sig_matrix_path Path to a csv SigProfiler matrix with signature probabilities for a set of samples.
#' Sample name column should be named "Sample". Also requires Mutation.Type Mutation.Subtype columns
#' @param sig_matrix_obj Data frame of SigProfiler matrix with signature probabilities for a set of samples.
#' Sample name column should be named "Sample". Also requires Mutation.Type Mutation.Subtype columns
#' @param mutations_path Path to mutations vcf file. Expected columns are: "CHROM",  "POS",     "ID",      "REF",     "ALT",     "QUAL",    "FILTER",  "INFO"
#' @param sample_id The sample ID used in the sig_matrix object so the mutations can be matched
#' @param output_id the sample ID that you want to use to name your output file. 
#' @param output_type must be either "data.frame" or "file"
#' @param output_path Path to output folder if writing to file
#' 
#' Use either sig_matrix_path or sig_matrix_obj, not both. 
get_mutation_signatures <- function(sig_matrix_path=NA, sig_matrix_obj=NA, mutations_path, sample_id, output_id, output_type, output_path=NA){
  # Read in or rename signature matrix
  if(!is.na(sig_matrix_path)){
    sig_matrix <- read.csv(sig_matrix_path, header = TRUE)
  }else if(!is.na(sig_matrix_obj)){
    sig_matrix <- sig_matrix_obj
  }else{
    stop("Either sig_matrix_path or sig_matrix_obj needs to be supplied")
  }
  
  if(!(output_type %in% c("data.frame","file"))){
    stop("output_type must be either 'data.frame' or 'file'.")
  }
  
  # subset signature matrix for relevant sample
  sample_sig_matrix <- subset(sig_matrix, Sample == sample_id)
  
  # Read in mutations vcf
  muts <- read.delim(mutations_path, header = FALSE, comment.char = "#")
  colnames(muts) <- c("CHROM",  "POS",     "ID",      "REF",     "ALT",     "QUAL",    "FILTER",  "INFO")
  
  #!!!!! IF USING VCFS THAT HAVEN'T BEEN FILTERED FOR PASSING VARIANTS YOU WILL NEED TO DO THIS HERE !!!!!
  
  ### Get trinucleotide context of mutations
  
  # Make sure we're only considering SNVs
  nucleotides <- c("A","T","G","C")
  muts <- muts[muts$REF %in% nucleotides & muts$ALT %in% nucleotides,]
  
  # create a VRanges object with your mutations to check against the reference genome
  if(nrow(muts[grep("chr", muts$CHROM),]) == 0){
    vr <- VariantAnnotation::VRanges(seqnames = paste0("chr",muts$CHROM),
                                     ranges = IRanges(muts$POS, muts$POS),
                                     ref = muts$REF, 
                                     alt = muts$ALT)
  }else{
    vr <- VariantAnnotation::VRanges(seqnames = muts$CHROM,
                                     ranges = IRanges(muts$POS, muts$POS),
                                     ref = muts$REF, 
                                     alt = muts$ALT)
  }
  
  # get mutation context using SomaticSignatures function
  context <- SomaticSignatures::mutationContext(vr, BSgenome.Hsapiens.UCSC.hg19, k=3, unify = TRUE)

  # turn this into a data frame 
  context_df <- as.data.frame(context)
  
  # Turn context data into input that matches SigProfiler output
  
  # split mutation for every mutation entry
  context_df$Mutation.Type <- lapply(context_df$alteration, split_mutation)
  
  # Get full context for all mutations
  context_df$Mutation.Subtype <- mapply(context_mutation, context_df$alteration, context = context_df$context)
  
  # Finally, get the signature and probability for each of the mutations
  sig_attributions <- as.data.frame(t(apply(context_df[,c('Mutation.Type','Mutation.Subtype')], 1, highest_mut_prob, prob_matrix = sample_sig_matrix)))
  colnames(sig_attributions) <- c("Signature","Probability")
  
  # tack onto the end of the mutation file
  context_df <- cbind(context_df, sig_attributions)
  
  #prepare final output
  output <- context_df[,c("seqnames","start","end","Mutation.Type","Mutation.Subtype","Signature","Probability")]
  colnames(output) <- c("chrom","start","end","Mutation.Type","Mutation.Subtype","Signature","Probability")
  
  # write output to a file or data.frame
  if(output_type == "data.frame"){
    return(output)
  }else if(output_type == "file"){
    output <- apply(output,2,as.character)
    write.table(output, file = paste0(output_path,output_id,"_individual_mutation_signature_attributions_",Sys.Date(),".txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
    
    print(paste0("Output saved to: ",output_path,output_id,"_individual_mutation_signature_attributions_",Sys.Date(),".txt"))
  }
}

##### Data #####
# sigprofiler matrix
sigpro_matrix <- read.csv("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/signatures/SigProfiler/SigProfilier_PCAWG_WGS_probabilities_SBS.csv", header = TRUE)

nrow(sigpro_matrix)/96

# metadata to match patient id to sample id (required to match  matrix to mutations)
metadata <- read.delim("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/release_may2016.v1.4.tsv", header = T, sep = "\t")

# some of the rows have multiple samples associated with the same donor
# for these rows we'll just keep the first sample listed
metadata$tumor_wgs_aliquot_id2 <- gsub(",.*","",metadata$tumor_wgs_aliquot_id)

metadata$mutation_file <- paste0("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/final_consensus_12oct_passonly/snv_mnv/",
                                 metadata$tumor_wgs_aliquot_id2,".consensus.20160830.somatic.snv_mnv.vcf.gz")
metadata$output_id <- paste0(metadata$icgc_donor_id,"_",metadata$tumor_wgs_aliquot_id2)

# check if mutation file exists and remove if it doesn't
metadata$mutations_exist <- lapply(metadata$mutation_file, file.exists)
metadata <- subset(metadata, mutations_exist == "TRUE")

# do the same to check that metadata samples match signature attribution samples
metadata$matched_to_signatures <- metadata$tumor_wgs_icgc_specimen_id %in% sigpro_matrix$Sample
metadata <- subset(metadata, matched_to_signatures == "TRUE")

nrow(metadata)

mapply(get_mutation_signatures, 
       mutations_path = metadata$mutation_file, 
       sample_id = metadata$tumor_wgs_icgc_specimen_id, 
       output_id = metadata$output_id,
       MoreArgs = list(sig_matrix_obj = sigpro_matrix,
                       output_type = "file", 
                       output_path = "/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/single_mutation_signature_attribution/"))

print("script finished successfully")

sessionInfo()