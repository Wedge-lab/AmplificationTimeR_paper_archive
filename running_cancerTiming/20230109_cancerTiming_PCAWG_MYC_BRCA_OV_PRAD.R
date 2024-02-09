# Code to calculate gain timing using cancerTiming
# For AmplificationTimeR manuscript

# Maria Jakobsdottir
# 09.01.23

# Run on PCAWG BRCA-UK, OV-AU, and PRAD-UK
# Use all mutations, C>T at CpG, and SBS1 and SBS5

### Libraries
library(cancerTiming)
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
# Prepare output tables
output_table <- as.data.frame(matrix(data = NA, nrow = nrow(samples), ncol = 51))
colnames(output_table) <- c("sample", "highest_gain_class","t1","t2","t3","t4","t5","t6","t7","t8","t9","t10",
                            "t1_2.5","t2_2.5","t3_2.5","t4_2.5","t5_2.5","t6_2.5","t7_2.5","t8_2.5","t9_2.5","t10_2.5",
                            "t1_97.5","t2_97.5","t3_97.5","t4_97.5","t5_97.5","t6_97.5","t7_97.5","t8_97.5","t9_97.5","t10_97.5",
                            "q1","q2","q3","q4","q5","q6","q7","q8","q9","q10",
                            "Total_muts","Qualifying_muts","Success","Fail_reason","Method","Identifiable_history", 
                            "tumor_wgs_aliquot_id","tumor_wgs_icgc_specimen_id","wgd_status")
output_table$sample <- samples$icgc_donor_id
output_table$tumor_wgs_aliquot_id <- samples$tumor_wgs_aliquot_id
output_table$tumor_wgs_icgc_specimen_id <- samples$tumor_wgs_icgc_specimen_id
output_table$wgd_status <- samples$wgd_status

output_table_fullMLE_ident <- output_table
output_table_fullMLE_Nident <- output_table
output_table_fullMLE_ident_ctcpg <- output_table
output_table_fullMLE_Nident_ctcpg <- output_table
output_table_fullMLE_ident_sbs <- output_table
output_table_fullMLE_Nident_sbs <- output_table


output_table_fullMLE_ident$Identifiable_history <- "TRUE"
output_table_fullMLE_Nident$Identifiable_history <- "FALSE"
output_table_fullMLE_ident_ctcpg$Identifiable_history <- "TRUE"
output_table_fullMLE_Nident_ctcpg$Identifiable_history <- "FALSE"
output_table_fullMLE_ident_sbs$Identifiable_history <- "TRUE"
output_table_fullMLE_Nident_sbs$Identifiable_history <- "FALSE"

output_table_fullMLE_ident$Method <- "fullMLE"
output_table_fullMLE_Nident$Method <- "fullMLE"
output_table_fullMLE_ident_ctcpg$Method <- "fullMLE"
output_table_fullMLE_Nident_ctcpg$Method <- "fullMLE"
output_table_fullMLE_ident_sbs$Method <- "fullMLE"
output_table_fullMLE_Nident_sbs$Method <- "fullMLE"


for(i in 1:nrow(samples)){
  if(all(c(file.exists(paste0("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/Battenberg/rerun005_final_subclones/",
                              samples$tumor_wgs_aliquot_id[i],"_subclones.txt")),
           file.exists(paste0("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/DPClust/",
                              samples$tumor_wgs_aliquot_id[i],"/",
                              samples$tumor_wgs_aliquot_id[i],"_allDirichletProcessInfo.txt")),
           file.exists(paste0("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/Battenberg/rerun005_final_subclones/",
                              samples$tumor_wgs_aliquot_id[i],"_rho_and_psi.txt"))
           ))){
    
    # Read in files
    bb <- read.delim(paste0("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/Battenberg/rerun005_final_subclones/",
                            samples$tumor_wgs_aliquot_id[i],"_subclones.txt"),
                     sep = "\t", header = T)
    mult <- read.delim(paste0("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/DPClust/",
                              samples$tumor_wgs_aliquot_id[i],"/",
                              samples$tumor_wgs_aliquot_id[i],"_allDirichletProcessInfo.txt"))
    
    purity <- read.delim(paste0("/mnt/bmh01-rds/UoOxford_David_W/shared/PCAWG/Battenberg/rerun005_final_subclones/",
                                samples$tumor_wgs_aliquot_id[i],"_rho_and_psi.txt"),
                         sep = "\t", header = T)
    
    
    # work out normal contamination based on equation from David
    norm_contamination <- (2*(1-purity$rho[1]))/(purity$ploidy[1]*purity$rho[1] + 2*(1-purity$rho[1]))
    
    if(output_table$wgd_status[i] == "no_wgd"){
      wgd_status <- FALSE
    }else if(output_table$wgd_status[i] == "wgd"){
      wgd_status <- TRUE
    }
    
    # Paramters of amplification being investigated
    amp_chr <- 8 # MYC chromosome
    amp_start <- 128748315 # MYC start hg19
    amp_stop <-  128753680 # MYC stop hg19
    amp_wgd <- wgd_status
    amp_name <- samples$sample[i]
    
    # Subset copy number for segment spanning desired region
    tmp_bb <- subset(bb, chr == amp_chr &
                       startpos <= amp_stop &
                       endpos >= amp_start)
    # Subset multiplicity file for copy number segment (all mutations)
    tmp_mult <- subset(mult, chr == amp_chr & 
                         end >= tmp_bb$startpos & 
                         end <= tmp_bb$endpos )
    
    
    # subset mutation data
    onlyMuts <- tmp_mult
    onlyMuts$t_depth <- onlyMuts$WT.count+onlyMuts$mut.count
    
    # code to deal with individual signature attributions etc whether they exist or not
    if(file.exists(paste0("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/single_mutation_signature_attribution/",
                          samples$icgc_donor_id[i],"_",
                          samples$tumor_wgs_aliquot_id[i],"_individual_mutation_signature_attributions_2022-11-21.txt"))){
      muts <- read.delim(paste0("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/single_mutation_signature_attribution/",
                                samples$icgc_donor_id[i],"_",
                                samples$tumor_wgs_aliquot_id[i],"_individual_mutation_signature_attributions_2022-11-21.txt"),
                         header = TRUE, sep = "\t")
      muts$chrom <- gsub("chr","",muts$chrom)
      
      # subset for clocklike mutations
      muts_sbs <- subset(muts, Signature %in% c("SBS1","SBS5"))
      
      # Get C>T at CpG mutations
      muts_cpg <- subset(muts, Mutation.Type == "C>T" & Mutation.Subtype %in% c("ACG","TCG","GCG","CCG"))
      
      # Subset multiplicity file for copy number segment and SBS1/5 mutations
      tmp_mult_sbs <- merge(tmp_mult, muts_sbs, by.x = c("chr","end"), by.y = c("chrom","end"))
      
      # Subset multiplicity file for copy number segment and C>T at CpG mutations
      tmp_mult_ctcpg <- merge(tmp_mult, muts_cpg, by.x = c("chr","end"), by.y = c("chrom","end"))
      
      onlyMuts_sbs <- tmp_mult_sbs
      onlyMuts_sbs$t_depth <- onlyMuts_sbs$WT.count+onlyMuts_sbs$mut.count
      onlyMuts_ctcpg <- tmp_mult_ctcpg
      onlyMuts_ctcpg$t_depth <- onlyMuts_ctcpg$WT.count+onlyMuts_ctcpg$mut.count
    }else{
      onlyMuts_sbs <- as.data.frame(matrix(nrow = 0, ncol = 1))
      onlyMuts_ctcpg <- as.data.frame(matrix(nrow = 0, ncol = 1))
    }
    
    # work out the highest copy number state in region
    if(nrow(tmp_bb) > 0){
      # Get total copy number
      tmp_bb$n1A.sum <- (tmp_bb$nMaj1_A + tmp_bb$nMin1_A)
      tmp_bb$n2A.sum <- (tmp_bb$nMaj2_A + tmp_bb$nMin2_A)
      
      tmp_bb$n1A <- paste0(tmp_bb$nMaj1_A, "+" ,tmp_bb$nMin1_A)
      tmp_bb$n2A <- paste0(tmp_bb$nMaj2_A, "+" ,tmp_bb$nMin2_A)
      if(!is.na(tmp_bb$nMaj2_A)){
        highest_gain <- c("n1A","n2A")[which.max(c(tmp_bb$n1A.sum,tmp_bb$n2A.sum))]
      }else{
        highest_gain <- "n1A"
      }
    }else{
      highest_gain <- "no segments"
    }
    
    if(norm_contamination <0){
      highest_gain <- "purity > 1"
    }
    
    if(highest_gain %in% c("no segments","purity > 1")){
      print("Skipped because there are no segments or purity > 1")
    }else if(((highest_gain == "n1A") & (tmp_bb$n1A.sum <= 5) & (tmp_bb$n1A.sum > 2)) |
             ((highest_gain == "n1A") & tmp_bb$n1A == "2+0")){
      
      if(tmp_bb$n1A == "2+0"){
        event_type <- "CNLOH"
      }else{
        event_type <- "gain"
      }
      highest_gain_class <- tmp_bb$n1A
      
      # get history
      h_ident <- makeEventHistory(type = gsub("CN","",event_type), 
                                  copies = c(tmp_bb$nMaj1_A,tmp_bb$nMin1_A), 
                                  totalCopy = (tmp_bb$nMaj1_A + tmp_bb$nMin1_A), 
                                  onlyIdentifiable = TRUE)
      h_all <- makeEventHistory(type = gsub("CN","",event_type), 
                                copies = c(tmp_bb$nMaj1_A,tmp_bb$nMin1_A), 
                                totalCopy = (tmp_bb$nMaj1_A + tmp_bb$nMin1_A), 
                                onlyIdentifiable = FALSE)
      
      if(length(h_ident) > 0){
        ### Time ###
        ### All mutations ###
        ## fullMLE with identifiable history
        x_ident_fullMLE <- eventTiming(x = onlyMuts$mut.count, 
                                       m = onlyMuts$t_depth, 
                                       history = h_ident[[1]], 
                                       method = "fullMLE",
                                       totalCopy = (tmp_bb$nMaj1_A + tmp_bb$nMin1_A), 
                                       type = event_type, 
                                       normCont = norm_contamination,
                                       bootstrapCI = "nonparametric",
                                       minMutations = 3)
        
        for(p in 1:length(x_ident_fullMLE$pi)){
          output_table_fullMLE_ident[i,2+p] <- x_ident_fullMLE$pi[[p]]
          output_table_fullMLE_ident[i,12+p] <- x_ident_fullMLE$piCI[[p,1]]
          output_table_fullMLE_ident[i,22+p] <- x_ident_fullMLE$piCI[[p,2]]
        }
        output_table_fullMLE_ident$Total_muts[i] <- x_ident_fullMLE$summaryTable[1]
        output_table_fullMLE_ident$Qualifying_muts[i] <- x_ident_fullMLE$summaryTable[2]
        output_table_fullMLE_ident$Success[i] <- x_ident_fullMLE$success
        if(!is.null(x_ident_fullMLE$failReason)){
          output_table_fullMLE_ident$Fail_reason[i] <- x_ident_fullMLE$failReason
        }else{
          output_table_fullMLE_ident$Fail_reason[i] <- "None"
        }
        output_table_fullMLE_ident$highest_gain_class[i] <- highest_gain_class
        ## fullMLE with all histories
        x_all_fullMLE <- eventTiming(x = onlyMuts$mut.count, 
                                     m = onlyMuts$t_depth, 
                                     history = h_all[[1]], 
                                     method = "fullMLE",
                                     totalCopy = (tmp_bb$nMaj1_A + tmp_bb$nMin1_A), 
                                     type = event_type, 
                                     normCont = norm_contamination,
                                     bootstrapCI = "nonparametric",
                                     minMutations = 3)
        for(p in 1:length(x_all_fullMLE$pi)){
          output_table_fullMLE_Nident[i,2+p] <- x_all_fullMLE$pi[[p]]
          output_table_fullMLE_Nident[i,12+p] <- x_all_fullMLE$piCI[[p,1]]
          output_table_fullMLE_Nident[i,22+p] <- x_all_fullMLE$piCI[[p,2]]
        }
        output_table_fullMLE_Nident$Total_muts[i] <- x_all_fullMLE$summaryTable[1]
        output_table_fullMLE_Nident$Qualifying_muts[i] <- x_all_fullMLE$summaryTable[2]
        output_table_fullMLE_Nident$Success[i] <- x_all_fullMLE$success
        if(!is.null(x_all_fullMLE$failReason)){
          output_table_fullMLE_Nident$Fail_reason[i] <- x_all_fullMLE$failReason
        }else{
          output_table_fullMLE_Nident$Fail_reason[i] <- "None"
        }
        output_table_fullMLE_Nident$highest_gain_class[i] <- highest_gain_class
        if(nrow(onlyMuts_sbs) > 0){
          ### SNS1 and SBS5 mutations ###
          ## fullMLE with identifiable history
          x_ident_fullMLE_sbs <- eventTiming(x = onlyMuts_sbs$mut.count, 
                                             m = onlyMuts_sbs$t_depth, 
                                             history = h_ident[[1]], 
                                             method = "fullMLE",
                                             totalCopy = (tmp_bb$nMaj1_A + tmp_bb$nMin1_A), 
                                             type = event_type, 
                                             normCont = norm_contamination,
                                             bootstrapCI = "nonparametric",
                                             minMutations = 3)
          
          for(p in 1:length(x_ident_fullMLE_sbs$pi)){
            output_table_fullMLE_ident_sbs[i,2+p] <- x_ident_fullMLE_sbs$pi[[p]]
            output_table_fullMLE_ident_sbs[i,12+p] <- x_ident_fullMLE_sbs$piCI[[p,1]]
            output_table_fullMLE_ident_sbs[i,22+p] <- x_ident_fullMLE_sbs$piCI[[p,2]]
          }
          output_table_fullMLE_ident_sbs$Total_muts[i] <- x_ident_fullMLE_sbs$summaryTable[1]
          output_table_fullMLE_ident_sbs$Qualifying_muts[i] <- x_ident_fullMLE_sbs$summaryTable[2]
          output_table_fullMLE_ident_sbs$Success[i] <- x_ident_fullMLE_sbs$success
          if(!is.null(x_ident_fullMLE_sbs$failReason)){
            output_table_fullMLE_ident_sbs$Fail_reason[i] <- x_ident_fullMLE_sbs$failReason
          }else{
            output_table_fullMLE_ident_sbs$Fail_reason[i] <- "None"
          }
          output_table_fullMLE_ident_sbs$highest_gain_class[i] <- highest_gain_class
          ## fullMLE with all histories
          x_all_fullMLE_sbs <- eventTiming(x = onlyMuts_sbs$mut.count, 
                                           m = onlyMuts_sbs$t_depth, 
                                           history = h_all[[1]], 
                                           method = "fullMLE",
                                           totalCopy = (tmp_bb$nMaj1_A + tmp_bb$nMin1_A), 
                                           type = event_type, 
                                           normCont = norm_contamination,
                                           bootstrapCI = "nonparametric",
                                           minMutations = 3)
          for(p in 1:length(x_all_fullMLE_sbs$pi)){
            output_table_fullMLE_Nident_sbs[i,2+p] <- x_all_fullMLE_sbs$pi[[p]]
            output_table_fullMLE_Nident_sbs[i,12+p] <- x_all_fullMLE_sbs$piCI[[p,1]]
            output_table_fullMLE_Nident_sbs[i,22+p] <- x_all_fullMLE_sbs$piCI[[p,2]]
          }
          output_table_fullMLE_Nident_sbs$Total_muts[i] <- x_all_fullMLE_sbs$summaryTable[1]
          output_table_fullMLE_Nident_sbs$Qualifying_muts[i] <- x_all_fullMLE_sbs$summaryTable[2]
          output_table_fullMLE_Nident_sbs$Success[i] <- x_all_fullMLE_sbs$success
          if(!is.null(x_all_fullMLE_sbs$failReason)){
            output_table_fullMLE_Nident_sbs$Fail_reason[i] <- x_all_fullMLE_sbs$failReason
          }else{
            output_table_fullMLE_Nident_sbs$Fail_reason[i] <- "None"
          }
          output_table_fullMLE_Nident_sbs$highest_gain_class[i] <- highest_gain_class
        }
        if(nrow(onlyMuts_ctcpg) > 0){
          ### C>T at CpG mutations ###
          ## fullMLE with identifiable history
          x_ident_fullMLE_ctcpg <- eventTiming(x = onlyMuts_ctcpg$mut.count, 
                                               m = onlyMuts_ctcpg$t_depth, 
                                               history = h_ident[[1]], 
                                               method = "fullMLE",
                                               totalCopy = (tmp_bb$nMaj1_A + tmp_bb$nMin1_A), 
                                               type = event_type, 
                                               normCont = norm_contamination,
                                               bootstrapCI = "nonparametric",
                                               minMutations = 3)
          
          for(p in 1:length(x_ident_fullMLE_ctcpg$pi)){
            output_table_fullMLE_ident_ctcpg[i,2+p] <- x_ident_fullMLE_ctcpg$pi[[p]]
            output_table_fullMLE_ident_ctcpg[i,12+p] <- x_ident_fullMLE_ctcpg$piCI[[p,1]]
            output_table_fullMLE_ident_ctcpg[i,22+p] <- x_ident_fullMLE_ctcpg$piCI[[p,2]]
          }
          output_table_fullMLE_ident_ctcpg$Total_muts[i] <- x_ident_fullMLE_ctcpg$summaryTable[1]
          output_table_fullMLE_ident_ctcpg$Qualifying_muts[i] <- x_ident_fullMLE_ctcpg$summaryTable[2]
          output_table_fullMLE_ident_ctcpg$Success[i] <- x_ident_fullMLE_ctcpg$success
          if(!is.null(x_ident_fullMLE_ctcpg$failReason)){
            output_table_fullMLE_ident_ctcpg$Fail_reason[i] <- x_ident_fullMLE_ctcpg$failReason
          }else{
            output_table_fullMLE_ident_ctcpg$Fail_reason[i] <- "None"
          }
          output_table_fullMLE_ident_ctcpg$highest_gain_class[i] <- highest_gain_class
          ## fullMLE with all histories
          x_all_fullMLE_ctcpg <- eventTiming(x = onlyMuts_ctcpg$mut.count, 
                                             m = onlyMuts_ctcpg$t_depth, 
                                             history = h_all[[1]], 
                                             method = "fullMLE",
                                             totalCopy = (tmp_bb$nMaj1_A + tmp_bb$nMin1_A), 
                                             type = event_type, 
                                             normCont = norm_contamination,
                                             bootstrapCI = "nonparametric",
                                             minMutations = 3)
          for(p in 1:length(x_all_fullMLE_ctcpg$pi)){
            output_table_fullMLE_Nident_ctcpg[i,2+p] <- x_all_fullMLE_ctcpg$pi[[p]]
            output_table_fullMLE_Nident_ctcpg[i,12+p] <- x_all_fullMLE_ctcpg$piCI[[p,1]]
            output_table_fullMLE_Nident_ctcpg[i,22+p] <- x_all_fullMLE_ctcpg$piCI[[p,2]]
          }
          output_table_fullMLE_Nident_ctcpg$Total_muts[i] <- x_all_fullMLE_ctcpg$summaryTable[1]
          output_table_fullMLE_Nident_ctcpg$Qualifying_muts[i] <- x_all_fullMLE_ctcpg$summaryTable[2]
          output_table_fullMLE_Nident_ctcpg$Success[i] <- x_all_fullMLE_ctcpg$success
          if(!is.null(x_all_fullMLE_ctcpg$failReason)){
            output_table_fullMLE_Nident_ctcpg$Fail_reason[i] <- x_all_fullMLE_ctcpg$failReason
          }else{
            output_table_fullMLE_Nident_ctcpg$Fail_reason[i] <- "None"
          }
          output_table_fullMLE_Nident_ctcpg$highest_gain_class[i] <- highest_gain_class
        }
      }
      
      ############################################################################
      ############################################################################
    }else if(((highest_gain == "n2A") & (tmp_bb$n2A.sum <= 5) & (tmp_bb$n2A.sum > 2)) |
             ((highest_gain == "n2A") & tmp_bb$n2A == "2+0")){
      if(tmp_bb$n2A == "2+0"){
        event_type <- "CNLOH"
      }else{
        event_type <- "gain"
      }
      
      highest_gain_class <- tmp_bb$n2A
      # get history
      h_ident <- makeEventHistory(type = gsub("CN","",event_type), 
                                  copies = c(tmp_bb$nMaj2_A,tmp_bb$nMin2_A), 
                                  totalCopy = (tmp_bb$nMaj2_A + tmp_bb$nMin2_A), 
                                  onlyIdentifiable = TRUE)
      h_all <- makeEventHistory(type = gsub("CN","",event_type), 
                                copies = c(tmp_bb$nMaj2_A,tmp_bb$nMin2_A), 
                                totalCopy = (tmp_bb$nMaj2_A + tmp_bb$nMin2_A), 
                                onlyIdentifiable = FALSE)
      
      if(length(h_ident) > 0){
        ### Time ###
        ## fullMLE with identifiable history
        x_ident_fullMLE <- eventTiming(x = onlyMuts$mut.count, 
                                       m = onlyMuts$t_depth, 
                                       history = h_ident[[1]], 
                                       method = "fullMLE",
                                       totalCopy = (tmp_bb$nMaj2_A + tmp_bb$nMin2_A), 
                                       type = event_type, 
                                       normCont = norm_contamination,
                                       bootstrapCI = "nonparametric",
                                       minMutations = 3)
        for(p in 1:length(x_ident_fullMLE$pi)){
          output_table_fullMLE_ident[i,2+p] <- x_ident_fullMLE$pi[[p]]
          output_table_fullMLE_ident[i,12+p] <- x_ident_fullMLE$piCI[[p,1]]
          output_table_fullMLE_ident[i,22+p] <- x_ident_fullMLE$piCI[[p,2]]
        }
        output_table_fullMLE_ident$Total_muts[i] <- x_ident_fullMLE$summaryTable[1]
        output_table_fullMLE_ident$Qualifying_muts[i] <- x_ident_fullMLE$summaryTable[2]
        output_table_fullMLE_ident$Success[i] <- x_ident_fullMLE$success
        if(!is.null(x_ident_fullMLE$failReason)){
          output_table_fullMLE_ident$Fail_reason[i] <- x_ident_fullMLE$failReason
        }else{
          output_table_fullMLE_ident$Fail_reason[i] <- "None"
        }
        output_table_fullMLE_ident$highest_gain_class[i] <- highest_gain_class
        ## fullMLE with all histories
        x_all_fullMLE <- eventTiming(x = onlyMuts$mut.count, 
                                     m = onlyMuts$t_depth, 
                                     history = h_all[[1]], 
                                     method = "fullMLE",
                                     totalCopy = (tmp_bb$nMaj2_A + tmp_bb$nMin2_A), 
                                     type = event_type, 
                                     normCont = norm_contamination,
                                     bootstrapCI = "nonparametric",
                                     minMutations = 3)
        for(p in 1:length(x_all_fullMLE$pi)){
          output_table_fullMLE_Nident[i,2+p] <- x_all_fullMLE$pi[[p]]
          output_table_fullMLE_Nident[i,12+p] <- x_all_fullMLE$piCI[[p,1]]
          output_table_fullMLE_Nident[i,22+p] <- x_all_fullMLE$piCI[[p,2]]
        }
        output_table_fullMLE_Nident$Total_muts[i] <- x_all_fullMLE$summaryTable[1]
        output_table_fullMLE_Nident$Qualifying_muts[i] <- x_all_fullMLE$summaryTable[2]
        output_table_fullMLE_Nident$Success[i] <- x_all_fullMLE$success
        if(!is.null(x_all_fullMLE$failReason)){
          output_table_fullMLE_Nident$Fail_reason[i] <- x_all_fullMLE$failReason
        }else{
          output_table_fullMLE_Nident$Fail_reason[i] <- "None"
        }
        output_table_fullMLE_Nident$highest_gain_class[i] <- highest_gain_class
        
        if(nrow(onlyMuts_sbs) > 0){
          ### SNS1 and SBS5 mutations ###
          ## fullMLE with identifiable history
          x_ident_fullMLE_sbs <- eventTiming(x = onlyMuts_sbs$mut.count, 
                                             m = onlyMuts_sbs$t_depth, 
                                             history = h_ident[[1]], 
                                             method = "fullMLE",
                                             totalCopy = (tmp_bb$nMaj2_A + tmp_bb$nMin2_A), 
                                             type = event_type, 
                                             normCont = norm_contamination,
                                             bootstrapCI = "nonparametric",
                                             minMutations = 3)
          
          for(p in 1:length(x_ident_fullMLE_sbs$pi)){
            output_table_fullMLE_ident_sbs[i,2+p] <- x_ident_fullMLE_sbs$pi[[p]]
            output_table_fullMLE_ident_sbs[i,12+p] <- x_ident_fullMLE_sbs$piCI[[p,1]]
            output_table_fullMLE_ident_sbs[i,22+p] <- x_ident_fullMLE_sbs$piCI[[p,2]]
          }
          output_table_fullMLE_ident_sbs$Total_muts[i] <- x_ident_fullMLE_sbs$summaryTable[1]
          output_table_fullMLE_ident_sbs$Qualifying_muts[i] <- x_ident_fullMLE_sbs$summaryTable[2]
          output_table_fullMLE_ident_sbs$Success[i] <- x_ident_fullMLE_sbs$success
          if(!is.null(x_ident_fullMLE_sbs$failReason)){
            output_table_fullMLE_ident_sbs$Fail_reason[i] <- x_ident_fullMLE_sbs$failReason
          }else{
            output_table_fullMLE_ident_sbs$Fail_reason[i] <- "None"
          }
          output_table_fullMLE_ident_sbs$highest_gain_class[i] <- highest_gain_class
          ## fullMLE with all histories
          x_all_fullMLE_sbs <- eventTiming(x = onlyMuts_sbs$mut.count, 
                                           m = onlyMuts_sbs$t_depth, 
                                           history = h_all[[1]], 
                                           method = "fullMLE",
                                           totalCopy = (tmp_bb$nMaj2_A + tmp_bb$nMin2_A), 
                                           type = event_type, 
                                           normCont = norm_contamination,
                                           bootstrapCI = "nonparametric",
                                           minMutations = 3)
          for(p in 1:length(x_all_fullMLE_sbs$pi)){
            output_table_fullMLE_Nident_sbs[i,2+p] <- x_all_fullMLE_sbs$pi[[p]]
            output_table_fullMLE_Nident_sbs[i,12+p] <- x_all_fullMLE_sbs$piCI[[p,1]]
            output_table_fullMLE_Nident_sbs[i,22+p] <- x_all_fullMLE_sbs$piCI[[p,2]]
          }
          output_table_fullMLE_Nident_sbs$Total_muts[i] <- x_all_fullMLE_sbs$summaryTable[1]
          output_table_fullMLE_Nident_sbs$Qualifying_muts[i] <- x_all_fullMLE_sbs$summaryTable[2]
          output_table_fullMLE_Nident_sbs$Success[i] <- x_all_fullMLE_sbs$success
          if(!is.null(x_all_fullMLE_sbs$failReason)){
            output_table_fullMLE_Nident_sbs$Fail_reason[i] <- x_all_fullMLE_sbs$failReason
          }else{
            output_table_fullMLE_Nident_sbs$Fail_reason[i] <- "None"
          }
          output_table_fullMLE_Nident_sbs$highest_gain_class[i] <- highest_gain_class 
        }
        if(nrow(onlyMuts_ctcpg) > 0){
          ### C>T at CpG mutations ###
          ## fullMLE with identifiable history
          x_ident_fullMLE_ctcpg <- eventTiming(x = onlyMuts_ctcpg$mut.count, 
                                               m = onlyMuts_ctcpg$t_depth, 
                                               history = h_ident[[1]], 
                                               method = "fullMLE",
                                               totalCopy = (tmp_bb$nMaj2_A + tmp_bb$nMin2_A), 
                                               type = event_type, 
                                               normCont = norm_contamination,
                                               bootstrapCI = "nonparametric",
                                               minMutations = 3)
          
          for(p in 1:length(x_ident_fullMLE_ctcpg$pi)){
            output_table_fullMLE_ident_ctcpg[i,2+p] <- x_ident_fullMLE_ctcpg$pi[[p]]
            output_table_fullMLE_ident_ctcpg[i,12+p] <- x_ident_fullMLE_ctcpg$piCI[[p,1]]
            output_table_fullMLE_ident_ctcpg[i,22+p] <- x_ident_fullMLE_ctcpg$piCI[[p,2]]
          }
          output_table_fullMLE_ident_ctcpg$Total_muts[i] <- x_ident_fullMLE_ctcpg$summaryTable[1]
          output_table_fullMLE_ident_ctcpg$Qualifying_muts[i] <- x_ident_fullMLE_ctcpg$summaryTable[2]
          output_table_fullMLE_ident_ctcpg$Success[i] <- x_ident_fullMLE_ctcpg$success
          if(!is.null(x_ident_fullMLE_ctcpg$failReason)){
            output_table_fullMLE_ident_ctcpg$Fail_reason[i] <- x_ident_fullMLE_ctcpg$failReason
          }else{
            output_table_fullMLE_ident_ctcpg$Fail_reason[i] <- "None"
          }
          output_table_fullMLE_ident_ctcpg$highest_gain_class[i] <- highest_gain_class
          ## fullMLE with all histories
          x_all_fullMLE_ctcpg <- eventTiming(x = onlyMuts_ctcpg$mut.count, 
                                             m = onlyMuts_ctcpg$t_depth, 
                                             history = h_all[[1]], 
                                             method = "fullMLE",
                                             totalCopy = (tmp_bb$nMaj2_A + tmp_bb$nMin2_A), 
                                             type = event_type, 
                                             normCont = norm_contamination,
                                             bootstrapCI = "nonparametric",
                                             minMutations = 3)
          for(p in 1:length(x_all_fullMLE_ctcpg$pi)){
            output_table_fullMLE_Nident_ctcpg[i,2+p] <- x_all_fullMLE_ctcpg$pi[[p]]
            output_table_fullMLE_Nident_ctcpg[i,12+p] <- x_all_fullMLE_ctcpg$piCI[[p,1]]
            output_table_fullMLE_Nident_ctcpg[i,22+p] <- x_all_fullMLE_ctcpg$piCI[[p,2]]
          }
          output_table_fullMLE_Nident_ctcpg$Total_muts[i] <- x_all_fullMLE_ctcpg$summaryTable[1]
          output_table_fullMLE_Nident_ctcpg$Qualifying_muts[i] <- x_all_fullMLE_ctcpg$summaryTable[2]
          output_table_fullMLE_Nident_ctcpg$Success[i] <- x_all_fullMLE_ctcpg$success
          if(!is.null(x_all_fullMLE_ctcpg$failReason)){
            output_table_fullMLE_Nident_ctcpg$Fail_reason[i] <- x_all_fullMLE_ctcpg$failReason
          }else{
            output_table_fullMLE_Nident_ctcpg$Fail_reason[i] <- "None"
          }
          output_table_fullMLE_Nident_ctcpg$highest_gain_class[i] <- highest_gain_class
        }
        
      }  
    }
  }
  print(i)
}



write.table(output_table_fullMLE_ident, paste("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/cancerTiming_MYC_BRCA_OV_PRAD/","PCAWG_MYC_BRCA-UK_OV-AU_PRAD-UK","_cancerTiming_fullMLE_identifiable_all_muts_",Sys.Date(),".txt", sep = ""),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(output_table_fullMLE_Nident, paste("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/cancerTiming_MYC_BRCA_OV_PRAD/","PCAWG_MYC_BRCA-UK_OV-AU_PRAD-UK","_cancerTiming_fullMLE_all_all_muts_",Sys.Date(),".txt", sep = ""),
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(output_table_fullMLE_ident_ctcpg, paste("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/cancerTiming_MYC_BRCA_OV_PRAD/","PCAWG_MYC_BRCA-UK_OV-AU_PRAD-UK","_cancerTiming_fullMLE_identifiable_CT_at_CpG_",Sys.Date(),".txt", sep = ""),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(output_table_fullMLE_Nident_ctcpg, paste("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/cancerTiming_MYC_BRCA_OV_PRAD/","PCAWG_MYC_BRCA-UK_OV-AU_PRAD-UK","_cancerTiming_fullMLE_all_CT_at_CpG_",Sys.Date(),".txt", sep = ""),
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(output_table_fullMLE_ident_sbs, paste("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/cancerTiming_MYC_BRCA_OV_PRAD/","PCAWG_MYC_BRCA-UK_OV-AU_PRAD-UK","_cancerTiming_fullMLE_identifiable_SBS1_5_",Sys.Date(),".txt", sep = ""),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(output_table_fullMLE_Nident_sbs, paste("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/cancerTiming_MYC_BRCA_OV_PRAD/","PCAWG_MYC_BRCA-UK_OV-AU_PRAD-UK","_cancerTiming_fullMLE_all_SBS1_5_",Sys.Date(),".txt", sep = ""),
            sep = "\t", row.names = FALSE, quote = FALSE)

writeLines(capture.output(sessionInfo()), paste0("/mnt/bmh01-rds/UoOxford_David_W/b05055gj/PCAWG_analysis/cancerTiming_MYC_BRCA_OV_PRAD/sessionInfo_",Sys.Date(),".txt"))
