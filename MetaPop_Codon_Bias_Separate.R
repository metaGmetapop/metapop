options <- commandArgs(trailingOnly = T)

check_for_default <- integer(0)

directory_name <- which(grepl("-dir", options))
library_location <- which(grepl("-lib", options))
specified_genes <- which(grepl("-genes", options))


quitNow <- 0

if(length(directory_name) > 1 | identical(check_for_default, directory_name)){
  print("Metapop requires the name of the directory where the folders output from preprocessing exist. Metapop will exit.")
  quitNow <- 1
}else{
  directory_name <- options[directory_name + 1]
}

if(quitNow == 1){
  quit(save = "no")
}

setwd(directory_name)

if(file.exists("metapop_run_settings/run_settings.tsv")){
  run_parameters <- read.csv("metapop_run_settings/run_settings.tsv", sep = "\t")  
}else{
  print("Run parameters not found in the specified directory. This is fine, but you will have to supply -lib and -genes for codonbias to work.")
  print("-lib is the location of your R libraries")
  print("-genes is the path to the genes you wish to perform codon bias on.")
}


if(identical(specified_genes, check_for_default)){
  if(file.exists("metapop_run_settings/run_settings.tsv")){
    geneFile <- as.character(run_parameters$setting[run_parameters$parameter=="Genes"])
  }else{
    print("You must supply a genes file in FASTA format with -genes or set the directory to the location of a metapop_preprocess run. Exiting.")
    quit(save = "no")
  }
  
}else{
  geneFile <- options[specified_genes+1]
  
  if(!file.exists(geneFile)){
    print("Genes not found! Provide either an absolute path, or a path relative to the specifed directory. Exiting.")
    quit(save = "no")
  }
  
}

if(identical(library_location, check_for_default)){
  if(file.exists("metapop_run_settings/run_settings.tsv")){
    library_location <- as.character(run_parameters$setting[run_parameters$parameter=="Library Location"])
  }else{
    print("You must supply a library location with -lib or set the directory to the location of a metapop_preprocess run. Exiting.")
    quit(save = "no")
  }
}else{
  library_location <- options[library_location+1]
}

suppressMessages(suppressWarnings(library(stringr, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(Biostrings, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(data.table, lib.loc = library_location)))

{
  construct <- c("AAA", "AAT", "AAC", "AAG",
                 "ATA", "ATT", "ATC", "ATG",
                 "ACA", "ACT", "ACC", "ACG",
                 "AGA", "AGT", "AGC", "AGG",
                 
                 "TAA", "TAT", "TAC", "TAG",
                 "TTA", "TTT", "TTC", "TTG",
                 "TCA", "TCT", "TCC", "TCG",
                 "TGA", "TGT", "TGC", "TGG",
                 
                 "CAA", "CAT", "CAC", "CAG",
                 "CTA", "CTT", "CTC", "CTG",
                 "CCA", "CCT", "CCC", "CCG",
                 "CGA", "CGT", "CGC", "CGG",
                 
                 "GAA", "GAT", "GAC", "GAG",
                 "GTA", "GTT", "GTC", "GTG",
                 "GCA", "GCT", "GCC", "GCG",
                 "GGA", "GGT", "GGC", "GGG")
  proteinConstruct <- c("K", "N", "N", "K",
                        "I", "I", "I", "M",
                        "T", "T", "T", "T",
                        "R", "S", "S", "R",
                        
                        "STOP", "Y", "Y", "STOP",
                        "L", "F", "F", "L",
                        "S", "S", "S", "S",
                        "STOP", "C", "C", "W",
                        
                        "Q", "H", "H", "Q",
                        "L", "L", "L", "L",
                        "P", "P", "P", "P",
                        "R", "R", "R", "R",
                        
                        "E", "D", "D", "E",
                        "V", "V", "V", "V",
                        "A", "A", "A", "A",
                        "G", "G", "G", "G")
  
  expNConstruct<- c(8/3, 8/3, 8/3, 8/3,
                    7/3, 7/3, 7/3, 3,
                    2,2,2,2,
                    7/3,8/3,8/3,7/3,
                    
                    7/3, 8/3, 8/3,8/3,
                    7/3,8/3,8/3,7/3,
                    2,2,2,2,
                    8/3,8/3,8/3,3,
                    
                    8/3,8/3,8/3,8/3,
                    5/3,1,1,5/3,
                    2,2,2,2,
                    5/3,2,2,5/3,
                    
                    8/3,8/3,8/3,8/3,
                    2,2,2,2,
                    2,2,2,2,
                    2,2,2,2)
  
  expSConstruct <- 3-expNConstruct
}

genes <- readDNAStringSet(geneFile)

# Might as well do this here and now - prepares the sequences to be processed by awk outside of this file - the flat format is just better for this
codon_bias_input <- data.table(seq = unname(as.character(genes)))

genes <- names(genes)

genes <- data.table(contig_gene = unlist(lapply(genes, function(x){
  return(strsplit(x, split = " ")[[1]][1])
})))

genes[, parent_contig := gsub("_\\d+$", "", contig_gene)]

codon_bias_input$contig_gene <- genes$contig_gene
codon_bias_input$parent_contig <- genes$parent_contig

if(!dir.exists("metapop_codon_bias")){
  system("mkdir metapop_codon_bias")
}

fwrite(codon_bias_input, "metapop_codon_bias/codon_bias_sequences.tsv", sep = "\t", col.names = F)

codon_initialize <- paste(paste("codons[\"", construct, "\"]=0;", sep = ""), collapse = " ")
codon_output <- paste(paste("codons[\"", construct, "\"]", sep = ""), collapse="\"\\t\"")
  
# The direct system call to make the codon bias part happen should work, but just in case a system doesn't like issuing a script from a script, it'll be saved for the user to bash
out <- paste("cut -f1 -d$'\\t' codon_bias_sequences.tsv | awk '{gsub(/.{3}/,\"& \")}1' | awk 'BEGIN{", codon_initialize, "} { for (key in codons) codons[key]=0; for ( i=1; i<NF+1; i++ ) codons[$i]+=1; print", codon_output, "}' > preprocessed_genes.tsv")
writeLines(c(paste0("cd ", getwd(),"/metapop_codon_bias"), out), "metapop_codon_bias/launch_bash.sh")
system("bash metapop_codon_bias/launch_bash.sh")

genes <- fread(cmd=("cut -d'\t' -f2,3 metapop_codon_bias/codon_bias_sequences.tsv"), sep = "\t", header = F)

#Codon bias
if(dir.exists("metapop_codon_bias")){
  
  codon_data <- fread("metapop_codon_bias/preprocessed_genes.tsv", sep = "\t", header = F)
  colnames(codon_data) = construct
  
  #We have to get the proportion that each codon contributes to the AA by AA; split exaple on isoleucine ATT ATC ATA proportions;
  I <- codon_data$ATT + codon_data$ATC + codon_data$ATA
  L <- codon_data$CTT + codon_data$CTC + codon_data$CTA + codon_data$CTG + codon_data$TTA + codon_data$TTG
  V <- codon_data$GTT + codon_data$GTC + codon_data$GTA + codon_data$GTG
  FE <- codon_data$TTT + codon_data$TTC
  C <- codon_data$TGT + codon_data$TGC
  A <- codon_data$GCT + codon_data$GCC + codon_data$GCA + codon_data$GCG
  G <- codon_data$GGT + codon_data$GGC + codon_data$GGA + codon_data$GGG
  P <- codon_data$CCT + codon_data$CCC + codon_data$CCA + codon_data$CCG
  TH <- codon_data$ACT + codon_data$ACC + codon_data$ACA + codon_data$ACG
  S <- codon_data$TCT + codon_data$TCC + codon_data$TCA + codon_data$TCG + codon_data$AGT + codon_data$AGC
  Y <- codon_data$TAT + codon_data$TAC
  Q <- codon_data$CAA + codon_data$CAG
  N <- codon_data$AAT + codon_data$AAC
  H <- codon_data$CAT + codon_data$CAC
  E <- codon_data$GAA + codon_data$GAG
  D <- codon_data$GAT + codon_data$GAC
  K <- codon_data$AAA + codon_data$AAG
  R <- codon_data$CGT + codon_data$CGC + codon_data$CGA + codon_data$CGG + codon_data$AGA + codon_data$AGG
  ST<- codon_data$TAA + codon_data$TAG + codon_data$TGA
  
  #Now get their relative proportions
  
  rel_props <- data.table(codon_data$AAA/K,
                          codon_data$AAT/N,
                          codon_data$AAC/N,
                          codon_data$AAG/K,
                          
                          codon_data$ATA/I,
                          codon_data$ATT/I,
                          codon_data$ATC/I,
                          1,
                          
                          codon_data$ACA/TH,
                          codon_data$ACT/TH,
                          codon_data$ACC/TH,
                          codon_data$ACG/TH,
                          
                          codon_data$AGA/R,
                          codon_data$AGT/S,
                          codon_data$AGC/S,
                          codon_data$AGG/R,
                          
                          codon_data$TAA/ST,
                          codon_data$TAT/Y,
                          codon_data$TAC/Y,
                          codon_data$TAG/ST,
                          
                          codon_data$TTA/L,
                          codon_data$TTT/FE,
                          codon_data$TTC/FE,
                          codon_data$TTG/L,
                          
                          codon_data$TCA/S,
                          codon_data$TCT/S,
                          codon_data$TCC/S,
                          codon_data$TCG/S,
                          
                          codon_data$TGA/ST,
                          codon_data$TGT/C,
                          codon_data$TGC/C,
                          1,
                          
                          codon_data$CAA/Q,
                          codon_data$CAT/H,
                          codon_data$CAC/H,
                          codon_data$CAG/Q,
                          
                          codon_data$CTA/L,
                          codon_data$CTT/L,
                          codon_data$CTC/L,
                          codon_data$CTG/L,
                          
                          codon_data$CCA/P,
                          codon_data$CCT/P,
                          codon_data$CCC/P,
                          codon_data$CCG/P,
                          
                          codon_data$CGA/R,
                          codon_data$CGT/R,
                          codon_data$CGC/R,
                          codon_data$CGG/R,
                          
                          codon_data$GAA/E,
                          codon_data$GAT/D,
                          codon_data$GAC/D,
                          codon_data$GAG/E,
                          
                          codon_data$GTA/V,
                          codon_data$GTT/V,
                          codon_data$GTC/V,
                          codon_data$GTG/V,
                          
                          codon_data$GCA/A,
                          codon_data$GCT/A,
                          codon_data$GCC/A,
                          codon_data$GCG/A,
                          
                          codon_data$GGA/G,
                          codon_data$GGT/G,
                          codon_data$GGC/G,
                          codon_data$GGG/G)
  
  colnames(rel_props) = construct
  
  codon_data$parent_contig <- genes$V2
  
  is.nan.data.frame <- function(x)
    do.call(cbind, lapply(x, is.nan))
  
  rel_props[is.nan(rel_props)] <- 0
  
  rel_props$parent_contig <- genes$V2
  
  setkey(codon_data, "parent_contig")
  mean_genes <- codon_data[, lapply(.SD, sum), key(codon_data)]
  
  I <- mean_genes$ATT + mean_genes$ATC + mean_genes$ATA
  L <- mean_genes$CTT + mean_genes$CTC + mean_genes$CTA + mean_genes$CTG + mean_genes$TTA + mean_genes$TTG
  V <- mean_genes$GTT + mean_genes$GTC + mean_genes$GTA + mean_genes$GTG
  FE <- mean_genes$TTT + mean_genes$TTC
  C <- mean_genes$TGT + mean_genes$TGC
  A <- mean_genes$GCT + mean_genes$GCC + mean_genes$GCA + mean_genes$GCG
  G <- mean_genes$GGT + mean_genes$GGC + mean_genes$GGA + mean_genes$GGG
  P <- mean_genes$CCT + mean_genes$CCC + mean_genes$CCA + mean_genes$CCG
  TH <- mean_genes$ACT + mean_genes$ACC + mean_genes$ACA + mean_genes$ACG
  S <- mean_genes$TCT + mean_genes$TCC + mean_genes$TCA + mean_genes$TCG + mean_genes$AGT + mean_genes$AGC
  Y <- mean_genes$TAT + mean_genes$TAC
  Q <- mean_genes$CAA + mean_genes$CAG
  N <- mean_genes$AAT + mean_genes$AAC
  H <- mean_genes$CAT + mean_genes$CAC
  E <- mean_genes$GAA + mean_genes$GAG
  D <- mean_genes$GAT + mean_genes$GAC
  K <- mean_genes$AAA + mean_genes$AAG
  R <- mean_genes$CGT + mean_genes$CGC + mean_genes$CGA + mean_genes$CGG + mean_genes$AGA + mean_genes$AGG
  ST<- mean_genes$TAA + mean_genes$TAG + mean_genes$TGA
  
  
  rel_props_mean_gene <- data.table(mean_genes$AAA/K,
                                    mean_genes$AAT/N,
                                    mean_genes$AAC/N,
                                    mean_genes$AAG/K,
                                    
                                    mean_genes$ATA/I,
                                    mean_genes$ATT/I,
                                    mean_genes$ATC/I,
                                    1,
                                    
                                    mean_genes$ACA/TH,
                                    mean_genes$ACT/TH,
                                    mean_genes$ACC/TH,
                                    mean_genes$ACG/TH,
                                    
                                    mean_genes$AGA/R,
                                    mean_genes$AGT/S,
                                    mean_genes$AGC/S,
                                    mean_genes$AGG/R,
                                    
                                    mean_genes$TAA/ST,
                                    mean_genes$TAT/Y,
                                    mean_genes$TAC/Y,
                                    mean_genes$TAG/ST,
                                    
                                    mean_genes$TTA/L,
                                    mean_genes$TTT/FE,
                                    mean_genes$TTC/FE,
                                    mean_genes$TTG/L,
                                    
                                    mean_genes$TCA/S,
                                    mean_genes$TCT/S,
                                    mean_genes$TCC/S,
                                    mean_genes$TCG/S,
                                    
                                    mean_genes$TGA/ST,
                                    mean_genes$TGT/C,
                                    mean_genes$TGC/C,
                                    1,
                                    
                                    mean_genes$CAA/Q,
                                    mean_genes$CAT/H,
                                    mean_genes$CAC/H,
                                    mean_genes$CAG/Q,
                                    
                                    mean_genes$CTA/L,
                                    mean_genes$CTT/L,
                                    mean_genes$CTC/L,
                                    mean_genes$CTG/L,
                                    
                                    mean_genes$CCA/P,
                                    mean_genes$CCT/P,
                                    mean_genes$CCC/P,
                                    mean_genes$CCG/P,
                                    
                                    mean_genes$CGA/R,
                                    mean_genes$CGT/R,
                                    mean_genes$CGC/R,
                                    mean_genes$CGG/R,
                                    
                                    mean_genes$GAA/E,
                                    mean_genes$GAT/D,
                                    mean_genes$GAC/D,
                                    mean_genes$GAG/E,
                                    
                                    mean_genes$GTA/V,
                                    mean_genes$GTT/V,
                                    mean_genes$GTC/V,
                                    mean_genes$GTG/V,
                                    
                                    mean_genes$GCA/A,
                                    mean_genes$GCT/A,
                                    mean_genes$GCC/A,
                                    mean_genes$GCG/A,
                                    
                                    mean_genes$GGA/G,
                                    mean_genes$GGT/G,
                                    mean_genes$GGC/G,
                                    mean_genes$GGG/G,
                                    parent_contig = mean_genes$parent_contig)
  
  colnames(rel_props_mean_gene)[1:64] = construct
  rel_props_mean_gene[is.nan(rel_props_mean_gene)] <- 0
  
  rel_props_mean_gene <- data.frame(rel_props_mean_gene[match(rel_props$parent_contig, rel_props_mean_gene$parent_contig), 1:64])
  
  distance <- sqrt(rowSums((as.matrix(rel_props[,1:64], ncol = 64)-as.matrix(rel_props_mean_gene, ncol = 64))^2))
  
  fwrite(rel_props, "metapop_codon_bias/codon_usage_proportions.tsv", sep = "\t")
  
  rm(mean_genes, rel_props, rel_props_mean_gene)
  
  colnames(genes) = c("gene", "parent_contig")
  genes$euc_dist <- distance
  
  setkey(genes, "parent_contig")
  genes_iqr <- genes[, list((quantile(euc_dist, .75)-quantile(euc_dist, .25)), mean(euc_dist)), by = key(genes)]
  colnames(genes_iqr)[2:3] = c("iqr","mu")
  
  genes_iqr[,outlier_threshold_distance := mu+1.5*iqr]
  
  genes$outlier_status <- ifelse(genes$euc_dist>genes_iqr$outlier_threshold_distance[match(genes$parent_contig, genes_iqr$parent_contig)], "outlier", "not outlier")  
  
  fwrite(genes, "metapop_codon_bias/gene_euclidean_distances.tsv", sep = "\t")
  fwrite(genes_iqr, "metapop_codon_bias/gene_IQR_and_mean.tsv", sep = "\t")
  
}





