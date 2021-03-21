options <- commandArgs(trailingOnly = T)

#universal
directory_name <- options[1]
threads <- as.numeric(options[2])
library_location <- options[3]

if(library_location == ""){
   library_location <- .libPaths()
}

ref_fasta <- options[4]
ref_genes <- options[5]

#Component specific
min_cov <- as.numeric(options[6])
min_dep <- as.numeric(options[7])

sub_samp <- as.numeric(options[8])

setwd(directory_name)


suppressMessages(suppressWarnings(library(doParallel, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(data.table, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(stringr, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(Biostrings, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(Rsamtools, lib.loc = library_location)))



#These get used later, so it's good to declare them now.
{
   
   #use the geom. mean of the depths across subsampled loci to determine harmonic number. Ceiling to be conservative.
   #I'll need to translate this set into a function.
   
   tajima_d_parameterization <- function(sample_size){
      harmonic <- sum(1/(c(1:(sample_size-1))))
      harmonic2 <- sum(1/(c((1:(sample_size-1))^2)))
      b1<-(sample_size+1)/(3*(sample_size-1))
      b2 <- (2*(sample_size^2+sample_size+3))/(9*sample_size*(sample_size-1))
      c1 <- b1-(1/harmonic)
      c2 <- b2-((sample_size+2)/(harmonic*sample_size)) + harmonic2/(harmonic^2)
      e1 <- c1/harmonic
      e2<-c2/(harmonic^2+harmonic2)
      
      return(list(harmonic, harmonic2, b1, b2, c1, c2, e1, e2))
   }
   
   harmony <- unlist(lapply(1:(sub_samp + 5), tajima_d_parameterization))
   
   harmonic <- harmony[c(T,F,F,F,F,F,F,F)]
   harmonic2 <- harmony[c(F,T,F,F,F,F,F,F)]
   b1 <- harmony[c(F,F,T,F,F,F,F,F)]
   b2 <- harmony[c(F,F,F,T,F,F,F,F)]
   c1 <- harmony[c(F,F,F,F,T,F,F,F)]
   c2 <- harmony[c(F,F,F,F,F,T,F,F)]
   e1 <- harmony[c(F,F,F,F,F,F,T,F)]
   e2 <- harmony[c(F,F,F,F,F,F,F,T)]
   
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
   
  reverseComplement <- unlist(lapply(construct,function(x){
     
     reversed <- intToUtf8(rev(utf8ToInt(x)))
     complimented <- chartr("ATCG", "TAGC", reversed)
     
  }))
  
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



genes <- fread("MetaPop/08.Codon_Bias/codon_bias_sequences.tsv", sep = "\t", header = F, select = c(3,4))
colnames(genes) = c("V1", "V2")

#Codon bias
if(dir.exists("MetaPop/08.Codon_Bias") & !file.exists("MetaPop/08.Codon_Bias/gene_IQR_and_mean.tsv")){
   
  codon_data <- fread("MetaPop/08.Codon_Bias/codon_counts_by_gene.tsv", sep = "\t", header = F)
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
  #Based on the standard codon <-> AA table, not microbial codon usage
   rel_props <- data.table(codon_data$AAA/K,
   codon_data$AAT/N,
   codon_data$AAC/N,
   codon_data$AAG/K,
  
   codon_data$ATA/I,
   codon_data$ATT/I,
   codon_data$ATC/I,
   #start codon
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
   #tryptophan
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
   
   fwrite(rel_props, "MetaPop/08.Codon_Bias/codon_usage_proportions.tsv", sep = "\t")
   
   rm(mean_genes, rel_props, rel_props_mean_gene)
   
   colnames(genes) = c("gene", "parent_contig")
   genes$euc_dist <- distance
   
   setkey(genes, "parent_contig")
   genes_iqr <- genes[, list((quantile(euc_dist, .75)-quantile(euc_dist, .25)), mean(euc_dist)), by = key(genes)]
   colnames(genes_iqr)[2:3] = c("iqr","mu")
      
   genes_iqr[,outlier_threshold_distance := mu+1.5*iqr]
   
   genes$outlier_status <- ifelse(genes$euc_dist>genes_iqr$outlier_threshold_distance[match(genes$parent_contig, genes_iqr$parent_contig)], "outlier", "not outlier")
   
   fwrite(genes, "MetaPop/08.Codon_Bias/gene_euclidean_distances.tsv", sep = "\t")
   fwrite(genes_iqr, "MetaPop/08.Codon_Bias/gene_IQR_and_mean.tsv", sep = "\t")
   
}

#subsample and microdiversity
# Occasionally the subsample will take 1 more than the sub_samp value indicates it should. 
# This happens because the particular counts and scaling factor sample_prop round up more than once

{
   # The normal data
   microdiv_data <- fread("MetaPop/07.Cleaned_SNPs/genic_snps.tsv", sep = "\t", header = T)
   
   setkeyv(microdiv_data, c("source", "contig"))
   
   non_genic_data <- fread("MetaPop/07.Cleaned_SNPs/non_genic_snps.tsv", sep = "\t", header = T)
   
   microdiv_data[,link := NULL]
   microdiv_data[, depth := a_ct+t_ct+c_ct+g_ct]
   # Can't calculate pi/etc on only one obs
   microdiv_data <- microdiv_data[microdiv_data$depth > 1]
   
   #This is collecting the number of SNPs in first, second, and 3rd position in codon over the genic SNPs, stratified by sample and contig. 
   codon_pos_count <- microdiv_data[,list(first_pos = sum(pos_in_codon==1), second_pos = sum(pos_in_codon==2), third_pos = sum(pos_in_codon==3)), by = key(microdiv_data)]
   
   microdiv_data$sample_prop <- sub_samp/microdiv_data$depth
   microdiv_data$sub_samp_a <- ifelse(microdiv_data$sample_prop < 1, round(microdiv_data$a_ct * microdiv_data$sample_prop), microdiv_data$a_ct)
   microdiv_data$sub_samp_t <- ifelse(microdiv_data$sample_prop < 1, round(microdiv_data$t_ct * microdiv_data$sample_prop), microdiv_data$t_ct)
   microdiv_data$sub_samp_c <- ifelse(microdiv_data$sample_prop < 1, round(microdiv_data$c_ct * microdiv_data$sample_prop), microdiv_data$c_ct)
   microdiv_data$sub_samp_g <- ifelse(microdiv_data$sample_prop < 1, round(microdiv_data$g_ct * microdiv_data$sample_prop), microdiv_data$g_ct)
   #Get the new depth per pos for pi calc
   microdiv_data[, sub_samp_depth := sub_samp_a+sub_samp_t+sub_samp_c+sub_samp_g]
   
   microdiv_data[,pi := 2*(((sub_samp_a * (sub_samp_t+sub_samp_c+sub_samp_g))+(sub_samp_t*(sub_samp_c+sub_samp_g))+(sub_samp_c*sub_samp_g))/((sub_samp_depth*sub_samp_depth)-sub_samp_depth))]
   
   non_genic_data[, depth := a_ct+t_ct+c_ct+g_ct]
   # Can't calculate pi/etc on only one obs
   non_genic_data <- non_genic_data[non_genic_data$depth > 1]
   non_genic_data$sample_prop <- sub_samp/non_genic_data$depth
   non_genic_data$sub_samp_a <- ifelse(non_genic_data$sample_prop < 1, round(non_genic_data$a_ct * non_genic_data$sample_prop), non_genic_data$a_ct)
   non_genic_data$sub_samp_t <- ifelse(non_genic_data$sample_prop < 1, round(non_genic_data$t_ct * non_genic_data$sample_prop), non_genic_data$t_ct)
   non_genic_data$sub_samp_c <- ifelse(non_genic_data$sample_prop < 1, round(non_genic_data$c_ct * non_genic_data$sample_prop), non_genic_data$c_ct)
   non_genic_data$sub_samp_g <- ifelse(non_genic_data$sample_prop < 1, round(non_genic_data$g_ct * non_genic_data$sample_prop), non_genic_data$g_ct)
   
   #Get the new depth per pos for pi calc
   non_genic_data[, sub_samp_depth := sub_samp_a+sub_samp_t+sub_samp_c+sub_samp_g]
   
   non_genic_data[,pi := 2*(((sub_samp_a * (sub_samp_t+sub_samp_c+sub_samp_g))+(sub_samp_t*(sub_samp_c+sub_samp_g))+(sub_samp_c*sub_samp_g))/((sub_samp_depth*sub_samp_depth)-sub_samp_depth))]
   
   
   # For theta and pi, we need the lengths of the genomes and the genes
   assembled_contigs <- readDNAStringSet(ref_fasta)
   assembled_contigs <- data.table(contig = names(assembled_contigs), length = nchar(assembled_contigs))
   
   assembled_contigs$contig <- unlist(lapply(assembled_contigs$contig, function(y){return(strsplit(y, split = " ")[[1]][1])}))
   
   gene_assembly <- fread("MetaPop/08.Codon_Bias/codon_bias_sequences.tsv", sep = "\t")
   colnames(gene_assembly) = c("seq", "OC", "gene", "contig")
   gene_assembly$length <- nchar(gene_assembly$seq)

   #Likely needs to be separate
   depth_file_names <- list.files(path = "MetaPop/03.Breadth_and_Depth", full.names = T)
   
   #returns a 0 length data.table if there are no passing contigs
   depth_info <- lapply(depth_file_names, function(x){
      tmp <- fread(x, sep = "\t")
      tmp$source <- substr(x, 30, nchar(x)-22)
      tmp <- tmp[V3 >= min_cov & V4 >= min_dep,]
      return(tmp)
   })
   
   names(depth_info) = substr(depth_file_names, 30, nchar(depth_file_names)-22)
   
   #removes the 0 length data.tables, i.e. samples with no passing contigs
   depth_info <- depth_info[unlist(lapply(depth_info, nrow)) > 0]
   
   #We need to keep the genes from these files for filling out the gene table later.
   all_contigs <- unique(unlist(lapply(depth_info, function(x){
      
      return(unique(x$V1))
      
   })))
   
   #Keeps only the genes belonging to contigs that passed at least one sample. The genes are shared by all samples, so there's no need to have these per sample
   gene_assembly <- gene_assembly[gene_assembly$contig %in% all_contigs,]
   
   #figure out what position in each gene the base falls.
   #genes are reverse complemented, so this should be fine
   pos_in_gene <- ((microdiv_data$codon-1)*3) + microdiv_data$pos_in_codon
   

   split_codons <- unname(lapply(gene_assembly$seq, function(x){
      str <- strsplit(x, split = "")[[1]]
      codons <- paste0(str[c(T,F,F)],
                       str[c(F,T,F)],
                       str[c(F,F,T)])
      return(codons)
   }))
   
   exp_N_S <- sapply(split_codons, function(x){
      
      N <- sum(expNConstruct[match(x, construct)])
      S <- sum(expSConstruct[match(x, construct)])
      
      return(list(N, S))
   })
   
   gene_assembly$expN <- as.numeric(exp_N_S[1,])
   gene_assembly$expS <- as.numeric(exp_N_S[2,])
   
   microdiv_data$original_codon <- mapply(function(codons, index){
      
      return(codons[index])
      
   }, split_codons[match(microdiv_data$contig_gene, gene_assembly$gene)], microdiv_data$codon)
   
   #Some codons are illegitimate on account of being on edges. Very few, but they need to be gone.
   microdiv_data <- microdiv_data[microdiv_data$original_codon %in% construct,]
   
   microdiv_data$original_AA <- proteinConstruct[match(microdiv_data$original_codon, construct)]
   
   #Handling SNPs - some positions have 2 or 3, and these need special handling.
   snps_by_row <- nchar(microdiv_data$snps)
   #Initial setup
   microdiv_data[,obsN := 0]
   microdiv_data[,obsS := 0]

   

   #This first bit handles the first, usually only, SNP per row and sets up the overall N/S output.
   first_snp <- substr(microdiv_data$snps, 1, 1)
   OG_cod <- microdiv_data$original_codon
   substr(OG_cod, microdiv_data$pos_in_codon, microdiv_data$pos_in_codon) <- first_snp
   OG_cod <- proteinConstruct[match(OG_cod, construct)] == microdiv_data$original_AA
   
   isA <- first_snp == "A"
   isT <- first_snp == "T"
   isC <- first_snp == "C"
   isG <- first_snp == "G"
   
   microdiv_data$obsN[(isA & !OG_cod)] <- microdiv_data$sub_samp_a[(isA & !OG_cod)]
   microdiv_data$obsN[(isT & !OG_cod)] <- microdiv_data$sub_samp_t[(isT & !OG_cod)]
   microdiv_data$obsN[(isC & !OG_cod)] <- microdiv_data$sub_samp_c[(isC & !OG_cod)]
   microdiv_data$obsN[(isG & !OG_cod)] <- microdiv_data$sub_samp_g[(isG & !OG_cod)]
   
   microdiv_data$obsS[(isA & OG_cod)] <- microdiv_data$sub_samp_a[(isA & OG_cod)]
   microdiv_data$obsS[(isT & OG_cod)] <- microdiv_data$sub_samp_t[(isT & OG_cod)]
   microdiv_data$obsS[(isC & OG_cod)] <- microdiv_data$sub_samp_c[(isC & OG_cod)]
   microdiv_data$obsS[(isG & OG_cod)] <- microdiv_data$sub_samp_g[(isG & OG_cod)]
   
   
   #Second SNP where applicable - here meaning multiple alleles on a single locus.
   second_snp <- substr(microdiv_data$snps[snps_by_row > 1], 2, 2)
   OG_cod <- microdiv_data$original_codon[snps_by_row > 1]
   
   substr(OG_cod, microdiv_data$pos_in_codon[snps_by_row > 1], microdiv_data$pos_in_codon[snps_by_row > 1]) <- second_snp
   OG_cod <- proteinConstruct[match(OG_cod, construct)] == microdiv_data$original_AA[snps_by_row > 1]
   
   isA <- second_snp == "A"
   isT <- second_snp == "T"
   isC <- second_snp == "C"
   isG <- second_snp == "G"
   
   microdiv_data$obsN[snps_by_row > 1][(isA & !OG_cod)] <- microdiv_data$obsN[snps_by_row > 1][(isA & !OG_cod)]+microdiv_data$sub_samp_a[snps_by_row > 1][(isA & !OG_cod)]
   microdiv_data$obsN[snps_by_row > 1][(isT & !OG_cod)] <- microdiv_data$obsN[snps_by_row > 1][(isT & !OG_cod)]+microdiv_data$sub_samp_t[snps_by_row > 1][(isT & !OG_cod)]
   microdiv_data$obsN[snps_by_row > 1][(isC & !OG_cod)] <- microdiv_data$obsN[snps_by_row > 1][(isC & !OG_cod)]+microdiv_data$sub_samp_c[snps_by_row > 1][(isC & !OG_cod)]
   microdiv_data$obsN[snps_by_row > 1][(isG & !OG_cod)] <- microdiv_data$obsN[snps_by_row > 1][(isG & !OG_cod)]+microdiv_data$sub_samp_g[snps_by_row > 1][(isG & !OG_cod)]
   
   microdiv_data$obsS[snps_by_row > 1][(isA & OG_cod)] <- microdiv_data$obsS[snps_by_row > 1][(isA & OG_cod)]+microdiv_data$sub_samp_a[snps_by_row > 1][(isA & OG_cod)]
   microdiv_data$obsS[snps_by_row > 1][(isT & OG_cod)] <- microdiv_data$obsS[snps_by_row > 1][(isT & OG_cod)]+microdiv_data$sub_samp_t[snps_by_row > 1][(isT & OG_cod)]
   microdiv_data$obsS[snps_by_row > 1][(isC & OG_cod)] <- microdiv_data$obsS[snps_by_row > 1][(isC & OG_cod)]+microdiv_data$sub_samp_c[snps_by_row > 1][(isC & OG_cod)]
   microdiv_data$obsS[snps_by_row > 1][(isG & OG_cod)] <- microdiv_data$obsS[snps_by_row > 1][(isG & OG_cod)]+microdiv_data$sub_samp_g[snps_by_row > 1][(isG & OG_cod)]
   
   
   #Ibid. for 2, but with 3 alt alleles on one locus. This is the third allele.
   third_snp <- substr(microdiv_data$snps[snps_by_row > 2], 3, 3)
   OG_cod <- microdiv_data$original_codon[snps_by_row > 2]
   
   substr(OG_cod, microdiv_data$pos_in_codon[snps_by_row > 2], microdiv_data$pos_in_codon[snps_by_row > 2]) <- third_snp
   OG_cod <- proteinConstruct[match(OG_cod, construct)] == microdiv_data$original_AA[snps_by_row > 2]
   
   isA <- third_snp == "A"
   isT <- third_snp == "T"
   isC <- third_snp == "C"
   isG <- third_snp == "G"
   
   microdiv_data$obsN[snps_by_row > 2][(isA & !OG_cod)] <- microdiv_data$obsN[snps_by_row > 2][(isA & !OG_cod)]+microdiv_data$sub_samp_a[snps_by_row > 2][(isA & !OG_cod)]
   microdiv_data$obsN[snps_by_row > 2][(isT & !OG_cod)] <- microdiv_data$obsN[snps_by_row > 2][(isT & !OG_cod)]+microdiv_data$sub_samp_t[snps_by_row > 2][(isT & !OG_cod)]
   microdiv_data$obsN[snps_by_row > 2][(isC & !OG_cod)] <- microdiv_data$obsN[snps_by_row > 2][(isC & !OG_cod)]+microdiv_data$sub_samp_c[snps_by_row > 2][(isC & !OG_cod)]
   microdiv_data$obsN[snps_by_row > 2][(isG & !OG_cod)] <- microdiv_data$obsN[snps_by_row > 2][(isG & !OG_cod)]+microdiv_data$sub_samp_g[snps_by_row > 2][(isG & !OG_cod)]
   
   microdiv_data$obsS[snps_by_row > 2][(isA & OG_cod)] <- microdiv_data$obsS[snps_by_row > 2][(isA & OG_cod)]+microdiv_data$sub_samp_a[snps_by_row > 2][(isA & OG_cod)]
   microdiv_data$obsS[snps_by_row > 2][(isT & OG_cod)] <- microdiv_data$obsS[snps_by_row > 2][(isT & OG_cod)]+microdiv_data$sub_samp_t[snps_by_row > 2][(isT & OG_cod)]
   microdiv_data$obsS[snps_by_row > 2][(isC & OG_cod)] <- microdiv_data$obsS[snps_by_row > 2][(isC & OG_cod)]+microdiv_data$sub_samp_c[snps_by_row > 2][(isC & OG_cod)]
   microdiv_data$obsS[snps_by_row > 2][(isG & OG_cod)] <- microdiv_data$obsS[snps_by_row > 2][(isG & OG_cod)]+microdiv_data$sub_samp_g[snps_by_row > 2][(isG & OG_cod)]
   
   rm(isA, isT, isC, isG, OG_cod, first_snp, second_snp, third_snp, snps_by_row)
   
   gm_mean = function(x){
      ceiling(exp(sum(log(x[x > 0]), na.rm=T) / length(x)))
   }
   
   
   {
   setkey(microdiv_data, "contig_gene")
   total_snps_genic <- microdiv_data[, list(length(unique(contig_pos)), unique(contig)), by = key(microdiv_data)]
   colnames(total_snps_genic)[2:3] = c("count", "contig")
   
   setkeyv(microdiv_data, c("contig_gene", "source"))
   gene_level_microdiv <- microdiv_data[, list(unique(contig), sum(pi, na.rm = T), length(contig), gm_mean(sub_samp_depth)), by=key(microdiv_data)]
   colnames(gene_level_microdiv)[3:6] = c("contig", "pi", "num_snps", "harmonic_number")
   
   #Shouldn't calculate microdiv if the depth is too low
   gene_level_microdiv <- gene_level_microdiv[harmonic_number > 2,]
   
   gene_level_microdiv$gene_len <- gene_assembly$length[match(gene_level_microdiv$contig_gene, gene_assembly$gene)]
   gene_level_microdiv$theta <- (gene_level_microdiv$num_snps/harmonic[gene_level_microdiv$harmonic_number])/(gene_level_microdiv$gene_len-(total_snps_genic$count[match(gene_level_microdiv$contig_gene, total_snps_genic$contig_gene)]-gene_level_microdiv$num_snps))
   
   #Modify pi to follow SNP count corrections exactly like theta
   gene_level_microdiv$pi <- gene_level_microdiv$pi/(gene_level_microdiv$gene_len-(total_snps_genic$count[match(gene_level_microdiv$contig_gene, total_snps_genic$contig_gene)]-gene_level_microdiv$num_snps))
   
   gene_level_microdiv$taj_D <- (gene_level_microdiv$pi-(gene_level_microdiv$num_snps/harmonic[gene_level_microdiv$harmonic_number]))/sqrt((e1[gene_level_microdiv$harmonic_number]*gene_level_microdiv$num_snps)+(e2[gene_level_microdiv$harmonic_number]*gene_level_microdiv$num_snps*(gene_level_microdiv$num_snps-1)))
   
   
   #check it out to make sure this works
   gene_level_microdiv <- gene_level_microdiv[, -6]
   
   saved_total_snps_genic <- total_snps_genic
   
   # Repeated at the level of the codon; all of the data gets incorporated as well as possible
   total_snps_genic$contig_gene = NULL
   
   setkey(non_genic_data, "contig")
   total_snps <- non_genic_data[, length(unique(contig_pos)), by=key(non_genic_data)]
   names(total_snps)[2] = "count"
   
   total_snps <- rbind(total_snps_genic, total_snps)

   setkey(total_snps, "contig")
   total_snps <- total_snps[, sum(count), by=key(total_snps)]
   colnames(total_snps)[2] = "count"
   
   cleaned_data <- data.table(contig_pos = microdiv_data$contig_pos, contig = microdiv_data$contig, pi = microdiv_data$pi, source = microdiv_data$source, dep = microdiv_data$sub_samp_depth)
   cleaned_data <- rbind(cleaned_data, data.table(contig_pos = non_genic_data$contig_pos, contig = non_genic_data$contig, pi=non_genic_data$pi, source = non_genic_data$source, dep = non_genic_data$sub_samp_depth))
   
   setkeyv(cleaned_data, c("contig", "source"))
   
   extant_names <- unique(cleaned_data$source)
   contigs_with_real_observations_by_sc <- cleaned_data[,list(list(unique(contig))), by = source]$V1
   names(contigs_with_real_observations_by_sc) = extant_names
   
   #Finds the contigs that passed the cov/depth filter but had no passing SNPs. 
   #This is added to the contig_level_microdiv as the has_snps = FALSE section.
   contigs_with_no_snps_by_sample <- lapply(names(depth_info), function(x){
      
      #if a sample has passing contigs, but no SNPs, we want to return all passing contigs in the depth file
      if(is.null(contigs_with_real_observations_by_sc[[x]])){
         return(depth_info[[x]])
      }else{
         #Else return all of the contigs that aren't in the data already - i.e. the SNPless ones.
         return(depth_info[[x]][!(V1 %in% contigs_with_real_observations_by_sc[[x]])])
      }
      
   })
   names(contigs_with_no_snps_by_sample) = names(depth_info)
   
   #Empty ones are fully covered already
   contigs_with_no_snps_by_sample <- contigs_with_no_snps_by_sample[unlist(lapply(contigs_with_no_snps_by_sample, nrow)) > 0]
   
   gene_depth_info <- rbindlist(lapply(names(depth_info), function(x){
      
      tmp <- gene_assembly[contig %in% depth_info[[x]]$V1, list(OC, gene, contig, length, expN, expS)]
      tmp$source <- x
      return(tmp)
      
   }))
   
   contig_level_microdiv <- cleaned_data[,list(sum(pi, na.rm=T), length(contig_pos), gm_mean(dep)), by=key(cleaned_data)]
   colnames(contig_level_microdiv)[3:5] = c("pi", "num_snps", "harmonic_number")
   
   #can't calc. taj. D if sample depth < 2, and SHOULDN'T unless it's 3.
   contig_level_microdiv <- contig_level_microdiv[harmonic_number > 2,]
   
   contig_level_microdiv$contig_len <- assembled_contigs$length[match(contig_level_microdiv$contig, assembled_contigs$contig)]
   contig_level_microdiv$theta <- (contig_level_microdiv$num_snps/harmonic[contig_level_microdiv$harmonic_number])/(contig_level_microdiv$contig_len-(total_snps$count[match(contig_level_microdiv$contig, total_snps$contig)]-contig_level_microdiv$num_snps))
   
   #Modify pi to follow SNP cpount corrections exactly like theta
   contig_level_microdiv$pi <- contig_level_microdiv$pi/(contig_level_microdiv$contig_len-(total_snps$count[match(contig_level_microdiv$contig, total_snps$contig)]-contig_level_microdiv$num_snps))
   
   contig_level_microdiv <- contig_level_microdiv[, -5]
   
   #If linked SNPs was done, corrects the microdiv_data for the directly observed read results
   if(dir.exists("MetaPop/09.Linked_SNPs")){
     print("Correcting subsamples for linked snps")
      
     mined_reads <- fread("MetaPop/09.Linked_SNPs/linked_snp_results.tsv", sep = "\t")
     mined_reads[, three_snps := NULL]
     
     mined_reads <- mined_reads[(mined_reads$ref_first+mined_reads$ref_count+mined_reads$snp_count+mined_reads$ref_second)>0,]
     
     mined_reads$original_codon <- mapply(function(codons, index){
        return(codons[index])
     },
     split_codons[match(mined_reads$contig_gene, gene_assembly$gene)], mined_reads$codon)
     
     mined_reads <- mined_reads[mined_reads$original_codon %in% construct,]
     
     has_alts <- grepl("[*]", mined_reads$refs)
     
     #We need to get the partial alts here for N/S
     mined_reads$alt1 <- mapply(function(read, alt, correct, is_alt){
        
        #We neither need to nor should find alts for triples; they are taken care of as part of their 2-part subsets.
        if(!is_alt){
           return("")
        }
        
        str <- strsplit(read, split="")[[1]]
        snp <- strsplit(alt, split="")[[1]]
        cor <- strsplit(correct, split="")[[1]]
        
        pos <- which(str == "*")
           
           if(pos == 1){
              hybrid <- c(cor[pos], str[pos+1], snp[pos+2])
           }
        
           if(pos == 2){
              hybrid <- c(str[pos-1], cor[pos], snp[pos+1])
           }
        
           if(pos == 3){
              hybrid <- c(str[pos-2], snp[pos-1], cor[pos])
           }
           
           hybrid <- paste0(hybrid, collapse = "")
        
        
        return(hybrid)
        
     }, mined_reads$refs, mined_reads$snp, mined_reads$original_codon, has_alts)
     mined_reads$alt2 <- mapply(function(read, alt, correct, is_alt){
        
        #We neither need to nor should find alts for triples; they are taken care of as part of their 2-part subsets.
        if(!is_alt){
           return("")
        }
        
        str <- strsplit(read, split="")[[1]]
        snp <- strsplit(alt, split="")[[1]]
        cor <- strsplit(correct, split="")[[1]]
        
        pos <- which(str == "*")
        
        if(pos == 1){
           hybrid <- c(cor[pos], snp[pos+1], str[pos+2])
        }
        
        if(pos == 2){
           hybrid <- c(snp[pos-1], cor[pos], str[pos+1])
        }
        
        if(pos == 3){
           hybrid <- c(snp[pos-2], str[pos-1], cor[pos])
        }
        
        hybrid <- paste0(hybrid, collapse = "")
        
        
        return(hybrid)
        
     }, mined_reads$refs, mined_reads$snp, mined_reads$original_codon, has_alts)
     
     mined_reads$refs <- mapply(function(read, correct){
        
        str <- strsplit(read, split="")[[1]]
        if(any(str == "*")){
           pos <- which(str == "*")
           cor <- strsplit(correct, split="")[[1]]
           str[pos] <- cor[pos]
           read <- paste0(str, collapse = "")
        }
        return(read)
        
     }, mined_reads$refs, mined_reads$original_codon)
     #This is the codon with all SNPs present - i.e. 2 or 3 SNPs on the same codon simultaneously
     mined_reads$snp <- mapply(function(read, correct){
        
        str <- strsplit(read, split="")[[1]]
        if(any(str == "*")){
           pos <- which(str == "*")
           cor <- strsplit(correct, split="")[[1]]
           str[pos] <- cor[pos]
           read <- paste0(str, collapse = "")
        }
        return(read)
        
     }, mined_reads$snp, mined_reads$original_codon)
     
     setkeyv(mined_reads, c("source", "contig_gene", "codon"))
     mined_reads[, overall_prop := sum(ref_count+snp_count+ref_first+ref_second), by = key(mined_reads)]
     mined_reads$overall_prop <- sub_samp/mined_reads$overall_prop
     
     #If there were invalid codons, e.g. TAN, GNN, NNN, etc. where there is a missing base, mistakes may occur. This cleans them out.
     mined_reads <- mined_reads[nchar(alt1)>0 & nchar(alt2)>0 & nchar(refs)>0 & nchar(snp)>0,]
     
     mined_reads[, ref_count := ifelse(overall_prop < 1, round(ref_count*overall_prop), ref_count) , by = .I]
     mined_reads[, snp_count := ifelse(overall_prop < 1, round(snp_count*overall_prop), snp_count) , by = .I]
     mined_reads[, ref_first := ifelse(overall_prop < 1, round(ref_first*overall_prop), ref_first) , by = .I]
     mined_reads[, ref_second := ifelse(overall_prop < 1, round(ref_second*overall_prop), ref_second) , by = .I]
     
     #original codon doesn't contribute
     mined_reads$ref_count <- 0

     mined_reads$orignal_AA <- proteinConstruct[match(mined_reads$original_codon, construct)]
     
     mined_reads$ref_AA <- proteinConstruct[match(mined_reads$refs, construct)]
     mined_reads$snp_AA <- proteinConstruct[match(mined_reads$snp, construct)]
     mined_reads$alt1_AA <- proteinConstruct[match(mined_reads$alt1, construct)]
     mined_reads$alt2_AA <- proteinConstruct[match(mined_reads$alt2, construct)]
     
     mined_reads$obsN <- 0
     mined_reads$obsN <- ifelse(mined_reads$snp_AA==mined_reads$ref_AA, mined_reads$obsN, mined_reads$obsN+mined_reads$snp_count)
     mined_reads$obsN <- ifelse(mined_reads$alt1_AA==mined_reads$ref_AA, mined_reads$obsN, mined_reads$obsN+mined_reads$ref_first)
     mined_reads$obsN <- ifelse(mined_reads$alt2_AA==mined_reads$ref_AA, mined_reads$obsN, mined_reads$obsN+mined_reads$ref_second)
     
     mined_reads$obsS <- 0
     mined_reads$obsS <- ifelse(mined_reads$snp_AA!=mined_reads$ref_AA, mined_reads$obsS, mined_reads$obsS+mined_reads$snp_count)
     mined_reads$obsS <- ifelse(mined_reads$alt1_AA!=mined_reads$ref_AA, mined_reads$obsS, mined_reads$obsS+mined_reads$ref_first)
     mined_reads$obsS <- ifelse(mined_reads$alt2_AA!=mined_reads$ref_AA, mined_reads$obsS, mined_reads$obsS+mined_reads$ref_second)
     
     mined_reads <- mined_reads[, list(obsN = sum(obsN, na.rm=T), obsS = sum(obsS, na.rm=T), ref = unique(refs)), by = key(mined_reads)]
     
     mined_reads$obsN <- ifelse(mined_reads$obsN+mined_reads$obsS>sub_samp, round((mined_reads$obsN/(mined_reads$obsN+mined_reads$obsS))*sub_samp), mined_reads$obsN)
     mined_reads$obsS <- ifelse(mined_reads$obsN+mined_reads$obsS>sub_samp, round((mined_reads$obsS/(mined_reads$obsN+mined_reads$obsS))*sub_samp), mined_reads$obsS)
     
     microdiv_data$match_key <- paste0(microdiv_data$source, microdiv_data$contig_gene, microdiv_data$codon)
     mined_reads$match_key <- paste0(mined_reads$source, mined_reads$contig_gene, mined_reads$codon)
     
     mined_reads <- mined_reads[mined_reads$match_key %in% microdiv_data$match_key,]
     
     selector <- microdiv_data$match_key %in% mined_reads$match_key
     
     #sets counts on each codon for microdiv data to zero to prevent double counting; sets codons to blank for correcting
     
     microdiv_data[, obsN := ifelse(selector, 0, obsN)]
     microdiv_data[, obsS := ifelse(selector, 0, obsS)]
     
     microdiv_data[, affected_by_linked_snps := F]
     microdiv_data[, affected_by_linked_snps := ifelse(selector, T, F)]

     selector <- match(mined_reads$match_key, microdiv_data$match_key)
     
     #sets only the first row for each codon with linked data to have the linked data values
     microdiv_data$obsN[selector] <- mined_reads$obsN
     microdiv_data$obsS[selector] <- mined_reads$obsS
     #microdiv_data$original_codon[selector] <- mined_reads$ref
     
     microdiv_data[,match_key:=NULL]
     
   }else{
      print("Linked SNPs corrections not performed; pNpS results may NOT be fully accurate.")
   }
   
   setkeyv(microdiv_data, c("contig_gene", "source"))
   pnps <- microdiv_data[, list(sum(obsN), sum(obsS)), by = key(microdiv_data)]
   colnames(pnps)[3:4] = c("obsN", "obsS") 
   
   gene_level_microdiv$expN <- gene_assembly$expN[match(gene_level_microdiv$contig_gene, gene_assembly$gene)]
   gene_level_microdiv$expS <- gene_assembly$expS[match(gene_level_microdiv$contig_gene, gene_assembly$gene)]
   
   pnps[, match_key := paste0(contig_gene, source),]
   gene_level_microdiv[, key := paste0(contig_gene, source)]
   
   gene_level_microdiv$obsN <- pnps$obsN[match(gene_level_microdiv$key, pnps$match_key)]
   gene_level_microdiv$obsS <- pnps$obsS[match(gene_level_microdiv$key, pnps$match_key)]
   
   #Calculating pNpS ratio
   gene_level_microdiv$pNpS_ratio <- (gene_level_microdiv$obsN/gene_level_microdiv$obsS)/(gene_level_microdiv$expN/gene_level_microdiv$expS)
   
   contig_level_microdiv[,snps_present := T]
   gene_level_microdiv[,snps_present := T]
   
   #Fix up the contigs w/o snps to match the form of the contigs with snps
   contigs_with_no_snps_by_sample <- rbindlist(contigs_with_no_snps_by_sample)
   contigs_with_no_snps_by_sample <- contigs_with_no_snps_by_sample[,list(V1, source, V2)]
   colnames(contigs_with_no_snps_by_sample)[c(1,3)] = c("contig", "contig_len")
   contigs_with_no_snps_by_sample[, pi := 0]
   contigs_with_no_snps_by_sample[, num_snps := 0]
   contigs_with_no_snps_by_sample[, theta := 0]
   contigs_with_no_snps_by_sample[, snps_present := FALSE]
   contigs_with_no_snps_by_sample <- contigs_with_no_snps_by_sample[, list(contig, source, pi, num_snps, contig_len, theta, snps_present)]
   
   #Bind the samples together.
   contig_level_microdiv <- rbindlist(list(contig_level_microdiv, contigs_with_no_snps_by_sample))
   
   gene_depth_info[, key := paste0(gene, source)]
   gene_depth_info <- gene_depth_info[!(key %in% gene_level_microdiv$key),]
   
   gene_depth_info[, key := NULL]
   gene_level_microdiv[, key := NULL]
   
   #fix up the missing genes to match the form of the genes with snps
   gene_depth_info[, pi := 0]
   gene_depth_info[, num_snps := 0]
   gene_depth_info[, obsN := 0]
   gene_depth_info[, obsS := 0]
   gene_depth_info[, theta := 0]
   gene_depth_info[, taj_D := 0]
   gene_depth_info[, pNpS_ratio := 0]
   gene_depth_info[, snps_present := FALSE]
   
   gene_depth_info <- gene_depth_info[, list(gene, source, contig, pi, num_snps, length, theta, taj_D, expN, expS, obsN, obsS, pNpS_ratio, snps_present)]
   colnames(gene_depth_info) = colnames(gene_level_microdiv)
   
   gene_level_microdiv <- rbindlist(list(gene_level_microdiv, gene_depth_info))
   

   
   fwrite(microdiv_data, "MetaPop/10.Microdiversity/global_raw_microdiversity_data_snp_loci_only.tsv", sep = "\t")
   
   fwrite(gene_level_microdiv, "MetaPop/10.Microdiversity/global_gene_microdiversity.tsv", sep = "\t")
   fwrite(contig_level_microdiv, "MetaPop/10.Microdiversity/global_contig_microdiversity.tsv", sep = "\t")
   
   fwrite(codon_pos_count, "MetaPop/10.Microdiversity/global_codon_position_summary.tsv", sep = "\t")
   
   }
   
   if(sum(codon_pos_count$third_pos) < sum(codon_pos_count$second_pos) | sum(codon_pos_count$third_pos) < sum(codon_pos_count$first_pos)){
      print("More SNPs were observed in either the first or second positions of codons than in the third position. This is suspicious, and results should be treated cautiously.")
      flush.console()
   }
   
   
   #Local level stuff
   vcfs <- list.files(path = "MetaPop/05.Variant_Calls", full.names = T)
   
   cl <- makeCluster(min(threads, detectCores(), length(vcfs)))
   registerDoParallel(cl)
   clusterExport(cl, c("library_location", "vcfs"))
   clusterEvalQ(cl=cl, expr = library(data.table, lib.loc = library_location))

   # Read and process mpileup files in parallel
   vcf_sums <- foreach(i=vcfs) %dopar%{
      tmp <- fread(i, sep = "\t", col.names = c("contig", "pos", "ref", "alt"))
      tmp[, source := substr(i, 26, nchar(i)-17)]
      tmp[, key := paste0(contig,"_", pos, source)]
      
      return(tmp)
   }
   
   stopCluster(cl)
   
   vcf_sums <- rbindlist(vcf_sums)
   
   
   microdiv_data[, key := paste0(contig_pos, source)]
   
   microdiv_data_local <- microdiv_data[key %in% vcf_sums$key,]
   
   {
      setkey(microdiv_data_local, "contig_gene")
      total_snps_genic <- microdiv_data_local[, list(length(unique(contig_pos)), unique(contig)), by = key(microdiv_data_local)]
      colnames(total_snps_genic)[2:3] = c("count", "contig")
      
      setkeyv(microdiv_data_local, c("contig_gene", "source"))
      gene_level_microdiv <- microdiv_data_local[, list(unique(contig), sum(pi, na.rm = T), length(contig), gm_mean(sub_samp_depth)), by=key(microdiv_data_local)]
      colnames(gene_level_microdiv)[3:6] = c("contig", "pi", "num_snps", "harmonic_number")
      
      #Shouldn't calculate microdiv if the depth is too low
      gene_level_microdiv <- gene_level_microdiv[harmonic_number > 2,]
      
      gene_level_microdiv$gene_len <- gene_assembly$length[match(gene_level_microdiv$contig_gene, gene_assembly$gene)]
      gene_level_microdiv$theta <- (gene_level_microdiv$num_snps/harmonic[gene_level_microdiv$harmonic_number])/(gene_level_microdiv$gene_len-(total_snps_genic$count[match(gene_level_microdiv$contig_gene, total_snps_genic$contig_gene)]-gene_level_microdiv$num_snps))
      
      #Modify pi to follow SNP count corrections exactly like theta
      gene_level_microdiv$pi <- gene_level_microdiv$pi/(gene_level_microdiv$gene_len-(total_snps_genic$count[match(gene_level_microdiv$contig_gene, total_snps_genic$contig_gene)]-gene_level_microdiv$num_snps))
      
      gene_level_microdiv$taj_D <- (gene_level_microdiv$pi-(gene_level_microdiv$num_snps/harmonic[gene_level_microdiv$harmonic_number]))/sqrt((e1[gene_level_microdiv$harmonic_number]*gene_level_microdiv$num_snps)+(e2[gene_level_microdiv$harmonic_number]*gene_level_microdiv$num_snps*(gene_level_microdiv$num_snps-1)))
      
      
      #check it out to make sure this works
      gene_level_microdiv <- gene_level_microdiv[, -6]
      
      saved_total_snps_genic <- total_snps_genic
      
      # Repeated at the level of the codon; all of the data gets incorporated as well as possible
      total_snps_genic$contig_gene = NULL
      
      setkey(non_genic_data, "contig")
      total_snps <- non_genic_data[, length(unique(contig_pos)), by=key(non_genic_data)]
      names(total_snps)[2] = "count"
      
      total_snps <- rbind(total_snps_genic, total_snps)
      
      setkey(total_snps, "contig")
      total_snps <- total_snps[, sum(count), by=key(total_snps)]
      colnames(total_snps)[2] = "count"
      
      cleaned_data <- data.table(contig_pos = microdiv_data_local$contig_pos, contig = microdiv_data_local$contig, pi = microdiv_data_local$pi, source = microdiv_data_local$source, dep = microdiv_data_local$sub_samp_depth)
      cleaned_data <- rbind(cleaned_data, data.table(contig_pos = non_genic_data$contig_pos, contig = non_genic_data$contig, pi=non_genic_data$pi, source = non_genic_data$source, dep = non_genic_data$sub_samp_depth))
      
      setkeyv(cleaned_data, c("contig", "source"))
      
      extant_names <- unique(cleaned_data$source)
      contigs_with_real_observations_by_sc <- cleaned_data[,list(list(unique(contig))), by = source]$V1
      names(contigs_with_real_observations_by_sc) = extant_names
      
      #Finds the contigs that passed the cov/depth filter but had no passing SNPs. 
      #This is added to the contig_level_microdiv as the has_snps = FALSE section.
      contigs_with_no_snps_by_sample <- lapply(names(depth_info), function(x){
         
         #if a sample has passing contigs, but no SNPs, we want to return all passing contigs in the depth file
         if(is.null(contigs_with_real_observations_by_sc[[x]])){
            return(depth_info[[x]])
         }else{
            #Else return all of the contigs that aren't in the data already - i.e. the SNPless ones.
            return(depth_info[[x]][!(V1 %in% contigs_with_real_observations_by_sc[[x]])])
         }
         
      })
      names(contigs_with_no_snps_by_sample) = names(depth_info)
      
      #Empty ones are fully covered already
      contigs_with_no_snps_by_sample <- contigs_with_no_snps_by_sample[unlist(lapply(contigs_with_no_snps_by_sample, nrow)) > 0]
      
      gene_depth_info <- rbindlist(lapply(names(depth_info), function(x){
         
         tmp <- gene_assembly[contig %in% depth_info[[x]]$V1, list(OC, gene, contig, length, expN, expS)]
         tmp$source <- x
         return(tmp)
         
      }))
      
      contig_level_microdiv <- cleaned_data[,list(sum(pi, na.rm=T), length(contig_pos), gm_mean(dep)), by=key(cleaned_data)]
      colnames(contig_level_microdiv)[3:5] = c("pi", "num_snps", "harmonic_number")
      
      #can't calc. taj. D if sample depth < 2, and SHOULDN'T unless it's 3.
      contig_level_microdiv <- contig_level_microdiv[harmonic_number > 2,]
      
      contig_level_microdiv$contig_len <- assembled_contigs$length[match(contig_level_microdiv$contig, assembled_contigs$contig)]
      contig_level_microdiv$theta <- (contig_level_microdiv$num_snps/harmonic[contig_level_microdiv$harmonic_number])/(contig_level_microdiv$contig_len-(total_snps$count[match(contig_level_microdiv$contig, total_snps$contig)]-contig_level_microdiv$num_snps))
      
      #Modify pi to follow SNP cpount corrections exactly like theta
      contig_level_microdiv$pi <- contig_level_microdiv$pi/(contig_level_microdiv$contig_len-(total_snps$count[match(contig_level_microdiv$contig, total_snps$contig)]-contig_level_microdiv$num_snps))
      
      contig_level_microdiv <- contig_level_microdiv[, -5]
      
      #Fixed up to linked. Fantastic, time to do THIS...
      #If linked SNPs was done, corrects the microdiv_data_local for the directly observed read results
      if(dir.exists("MetaPop/09.Linked_SNPs")){
         print("Correcting local subsamples for linked snps")

         microdiv_data_local$match_key <- paste0(microdiv_data_local$source, microdiv_data_local$contig_gene, microdiv_data_local$codon)
         mined_reads$match_key <- paste0(mined_reads$source, mined_reads$contig_gene, mined_reads$codon)
         
         mined_reads <- mined_reads[mined_reads$match_key %in% microdiv_data_local$match_key,]
         
         selector <- microdiv_data_local$match_key %in% mined_reads$match_key
         
         #sets counts on each codon for microdiv data to zero to prevent double counting; sets codons to blank for correcting
         
         microdiv_data_local[, obsN := ifelse(selector, 0, obsN)]
         microdiv_data_local[, obsS := ifelse(selector, 0, obsS)]
         
         microdiv_data_local[, affected_by_linked_snps := F]
         
         selector <- match(mined_reads$match_key, microdiv_data_local$match_key)
         
         #sets only the first row for each codon with linked data to have the linked data values
         microdiv_data_local$obsN[selector] <- mined_reads$obsN
         microdiv_data_local$obsS[selector] <- mined_reads$obsS
         #microdiv_data_local$original_codon[selector] <- mined_reads$ref
         
         microdiv_data_local[,match_key:=NULL]
         
      }else{
         print("Linked SNPs corrections not performed; pNpS results may NOT be fully accurate.")
         print("We recommend running MetaPop_Mine_Reads.R before the microdiversity module to correct SNP data over reads.")
      }
      
      setkeyv(microdiv_data_local, c("contig_gene", "source"))
      pnps <- microdiv_data_local[, list(sum(obsN), sum(obsS)), by = key(microdiv_data_local)]
      colnames(pnps)[3:4] = c("obsN", "obsS") 
      
      gene_level_microdiv$expN <- gene_assembly$expN[match(gene_level_microdiv$contig_gene, gene_assembly$gene)]
      gene_level_microdiv$expS <- gene_assembly$expS[match(gene_level_microdiv$contig_gene, gene_assembly$gene)]
      
      pnps[, match_key := paste0(contig_gene, source),]
      gene_level_microdiv[, key := paste0(contig_gene, source)]
      
      gene_level_microdiv$obsN <- pnps$obsN[match(gene_level_microdiv$key, pnps$match_key)]
      gene_level_microdiv$obsS <- pnps$obsS[match(gene_level_microdiv$key, pnps$match_key)]
      
      #Calculating pNpS ratio
      gene_level_microdiv$pNpS_ratio <- (gene_level_microdiv$obsN/gene_level_microdiv$obsS)/(gene_level_microdiv$expN/gene_level_microdiv$expS)
      
      contig_level_microdiv[,snps_present := T]
      gene_level_microdiv[,snps_present := T]
      
      #Fix up the contigs w/o snps to match the form of the contigs with snps
      contigs_with_no_snps_by_sample <- rbindlist(contigs_with_no_snps_by_sample)
      contigs_with_no_snps_by_sample <- contigs_with_no_snps_by_sample[,list(V1, source, V2)]
      colnames(contigs_with_no_snps_by_sample)[c(1,3)] = c("contig", "contig_len")
      contigs_with_no_snps_by_sample[, pi := 0]
      contigs_with_no_snps_by_sample[, num_snps := 0]
      contigs_with_no_snps_by_sample[, theta := 0]
      contigs_with_no_snps_by_sample[, snps_present := FALSE]
      contigs_with_no_snps_by_sample <- contigs_with_no_snps_by_sample[, list(contig, source, pi, num_snps, contig_len, theta, snps_present)]
      
      #Bind the samples together.
      contig_level_microdiv <- rbindlist(list(contig_level_microdiv, contigs_with_no_snps_by_sample))
      
      gene_depth_info[, key := paste0(gene, source)]
      gene_depth_info <- gene_depth_info[!(key %in% gene_level_microdiv$key),]
      
      gene_depth_info[, key := NULL]
      gene_level_microdiv[, key := NULL]
      
      #fix up the missing genes to match the form of the genes with snps
      gene_depth_info[, pi := 0]
      gene_depth_info[, num_snps := 0]
      gene_depth_info[, obsN := 0]
      gene_depth_info[, obsS := 0]
      gene_depth_info[, theta := 0]
      gene_depth_info[, taj_D := 0]
      gene_depth_info[, pNpS_ratio := 0]
      gene_depth_info[, snps_present := FALSE]
      
      gene_depth_info <- gene_depth_info[, list(gene, source, contig, pi, num_snps, length, theta, taj_D, expN, expS, obsN, obsS, pNpS_ratio, snps_present)]
      colnames(gene_depth_info) = colnames(gene_level_microdiv)
      
      gene_level_microdiv <- rbindlist(list(gene_level_microdiv, gene_depth_info))
      
      fwrite(microdiv_data_local, "MetaPop/10.Microdiversity/local_raw_microdiversity_data_snp_loci_only.tsv", sep = "\t")
      
      fwrite(gene_level_microdiv, "MetaPop/10.Microdiversity/local_gene_microdiversity.tsv", sep = "\t")
      fwrite(contig_level_microdiv, "MetaPop/10.Microdiversity/local_contig_microdiversity.tsv", sep = "\t")
      
      fwrite(codon_pos_count, "MetaPop/10.Microdiversity/local_codon_position_summary.tsv", sep = "\t")
      
   }
   
}

#FST

if(length(vcfs) > 1){
   print("Calculating fixation index (Fst)")
   
   #FST is a relevant stat only when a contig is observed in at least two samples. Thus, only contigs seen 2+ times are used
   setkey(microdiv_data, "contig")
   contig_counts_by_sample <- microdiv_data[, list(length(unique(source)), length(unique(pos))), by=key(microdiv_data)]
   colnames(contig_counts_by_sample)[2:3] = c("num_occur", "total_snps")
   
   contig_counts_by_sample <- contig_counts_by_sample[contig_counts_by_sample$num_occur > 1,]
   
   #Remove the singletons
   microdiv_data <- microdiv_data[microdiv_data$contig %in% contig_counts_by_sample$contig,]
   
   all_samples <- unique(microdiv_data$source)
   
   samp_pairs <- rbindlist(lapply(1:length(all_samples), function(x){
      
   CJ(all_samples[x], all_samples[x:length(all_samples)])
   
   }))
   
   colnames(samp_pairs) = c("row_sample", "col_sample")
   
   samp_pairs <- samp_pairs[row_sample != col_sample,]
   
   contigs <- contig_counts_by_sample$contig
   
   #The repository for FST data. Contains each pairwise comparison of samples for each contig, including pairs where a contig was not observed in one or both samples.
   FST_pairs <- data.table(row_samp = rep(samp_pairs$row_sample, length(contigs)), col_samp = rep(samp_pairs$col_sample, length(contigs)), contig = rep(contigs, each=nrow(samp_pairs)), fst = NA)
   
   #Presplit dataframes
   namesave = unique(microdiv_data$contig)
   microdiv_data <- microdiv_data[, list(list(.SD)), by = contig]$V1
   names(microdiv_data) = namesave
   
   #Pre split internal data frames
   microdiv_data <- lapply(microdiv_data, function(x){
      namesave = unique(x$source)
      x <- x[, list(list(.SD)), by = source]$V1
      names(x) = namesave
      return(x)
   })
   
   cl <- makeCluster(min(threads, detectCores()))
   registerDoParallel(cl)
   
   fst_results <- unlist(foreach(x = contigs) %dopar% {
      sub <- microdiv_data[[x]]
      
      contig_len <- contig_level_microdiv$contig_len[match(x, contig_level_microdiv$contig)]
      tot_snp_ct <- contig_counts_by_sample$total_snps[match(x, contig_counts_by_sample$contig)]
      
      
      unlist(lapply(1:(length(sub)-1), function(i){
         
         t1 <- sub[[i]]
         cov1 <- contig_len-(tot_snp_ct-nrow(t1))
         pi1 <- sum(t1$pi, na.rm = T)/cov1
         
         unlist(lapply((i+1):length(sub), function(j){
            
            t2 <- sub[[j]]
            
            cov2 <- contig_len-(tot_snp_ct-nrow(t2))
            pi2 <- sum(t2$pi, na.rm = T)/cov2
            
            t1 <- t1[t1$pos %in% t2$pos,]
            t2 <- t2[t2$pos %in% t1$pos,]
            
            cov3 <- contig_len-(tot_snp_ct-nrow(t1))
            
            if(nrow(t1)>0 & nrow(t2)>0){
               
               #This should produce pi between
               t1<- t1[order(t1$pos),]
               t2 <- t2[order(t2$pos),]
               
               cc <- t1$sub_samp_depth*t2$sub_samp_depth
               
               at <- (t1$sub_samp_a*t2$sub_samp_t)/cc
               ac <- (t1$sub_samp_a*t2$sub_samp_c)/cc
               ag <- (t1$sub_samp_a*t2$sub_samp_g)/cc
               ta <- (t1$sub_samp_t*t2$sub_samp_a)/cc
               tc <- (t1$sub_samp_t*t2$sub_samp_c)/cc
               tg <- (t1$sub_samp_t*t2$sub_samp_g)/cc
               ca <- (t1$sub_samp_c*t2$sub_samp_a)/cc
               ct <- (t1$sub_samp_c*t2$sub_samp_t)/cc
               cg <- (t1$sub_samp_c*t2$sub_samp_g)/cc
               ga <- (t1$sub_samp_g*t2$sub_samp_a)/cc
               gt <- (t1$sub_samp_g*t2$sub_samp_t)/cc
               gc <- (t1$sub_samp_g*t2$sub_samp_c)/cc
               
               fst <- sum(at,ac,ag,ta,tc,tg,ca,ct,cg,ga,gt,gc)/cov3
               fst <- 1-(((pi1+pi2)/2)/fst)
               
               if(is.nan(fst)){
                  fst <- 1
               }
               
               if(fst < 0){
                  fst <- 0
               }
               
               
               
            } else {
               
               fst = 1
               
            }
            
            res <- list(names(sub)[i], names(sub)[j], x, fst)
            
            return(res)
            
         }))
         
         
         
      }))
      
   })
   
   stopCluster(cl)
   
   
   fst_results <- data.table(row_samp = fst_results[c(T,F,F,F)], col_samp = fst_results[c(F,T,F,F)], contig=fst_results[c(F,F,T,F)], fst = as.numeric(fst_results[c(F,F,F,T)]))
   
   fst_results <- rbind(fst_results, FST_pairs)
   fst_results <- unique(fst_results, by=key(fst_results))
   fst_results <- fst_results[order(contig, row_samp, col_samp),]
   
   fwrite(fst_results, "MetaPop/10.Microdiversity/fixation_index.tsv", sep = "\t")
   
}else{
   print("Fixation index cannot be calculated for only one sample.")
}
















