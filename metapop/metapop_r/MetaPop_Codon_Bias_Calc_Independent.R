options <- commandArgs(trailingOnly = T)

dir = options[1]

gene_file <- options[2]

if(length(options)>2){
library_location <- options[3]
}

if(library_location == ""){
  library_location <- .libPaths()
}

suppressMessages(suppressWarnings(library(data.table, lib.loc = library_location)))

genes = readLines(gene_file)

setwd(dir)

gene_locs = which(substr(genes, 1, 1) == ">")
gene_headers = genes[gene_locs]
genes[gene_locs] = ""
gene_headers = unlist(lapply(gene_headers, function(x){
  name = strsplit(x, " # ")[[1]][1]
  name = substr(name, 2, nchar(name))
  return(name)
}))

sizes = (gene_locs[2:length(gene_locs)] - gene_locs[1:(length(gene_locs)-1)])
#Append final index
sizes = c(sizes, length(genes)-gene_locs[length(gene_locs)]+1)

splits = rep(1:length(gene_headers), times = sizes)

gsplit = split(genes, splits)

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

codon_data = data.table(matrix(unlist(lapply(gsplit, (function(x){
  #Combine the seq into a single string
  seq = paste(x, collapse = "")
  #split to chars
  seq = unlist(strsplit(seq, split = ""))
  #recombine to codons
  seq = paste0(seq[c(T,F,F)], seq[c(F,T,F)], seq[c(F,F,T)])
  #count the occurrences of each codon
  seq = table(seq)
  
  #create a named, ordered list of 64 counts exactly matching the order of construct
  counts = rep(0, 64)
  names(counts) = construct
  
  #fill the appropriate, ordred counts with the real data, leaving zeroes elsewhere.
  counts[names(seq)] = seq
  
  #print(length(counts))
  
  #return the ordered count data
  return(counts)
  
}))), ncol = 64))


genes = data.table(V1 = gene_headers)
rm(gene_locs, sizes, splits, gsplit, gene_headers)
genes[, V2 := gsub('_\\d+$', '', genes$V1), ]


#genes <- fread("MetaPop/08.Codon_Bias/codon_bias_sequences.tsv", sep = "\t", header = F, select = c(3,4))
colnames(genes) = c("V1", "V2")

#codon_data <- fread("MetaPop/08.Codon_Bias/codon_counts_by_gene.tsv", sep = "\t", header = F)
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

if(!dir.exists("MetaPop")){
  dir.create("MetaPop")
}
if(!dir.exists("MetaPop/08.Codon_Bias")){
  dir.create("MetaPop/08.Codon_Bias")
}

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