options <- commandArgs(trailingOnly = T)

check_for_default <- integer(0)

directory_name <- which(grepl("-dir", options))

quitNow <- 0

first_sample <- which(grepl("-first", options))

snp_obs_cutoff <- which(grepl("-obs", options))

min_rep <- which(grepl("-rep", options))

var_quality <- which(grepl("-var_qual", options))

if(identical(snp_obs_cutoff, check_for_default)){
  print("Defaulting number of observations of a SNV to be called a SNP to 4")
  snp_obs_cutoff <- 4
}else{
  snp_obs_cutoff <- as.numeric(options[snp_obs_cutoff + 1])
}

if(identical(min_rep, check_for_default)){
  print("Defaulting minimum pop. representation for a SNV to be called a SNP to 1%")
  min_rep <- 0.01
}else{
  min_rep <- as.numeric(options[min_rep + 1])/100
}

if(identical(var_quality, check_for_default)){
  print("Defaulting vcf min phred call quality to 20")
  var_quality <- 20
}else{
  var_quality <- as.numeric(options[var_quality + 1])
}

time_series <- F

if(identical(first_sample, check_for_default)){
  print("No sample set as the point of reference for nucleotide calls. This is an optional setting used for time series.")
  print("Nucleotide defined as reference will be the most abundant nucleotide at each SNV locus across all samples.")
  first_sample <- ""
}else{
  first_sample <- as.character(options[first_sample+1])
  print(paste(first_sample, "will be used as the point of reference for nucleotide bases in SNV calls."))
  print("This should be the name metapop uses for the sample, internally. Check metapop_run_settings/SNV_call_base_names.tsv if this step fails.")
  print("Caution: Only positions found within the sample of reference will be retained going forward.")
  time_series <- T
}

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

run_parameters <- read.csv("metapop_run_settings/run_settings.tsv", sep = "\t", stringsAsFactors = F)

library_location <- as.character(run_parameters$setting[run_parameters$parameter=="Library Location"])

suppressMessages(suppressWarnings(library(doParallel, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(data.table, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(stringr, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(Biostrings, lib.loc = library_location)))

run_summaries <- fread("metapop_run_settings/run_dates.tsv", sep = "\t")

run_summaries$begin[run_summaries$stage=="SNP Calling"] <- as.character(Sys.time())

fwrite(run_summaries, "metapop_run_settings/run_dates.tsv", sep = "\t")

if(time_series){
  if(identical(which(run_parameters$parameter == "First Sample"), check_for_default)){
    run_parameters <- rbind(run_parameters, data.frame(parameter = "First Sample", setting = first_sample))
  }else{
    run_parameters$setting[which(run_parameters$parameter == "First Sample")] <- first_sample
  }
  fwrite(run_parameters, "metapop_run_settings/run_settings.tsv", sep = "\t")
}

run_parameters <- rbind(run_parameters, data.frame(parameter = c("Variant Base Call Cutoff", "Min SNP Obs.", "Min Pop. Rep."), setting = c(var_quality, snp_obs_cutoff, as.character(min_rep))))
fwrite(run_parameters, "metapop_run_settings/run_settings.tsv", sep = "\t")

threads <- as.numeric(run_parameters$setting[run_parameters$parameter=="Threads"])

samtools_location <- as.character(run_parameters$setting[run_parameters$parameter=="Samtools Location"])
bcftools_location <- as.character(run_parameters$setting[run_parameters$parameter=="BCFTools Location"])

original_assembly <- as.character(run_parameters$setting[run_parameters$parameter=="Assembly"])
geneFile <- as.character(run_parameters$setting[run_parameters$parameter=="Genes"])

print(run_parameters)
if(time_series){
  print(paste("First sample:", first_sample))
}
  
print("Making initial SNV calls... ")

#Reacquire reads with actual info
newAlignments <- list.files(pattern = ".bam", path = "metapop_preprocessed_reads", full.names = T)
newAlignments <- newAlignments[!grepl(".bai", newAlignments)]

#If any were lost, select the base names that lived
base_names <- substr(newAlignments, 28, nchar(newAlignments)-17)

if(time_series){
  fwrite(data.table(original_name = newAlignments, first_sample_name = base_names), "metapop_run_settings/SNV_call_base_names.tsv", sep = "\t")  
}


if(!dir.exists("metapop_variants")){
  system("mkdir metapop_variants")
}
if(!dir.exists("metapop_pileups")){
  system("mkdir metapop_pileups")
}
if(!dir.exists("metapop_temp")){
  system("mkdir metapop_temp")
}

plf <- c(a= "*", b = "*", c = "*", d = "*", e = 1)
writeLines(plf, "metapop_temp/ploidy_file.tsv", sep = "\t")

ploidy_file <- "metapop_temp/ploidy_file.tsv"

cl <- makeCluster(min(threads, detectCores(), length(newAlignments)))
registerDoParallel(cl)

#Calls SNPs with no indels and quality filters them, reads in the results
silent <- foreach(i=1:length(newAlignments)) %dopar%{
  #Produces SNV calls, filters these for SNV positions only (remove indels)  
  system(paste(bcftools_location, "bcftools mpileup -Ob -I -f ", original_assembly, " ", newAlignments[i], " | ", bcftools_location, "bcftools call -mv -v -Ob --ploidy-file ", ploidy_file, " | ", 
               bcftools_location, "bcftools filter -i'QUAL>", var_quality, "' | ", 
               bcftools_location, "bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' > metapop_variants/", base_names[i], "_variants.tsv", sep =""))
  
  system(paste("cut -f1,2 metapop_variants/", base_names[i], "_variants.tsv | sed 's/\\t/_/' > metapop_temp/", base_names[i], "_contig_pos.tsv", sep =""))
  
}

stopCluster(cl)

system("rm metapop_temp/ploidy_file.tsv")
#retains only unique contig_pos combos and puts them into a file for filtering the mpileups next
system(paste("cat metapop_temp/* | sort -u > metapop_temp/unique_variant_contig_pos.txt", sep =""))
#cleans up the now redundant files
system("rm metapop_temp/*.tsv")

cl <- makeCluster(min(threads, detectCores(), length(newAlignments)))
registerDoParallel(cl)

#Make normal samtools pileups, filter them by the called SNV positions
silent <- foreach(i=1:length(newAlignments)) %dopar%{
  #Produces mpileups for each sample, concatenates the contig and position for joining purposes and to avoid having to do this later, and removes an uneccessary field
  
  system(paste(samtools_location, "samtools mpileup -I -f ", original_assembly, " ", newAlignments[i], 
               " | awk -F'\t' 'BEGIN{OFS=\"\t\";} {print $1\"_\"$2,$1,$2,$3,$4,$5}' | sort -k 1,1 | join metapop_temp/unique_variant_contig_pos.txt - > metapop_pileups/", 
               base_names[i], "_cleaned_pile.tsv",  sep =""))
  
}

stopCluster(cl)

system("rm -r metapop_temp")

print("done!")

cat("Refining SNVs into SNPs... ")

#Since SNPs are now called, they can be quality filtered, and base correction may be performed.

#These get used later, so it's good to declare them now.
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

pile_short <- list.files(path="metapop_pileups")
piles <- list.files(path="metapop_pileups", full.names = T)

useless <- "([^ATCG.,])"
refReplace <- "([.,])"

cl <- makeCluster(min(threads, detectCores(), length(piles)))
registerDoParallel(cl)
clusterExport(cl, c("library_location", "useless", "refReplace", "piles"))
clusterEvalQ(cl=cl, expr = library(data.table, lib.loc = library_location))
clusterEvalQ(cl=cl, expr = library(stringr, lib.loc = library_location))

# Read and process mpileup files in parallel
pile_counts <- foreach(i=piles) %dopar%{
  tmp <- fread(i, col.names = c("contig_pos", "contig", "pos", "ref_base", "depth", "pileup"))
  tmp$pileup <- str_to_upper(tmp$pileup)
  tmp$pileup <- str_replace_all(tmp$pileup, useless, "")
  tmp$pileup <- str_replace_all(tmp$pileup, refReplace, tmp$ref_base)
  tmp$a_ct <- str_count(tmp$pileup, "A")
  tmp$t_ct <- str_count(tmp$pileup, "T")
  tmp$c_ct <- str_count(tmp$pileup, "C")
  tmp$g_ct <- str_count(tmp$pileup, "G")
  tmp <- tmp[,-6]
  return(tmp)
}

stopCluster(cl)

rm(useless, refReplace)

#Fix sources before combining. They'll be used to split things later.
pile_counts <- lapply(1:length(pile_counts), function(x){
  tmp <- pile_counts[[x]]
  tmp$source <- unlist(substring(pile_short[x], 1, nchar(pile_short[x])-17))
  return(tmp)
})

pile_counts <- rbindlist(pile_counts)

if(time_series){
  #attempts to correct for issues with first sample calling names.
  if(!(first_sample %in% pile_counts$source)){
    cat("Exact first sample name not found. Attempting to correct...")
    if(sum(grepl(first_sample, newAlignments)) > 0){
      first_sample <- base_names[grepl(first_sample, newAlignments)]
      cat("done!\n")
    }else{
      cat("\n")
      print("Sample not found. Try checking metapop_run_settings/SNV_call_base_names.tsv to find the correct name to use. Use the name found in the second column, please.")
      quit(save="no")
    }
  }
  
  #Aggregate observations of each base from each pileup, per position, according to the original sample specified by the user.
  first_sample_corrective <- pile_counts[source==first_sample, lapply(list(depth, a_ct, t_ct, c_ct, g_ct), sum), by=contig_pos]
  colnames(first_sample_corrective)[2:6] = c("depth", "a_ct", "t_ct", "c_ct", "g_ct")
  first_sample_corrective$refBase <- c("A", "T", "C", "G")[max.col(first_sample_corrective[,3:6], ties.method = "random")]
  
  #We only care about the SNV positions from the first sample
  positions <- unique(first_sample_corrective$contig_pos)
  
  #Restrict things to those positions before aggregating
  pile_counts <- pile_counts[pile_counts$contig_pos %in% positions,]
  
  #Aggregate over these positions. Reference will be declared according to first sample corrective.
  aggregate_piles <- pile_counts[, lapply(list(depth, a_ct, t_ct, c_ct, g_ct), sum), by=contig_pos]

  #aggregate needs info on these positions from other samples, too, though
}else{
  #Aggregate observations of each base from each pileup, per position
  aggregate_piles <- pile_counts[, lapply(list(depth, a_ct, t_ct, c_ct, g_ct), sum), by=contig_pos]
}

#legibility
colnames(aggregate_piles)[2:6] = c("depth", "a_ct", "t_ct", "c_ct", "g_ct")

# Only retain positions matching human genome project standards - >=4  obs snps in >=2 pos and >= 1% rep for the SNP
# >3 obs SNPs
aggregate_piles <- aggregate_piles[(rowSums(aggregate_piles[,3:6] >= snp_obs_cutoff) > 1),]
# >= 1% representation - get proportions
aggregate_piles <- aggregate_piles[(rowSums((aggregate_piles[,3:6]/unlist(aggregate_piles[,2])) > min_rep) > 1),]

# MetaPop corrects the reference base according to highest abundance per position. This is for the purpose of the biostats later on.
aggregate_piles$refBase <- c("A", "T", "C", "G")[max.col(aggregate_piles[,3:6], ties.method = "random")]

pile_counts <- pile_counts[pile_counts$contig_pos %in% aggregate_piles$contig_pos,]

if(time_series){
  pile_counts$ref_base <- first_sample_corrective$refBase[match(pile_counts$contig_pos, first_sample_corrective$contig_pos)]
}else{
  pile_counts$ref_base <- aggregate_piles$refBase[match(pile_counts$contig_pos, aggregate_piles$contig_pos)] 
}

######Changing gears - time to figure out where genes are and operate on them

genes <- readDNAStringSet(geneFile)

s <- strsplit(names(genes), "[# \t]+") # split names by tab/space
genes <- data.table(matrix(unlist(s), ncol=5, byrow=T))

names(genes)[1:4] = c("contig_gene", "start", "end", "OC")

genes$start <- as.numeric((genes$start))
genes$end <- as.numeric((genes$end))

genes <- genes[,-5]

# Figure out what contig they come from, mostly for cleaning purposes
genes$parent_contig <- gsub("_\\d+$", "", genes$contig_gene)


# Might as well do this here and now - prepares the sequences to be processed by awk outside of this file - the flat format is just better for this
codon_bias_input <- data.table(seq = NA, OC = as.numeric(genes$OC), contig_gene = genes$contig_gene, parent_contig = genes$parent_contig)

# Figure out what bases each position has as SNPs
a_snp <- t_snp <- c_snp <- g_snp <- rep("", nrow(aggregate_piles))


a_snp[(aggregate_piles$a_ct > snp_obs_cutoff & aggregate_piles$a_ct/aggregate_piles$depth > min_rep & aggregate_piles$refBase != "A")] <- "A"
t_snp[(aggregate_piles$t_ct > snp_obs_cutoff & aggregate_piles$t_ct/aggregate_piles$depth > min_rep & aggregate_piles$refBase != "T")] <- "T"
c_snp[(aggregate_piles$c_ct > snp_obs_cutoff & aggregate_piles$c_ct/aggregate_piles$depth > min_rep & aggregate_piles$refBase != "C")] <- "C"
g_snp[(aggregate_piles$g_ct > snp_obs_cutoff & aggregate_piles$g_ct/aggregate_piles$depth > min_rep & aggregate_piles$refBase != "G")] <- "G"

aggregate_piles$snps <- paste0(a_snp, t_snp, c_snp, g_snp)

rm(a_snp, t_snp, c_snp, g_snp)

aggregate_piles$pos <- pile_counts$pos[match(aggregate_piles$contig_pos, pile_counts$contig_pos)]
aggregate_piles$contig <- pile_counts$contig[match(aggregate_piles$contig_pos, pile_counts$contig_pos)]


# This retains all the genes lacking a SNP, but on the same contig. This is on purpose.
genes <- genes[genes$parent_contig %in% aggregate_piles$contig,]

#Get the plausible snps into the split pile
pile_counts$snps <- aggregate_piles$snps[match(pile_counts$contig_pos, aggregate_piles$contig_pos)]

#Split genes into their parent contigs ahead of time
gspl <- genes[, list(list(.SD)), by = parent_contig]$V1
names(gspl) = unique(genes$parent_contig)
rm(genes)

#Positions are given relative to the contig.
gene_match <- mapply(function(pos, genes, CP){
  
  rows <- genes[genes$start <= pos & genes$end >= pos,]
  
  retval <- nrow(rows) > 0
  
  #An odd case where the contig has no predicted genes. Popped up in testing.
  if(identical(retval, logical(0))){
    return(NA)
  }
  
  if(retval){
    ##Figure out the codon, depending on the strand - count forwards if forwards, backwards otw.
    rows$codon <- ifelse(rows$OC > 0,((pos-rows$start)%/%3)+1, ((rows$end-pos)%/%3)+1)
    #Figure out the position within the codon, depending on the strand
    rows$pos_in_codon <- ifelse(rows$OC > 0,((pos-rows$start)%%3)+1, ((rows$end-pos)%%3)+1)
    
    rows$contig_pos <- CP
    
    return(list(rows))
  }
  else{
    return(NA)
  }
  
}, aggregate_piles$pos, gspl[aggregate_piles$contig], aggregate_piles$contig_pos)

#end testing

gene_match <- rbindlist(gene_match[!is.na(gene_match)])

# If there is a plausible match within a single sample, this flags the rows as having them.

# Set things up for fast joins
setkeyv(gene_match, c("contig_pos"))
setkeyv(pile_counts, c("contig_pos"))

# merges these tables with the info for entries in gene_match
genic_snps <- pile_counts[gene_match, nomatch=0]

setkeyv(genic_snps, c("source", "contig_gene", "OC", "codon"))
genic_snps[, link := .N > 1, by = key(genic_snps)]


if(!dir.exists("metapop_called_snps")){
  system("mkdir metapop_called_snps")
}

# Remove any missing cases from the join, as they are not useful.
genic_snps <- genic_snps[nchar(genic_snps$snps)>0,]

# We will need the non-genic snps later for diversity stats
pile_counts <- pile_counts[!pile_counts$contig_pos %in% genic_snps$contig_pos,]
# Also have to remove the missing cases
pile_counts <- pile_counts[nchar(pile_counts$snps)>0,]

fwrite(pile_counts, "metapop_called_snps/non_genic_snps.tsv", sep = "\t")

#sam format indicates that all reads in the file are and must be reported as primary strand match, with a flag indicating that a read was originally a match to a complement
#Although pileups also indicate this primary/complementary strand information, metapop counts all hits as belonging to the primary strand - this produces accurate counts of coverage and mutations
#However, it results in the relevant loci on reverse-strand genes seeing the opposite of the observed counts in the case of the recorded data, and thereby the observed reference base and SNP
#Therefore, these must be changed to their complements.
genic_snps[, ref_base := ifelse(OC == -1, chartr("ATCG", "TAGC", ref_base), ref_base)]
genic_snps[, snps := ifelse(OC == -1, chartr("ATCG", "TAGC", snps), snps)]

#only the relevant old counts
complement_swap <- genic_snps[OC < 0, list(a_ct, t_ct, c_ct, g_ct)]

#update the counts in the overall file
genic_snps[OC == -1, a_ct := complement_swap$t_ct]
genic_snps[OC == -1, t_ct := complement_swap$a_ct]
genic_snps[OC == -1, c_ct := complement_swap$g_ct]
genic_snps[OC == -1, g_ct := complement_swap$c_ct]

fwrite(genic_snps, "metapop_called_snps/called_snps.tsv", sep = "\t")

genic_snps[,pos_in_gene := ((genic_snps$codon-1)*3) + genic_snps$pos_in_codon]

genes <- readDNAStringSet(geneFile)

match_name <- unlist(lapply(names(genes), function(x){
  strsplit(x, split = " ")[[1]][1]
}))

gene_base_replace <- genic_snps[, list("locs" = list(as.integer(pos_in_gene)), "reps" = list(ref_base)), by = contig_gene]

#genes[[match(gene_base_replace$contig[i], match_name)]] <<- 

cl <- makeCluster(min(threads, detectCores(), length(newAlignments)))

clusterExport(cl, varlist = c("genes", "match_name", "gene_base_replace", "library_location"))
clusterEvalQ(cl, expr = library(Biostrings, lib.loc = library_location))

registerDoParallel(cl)

#Calls SNPs with no indels and quality filters them, reads in the results
replacements <- foreach(i=1:nrow(gene_base_replace)) %dopar%{
  #Produces SNV calls, filters these for SNV positions only (remove indels)  
  
  return(replaceLetterAt(genes[[match(gene_base_replace$contig_gene[i], match_name)]], gene_base_replace$locs[[i]], gene_base_replace$reps[[i]])) 
}

stopCluster(cl)

genes[match(gene_base_replace$contig, match_name)] <- replacements

codon_bias_input[, seq := as.character(genes)]

writeXStringSet(genes, filepath = "Metapop_Assembly/base_corrected_genes.fna")

run_parameters$setting[run_parameters$parameter == "Genes"] <- "Metapop_Assembly/base_corrected_genes.fna"

if(!dir.exists("metapop_codon_bias")){
  system("mkdir metapop_codon_bias")
}

fwrite(codon_bias_input, "metapop_codon_bias/codon_bias_sequences.tsv", sep = "\t", col.names = F)
rm(codon_bias_input)

#just for corrections
original_assembly <- as.character(run_parameters$setting[run_parameters$parameter=="Assembly"])
assembled_contigs <- readDNAStringSet(original_assembly)

base_replace <- aggregate_piles[, list("locs" = list(as.integer(pos)), "reps" = list(refBase)), by = contig]

match_name <- unlist(lapply(names(assembled_contigs), function(x){
  return(strsplit(x, split = " ")[[1]][1])
}))


cl <- makeCluster(min(threads, detectCores(), length(newAlignments)))

clusterExport(cl, varlist = c("assembled_contigs", "match_name", "base_replace", "library_location"))
clusterEvalQ(cl, expr = library(Biostrings, lib.loc = library_location))

registerDoParallel(cl)

#Calls SNPs with no indels and quality filters them, reads in the results
replacements <- foreach(i=1:nrow(base_replace)) %dopar%{
  #Produces SNV calls, filters these for SNV positions only (remove indels)  
  return(replaceLetterAt(assembled_contigs[[match(base_replace$contig[i], match_name)]], base_replace$locs[[i]], base_replace$reps[[i]]))
  
}

stopCluster(cl)

assembled_contigs[match(base_replace$contig, match_name)] <- replacements

writeXStringSet(assembled_contigs, filepath = "Metapop_Assembly/base_corrected_contigs.fna")

rm(match_name)

run_parameters$setting[run_parameters$parameter == "Assembly"] <- "Metapop_Assembly/base_corrected_contigs.fna"

fwrite(run_parameters, "metapop_run_settings/run_settings.tsv", sep = "\t")

#Way easier here than manually in awk
codon_initialize <- paste(paste("codons[\"", construct, "\"]=0;", sep = ""), collapse = " ")
codon_output <- paste(paste("codons[\"", construct, "\"]", sep = ""), collapse="\"\\t\"")

# The direct system call to make the codon bias part happen should work, but just in case a system doesn't like issuing a script from a script, it'll be saved for the user to bash
out <- paste("cut -f1 -d$'\\t' codon_bias_sequences.tsv | awk '{gsub(/.{3}/,\"& \")}1' | awk 'BEGIN{", codon_initialize, "} { for (key in codons) codons[key]=0; for ( i=1; i<NF+1; i++ ) codons[$i]+=1; print", codon_output, "}' > preprocessed_genes.tsv")
writeLines(c(paste0("cd ", getwd(),"/metapop_codon_bias"), out), "metapop_codon_bias/launch_bash.sh")
system("bash metapop_codon_bias/launch_bash.sh")

print("Done!")

run_summaries$end[run_summaries$stage=="SNP Calling"] <- as.character(Sys.time())

fwrite(run_summaries, "metapop_run_settings/run_dates.tsv", sep = "\t")
