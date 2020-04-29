options <- commandArgs(trailingOnly = T)

check_for_default <- integer(0)

directory_name <- which(grepl("-dir", options))

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


run_parameters <- read.csv("metapop_run_settings/run_settings.tsv", sep = "\t", stringsAsFactors = F)

library_location <- as.character(run_parameters$setting[run_parameters$parameter=="Library Location"])

suppressMessages(suppressWarnings(library(doParallel, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(data.table, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(stringr, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(Biostrings, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(Rsamtools, lib.loc = library_location)))

run_summaries <- fread("metapop_run_settings/run_dates.tsv", sep = "\t")

run_summaries$begin[run_summaries$stage=="SNP Linkage"] <- as.character(Sys.time())

threads <- as.numeric(run_parameters$setting[run_parameters$parameter=="Threads"])

linked_snps <- fread("metapop_called_snps/called_snps.tsv", sep = "\t")
linked_snps <- linked_snps[link == T,]


#Use this for cleaning out invalid codons containing Ns.
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

if(nrow(linked_snps)==0){
  
  print("There are no potentially linked snps in your data. This process is not applicable, but other steps will still function.")
  quit("no")
  
}

colnames(linked_snps) = c("contig_pos",	"contig",	"pos",	"ref_base",	"depth",	"a_ct",	"t_ct",	"c_ct",	"g_ct",	"source",	"snps",	"contig_gene",	"start",	"end",	"OC",	"codon",	"pos_in_codon",	"link")

#For this step, most of this data is useless. Removing it now is good.
linked_snps <- linked_snps[,c(1:4, 10:12, 15:17)]

#Order things for consistent behavior later
linked_snps <- linked_snps[order(source, contig_gene, OC, codon, pos),]

#Gather information into a wider format so that shaping to long is easier later

#Better group_by code
setkeyv(linked_snps, c("source", "contig_gene", "OC", "codon"))
collections <- linked_snps[, list(pos=list(pos), refs=list(ref_base), snps=list(snps), pos_in_codon=list(pos_in_codon)), by = key(linked_snps)]

# Based on the method of IDing putative linked SNPs, some false positives get through. This corrects that.
collections <- collections[lengths(collections$pos)>1,]

setkeyv(collections, c("source", "contig_gene", "OC", "codon"))

# Clean out the false positives here, too
linked_snps[, hasMatch := length(ref_base) > 1, by = key(linked_snps)]
linked_snps <- linked_snps[linked_snps$hasMatch,]
linked_snps[,hasMatch := NULL]

# no need to multi-read data
unique_codons <- unique(linked_snps, by = key(linked_snps))

# If it's original strand, the max pos. in genome is >= pos, by 2 is pos in codon = 1, 1 if 2, 0 if 3. 
#For the compliment, the max pos. in gen. is when the pos. in codon is 1, so we add when pos. in contig is 3 instead
unique_codons$max_pos <- ifelse(unique_codons$OC > 0, 
                                unique_codons$pos+(c(0,1,2)[match(unique_codons$pos_in_codon, c(3,2,1))]),
                                unique_codons$pos+(c(2,1,0)[match(unique_codons$pos_in_codon, c(3,2,1))]))

unique_codons$min_pos <- ifelse(unique_codons$OC > 0, 
                                unique_codons$pos-(c(0,1,2)[match(unique_codons$pos_in_codon, c(1,2,3))]),
                                unique_codons$pos-(c(2,1,0)[match(unique_codons$pos_in_codon, c(1,2,3))]))

unique_unique_sources <- unique(unique_codons$source)

# Extracts only the seqs from each sample which correspond to covering a region where a potentially linked snp might be occurring
# Then, processes the reads to extract the relevant codon (where it is actually found), and returns counts of each unique codon occurrence

# Pileup automatically considers read call quality in generating its counts. Typically, there will be a few more reads returned by this than the pileup indicates
cl <- makeCluster(min(detectCores(), threads, length(unique_unique_sources)))
registerDoParallel(cl)
clusterExport(cl=cl, c("library_location", "unique_codons", "unique_unique_sources", "construct"))
suppressMessages(clusterEvalQ(cl=cl, expr = suppressMessages(suppressWarnings(library(doParallel, lib.loc = library_location)))))
suppressMessages(clusterEvalQ(cl=cl, expr = suppressMessages(suppressWarnings(library(stringr, lib.loc = library_location)))))
suppressMessages(clusterEvalQ(cl=cl, expr = suppressMessages(suppressWarnings(library(Rsamtools, lib.loc = library_location)))))
suppressMessages(clusterEvalQ(cl=cl, expr = suppressMessages(suppressWarnings(library(Biostrings, lib.loc = library_location)))))


#Reverse complements are not being done here - this is strictly observing counts of occurrences of particular SNPs relative to their position in the genome according to the primary strand.
#The task of correcting these codons happens in microdiversity
readSets <- foreach(i=unique_unique_sources) %dopar% {
  
  subset <- unique_codons[unique_codons$source == i,]
  
  gr <- GRanges(seqnames = subset$contig, ranges = IRanges(start = subset$min_pos, end = subset$max_pos))
  
  param <- ScanBamParam(what = c("pos", "cigar", "seq"), which = gr)
  
  reads <- scanBam(paste("metapop_preprocessed_reads/", i, "_preprocessed.bam", sep = ""), param = param)
  
  #end N cleaning
  
  collected_reads <- lapply(1:length(reads), function(y){
      
      if(length(reads[[y]]$cigar) > 0){
      
      #The set of sequences without insertions or deletions
      easy_cases <- lengths(str_extract_all(reads[[y]]$cigar, "[A-Z]")) == 1
      
      #We need the +1 because this would end up zero-indexing the read, and R 1-indexes
      easy_subs <- substring(reads[[y]]$seq[easy_cases], subset$min_pos[y]-reads[[y]]$pos[easy_cases]+1, subset$max_pos[y]-reads[[y]]$pos[easy_cases]+1)
      
      #TODO - here would be where a rev-comp is done for a negative strand read
      
      #The way that reads are extracted, if the start is less than max_pos and the end is greater than min_pos, the reads comes in. This little cleanup takes care of that.
      easy_subs <- easy_subs[nchar(easy_subs) == 3]
      
      if(sum(!easy_cases) > 0){  
      hard_cases_cig <- reads[[y]]$cigar[!easy_cases]
      hard_cases_reads <- as.character(reads[[y]]$seq[!easy_cases])
      hard_cases_pos <- reads[[y]]$pos[!easy_cases]
      
      cigar_nums <- lapply(strsplit(hard_cases_cig, split = "[A-Z]"), as.numeric)
      cigar_letters <- lapply(strsplit(hard_cases_cig, split = "[0-9]+"), function(w){return(w[2:length(w)])})
      
      hard_subs <- mapply(function(min_pos, read, pos, nums, letters){
        
        #These are the positions, relative to the start of the read, that need returned
        pulls <- (min_pos-pos+1):(min_pos-pos+3)
        
        #expands the cigar string into a sequence of letters matching the counts of each subsequent cigar match
        cig_by_pos <- rep(letters, times=nums)
        
        #Actual read bases are either a match or an insertion. This is the set that matches those
        matches_read <- cig_by_pos[cig_by_pos == "M" | cig_by_pos == "I"]
        
        #This is the set of positions that contribute to the pos. in genome
        cig_by_pos <- cig_by_pos[cig_by_pos == "M"|cig_by_pos=="D"]
        
        #Reduce the read to matches only, at the correct positions
        split_read <- strsplit(read, split="")[[1]]
        split_read <- split_read[matches_read=="M"]
        
        #Get all the positions, no skips
        corrected_position <- 1:sum(nums[letters=="M"|letters=="D"])
        corrected_position <- corrected_position[cig_by_pos == "M"]
        
        
        read <- paste0(split_read[match(pulls, corrected_position)], collapse = "")
        
        return(read)
        
        
      }, subset$min_pos[y], hard_cases_reads, hard_cases_pos, cigar_nums, cigar_letters)
      
      #Same cleanup as before. This also handles cases where there was a deletion in the codon; 
      # functionally the codon is no longer operating in these cases, as a base is missing, and they should be excluded.
      hard_subs <- hard_subs[nchar(hard_subs) == 3]
      
      fullSet <- c(easy_subs, hard_subs)
      fullSet <- fullSet[fullSet %in% construct]
      
      # We only need the counts of each type of occurrence, not the full-fledged set.
      return(rle(sort(fullSet)))
      
      } else {
        
        # And if there are no hard cases, this is fine
        fullSet <- easy_subs
        fullSet <- fullSet[fullSet %in% construct]
        
        return(rle(sort(fullSet)))
        
      }
      
      } else {
        
        return(rle(integer(0)))
        
      }
      
      
      
    })
    
  return(collected_reads)
  
}

stopCluster(cl)

#Translate the list of lists returned by the above loop into a single list with the same length as the unique_codons - this is the raw data we're working from
readSets <- do.call(c, readSets)

# For each codon on the same gene/strand in linked_snps, there is a single reference base at each variant position, and between 1 and 3 snps.
# To be thorough, each pairwise combination of SNPs should be interrogated 
# - up to 27 comparisons if there are 3 snps at each position in a codon, though this is very unikely.

# For the pusposes of creating a readable file that is convenient for a user and for extracting data later, linked_snps is converted to a format where the reference codon is compared to a
# specific snp combination, per line. Thus, each pairwise combo takes up a single line.

ref_set <- function(refs, positions){
  
  if(length(refs)==3){
    return(paste0(refs, collapse=""))
  } else {
    
    if(all(positions == c(1, 3)) | all(positions == c(3,1))){
      
      return(paste0(refs[1], "*", refs[2], collapse = ""))
      
    }
    
    if( all(positions == c(1,2)) ){
      
      return(paste0(refs[1], refs[2], "*", collapse = ""))
      
    }
    
    if( all(positions == c(2,1)) ){
      
      return(paste0("*", refs[1], refs[2], collapse = ""))
      
    }
    
    if( all(positions == c(2,3)) ){
      
      return(paste0("*", refs[1], refs[2], collapse = ""))
      
    }
    
    if( all(positions == c(3,2)) ){
      
      return(paste0(refs[1], refs[2], "*", collapse = ""))
      
    }
    
  }
  
}

collections$refs <- mapply(ref_set, collections$refs, collections$pos_in_codon)

# Gets every combo of >= 2 snps on the same codon; for codons with 3 snps, this includes each subset of 2 snps
snp_expand <- function(snps, positions){
  
  if(length(snps) == 3){
    
    set1 <- strsplit(snps[1], split = "")[[1]]
    set2 <- strsplit(snps[2], split = "")[[1]]
    set3 <- strsplit(snps[3], split = "")[[1]]
    
    snps <- expand.grid(set1, set2, set3)
    
    snps <- c(paste0(snps$Var1, snps$Var2, snps$Var3),
              paste0("*", snps$Var2, snps$Var3),
              paste0(snps$Var1, "*", snps$Var3),
              paste0(snps$Var1, snps$Var2, "*"))
    
    return(list(snps))
    
  } else {
    
    set1 <- strsplit(snps[1], split = "")[[1]]
    set2 <- strsplit(snps[2], split = "")[[1]]
    snps <- expand.grid(set1, set2)
    
    if(all(positions == c(1, 3)) | all(positions == c(3,1))){
      
      snps <- paste0(snps$Var1, "*", snps$Var2)
      return(list(snps))
      
    }
    
    if( all(positions == c(1,2)) ){
      
      snps <- paste0(snps$Var1, snps$Var2, "*")
      return(list(snps))  
      
    }
    
    if( all(positions == c(2,1)) ){
      
      snps <- paste0("*",snps$Var1, snps$Var2)
      return(list(snps))
      
    }
    
    if( all(positions == c(2,3)) ){
      
      snps <- paste0("*",snps$Var1, snps$Var2)
      return(list(snps))
      
    }
    
    if( all(positions == c(3,2)) ){
      
      snps <- paste0(snps$Var1, snps$Var2, "*")
      return(list(snps))  
      
    }
    
    
  }
  
}

flattened_snps <- mapply(snp_expand, collections$snps, collections$pos_in_codon)

repetitions <- lengths(flattened_snps)

# This contains all the info needed to access the associated mined reads and perform a fisher's exact test on each set of candidate linked SNPs
long_snp_table <- data.table(data_ref = rep(1:length(readSets), times = repetitions),
                             source = rep(collections$source, times = repetitions),
                             contig_gene = rep(collections$contig_gene, times = repetitions),
                             OC = rep(collections$OC, times = repetitions),
                             codon = rep(collections$codon, times = repetitions),
                             refs = rep(collections$refs, times = repetitions),
                             snp = unlist(flattened_snps))

# lots of data we no longer need
rm(linked_snps, collections, unique_codons, unique_unique_sources)

#Occasionally mining snps will produce no reads covering the spot. This is OK, but must be cleaned
lacking_data <- which(unlist(lapply(readSets, function(x){length(x$values)>0})))

#We don't want to mess with the reference mined reads, because this would require further correction. Instead, the empty indices are simply skipped.
long_snp_table <- long_snp_table[long_snp_table$data_ref %in% lacking_data,]

match_ref <- strsplit(long_snp_table$refs, split = "")
match_snp <- strsplit(long_snp_table$snp, split = "")

calculate_linkage <- function(data_index, split_ref, split_snp, strand){
  
  #Get the associated data
  data <- readSets[[data_index]]
  
  if(strand < 0){
    data$values <- lapply(data$values, function(x){
        
        reversed <- intToUtf8(rev(utf8ToInt(x)))
        complemented <- chartr("ATCG", "TAGC", reversed)
        
        return(complemented)
    })  
  }
  
  # transform it into a form that matches the test ref/snps
  data$values <- lapply(data$values, function(x){return(strsplit(x, split="")[[1]])})
  
  #Check for 3 snp case
  if(sum(split_ref == "*") == 0){
    
    #Convenience
    total <- sum(data$lengths)
    
    all_ref <- sum(data$lengths[which(unlist(lapply(data$values, function(x){
      return(all(x == split_ref))
    })))])
    all_snp <- sum(data$lengths[which(unlist(lapply(data$values, function(x){
      return(all(x == split_snp))
    })))])
    
    if(length(all_ref)==0){all_ref <- 0}
    if(length(all_snp)==0){all_snp <- 0}
    
    # randomize whether an odd count of off-diag cases will favor top right or bottom left in contingency table.
    randomizer <- (runif(1,0,1) > .5)
    
    #There's not a good way to do this with 3; this is a worst-case scenario for demonstrating linkage
    off_diag <- total-(all_ref+all_snp)
    
    first_is_ref <- ifelse(randomizer, floor(off_diag/2), ceiling(off_diag/2))
    second_is_ref <- off_diag-first_is_ref
    
    
  } else {
  
    
    all_ref <- sum(data$lengths[which(unlist(lapply(data$values, function(x){
      return(sum(x == split_ref)==2)
    })))])
    all_snp <- sum(data$lengths[which(unlist(lapply(data$values, function(x){
      return(sum(x == split_snp)==2)
    })))])
    
    if(length(all_ref)==0){all_ref <- 0}
    if(length(all_snp)==0){all_snp <- 0}
    
    # for the comparison of first ref/second ref, knowing where the star is matters
    missing <- which(split_ref == "*")
    
    if(missing == 1){first_place=2; second_place=3}
    if(missing == 2){first_place=1; second_place=3}
    if(missing == 3){first_place=1; second_place=2}
  
    
    first_is_ref <- sum(data$lengths[which(unlist(lapply(data$values, function(x){
      return((x[first_place] == split_ref[first_place] & x[second_place] == split_snp[second_place]))
    })))])
    second_is_ref <- sum(data$lengths[which(unlist(lapply(data$values, function(x){
      return((x[first_place] == split_snp[first_place] & x[second_place] == split_ref[second_place]))
    })))])
    
    if(length(first_is_ref)==0){first_is_ref <- 0}
    if(length(second_is_ref)==0){second_is_ref <- 0}
    
  }
  


  return(matrix(c(all_ref, all_snp, first_is_ref, second_is_ref), nrow=1))
  
}

linkage_results <- data.table(t(mapply(calculate_linkage, long_snp_table$data_ref, match_ref, match_snp, long_snp_table$OC)))
names(linkage_results) = c("ref_count", "snp_count", "ref_first", "ref_second")

long_snp_table <- cbind(long_snp_table, linkage_results)

linkage_test <- function(all_ref, second_is_ref, first_is_ref, all_snp){

  contingency_table <- matrix(c(all_ref, second_is_ref, first_is_ref, all_snp), nrow = 2, ncol = 2)
  
  if(any(contingency_table < 0)){
	return(list(NA, NA))
  }


  p <- fisher.test(contingency_table)$p
  top <- (all_ref*all_snp)-(second_is_ref*first_is_ref)
  botA <- sqrt(all_ref+first_is_ref)
  botB <- sqrt(all_ref+second_is_ref)
  botC <- sqrt(all_snp+first_is_ref)
  botD <- sqrt(all_snp+second_is_ref)
  botDiv <- botA*botB*botC*botD
  
  phi <- top/botDiv
  
  return(list(p, phi))
  
}

linkage_results <- mapply(linkage_test, long_snp_table$ref_count, long_snp_table$ref_second, long_snp_table$ref_first, long_snp_table$snp_count)

long_snp_table$fisher_p <- unlist(linkage_results[c(T,F)])
long_snp_table$phi_coef <- unlist(linkage_results[c(F,T)])

long_snp_table <- long_snp_table[!is.na(long_snp_table$fisher_p),]

if(!dir.exists("metapop_linked_snps")){
  system("mkdir metapop_linked_snps")
}

#test this
long_snp_table[, three_snps := grepl("*", snp)]

fwrite(long_snp_table, "metapop_linked_snps/linked_snp_results.tsv", sep = "\t")

run_summaries$end[run_summaries$stage=="SNP Linkage"] <- as.character(Sys.time())

fwrite(run_summaries, "metapop_run_settings/run_dates.tsv", sep = "\t")
