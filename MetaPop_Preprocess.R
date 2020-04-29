options <- commandArgs(trailingOnly = T)

options <- options[!grepl("--args", options)]

print(paste("Start Preprocessing:", Sys.time()))

check_for_default <- integer(0)

directory_name <- which(grepl("-dir", options))
library_location <- which(grepl("-lib", options))
samtools_location <- which(grepl("-sam", options))
bcftools_location <- which(grepl("-bcf", options))
prodigal_location <- which(grepl("-prodigal", options))

quitNow <- 0

if(length(directory_name) > 1 | identical(check_for_default, directory_name)){
  print("Metapop requires the name of the directory where the BAM format aligned reads exist. This must be a single directory. Exiting")
  flush.console()
  quit(save="no")
}else{
  directory_name <- options[directory_name + 1]
}


if(length(samtools_location) > 1 | identical(check_for_default, samtools_location)){
  print("Metapop requires the installation location of samtools so that it may be called. Ignore this if samtools is in the path.")
  flush.console()
  samtools_location <- ""
  check_samtools <- system(paste(samtools_location, "samtools --version", sep =""), intern = T)[1]
  #
  if(grep("samtools 1.", check_samtools)){
    print("samtools found in environment. Proceeding")
  }else{
    quitNow <- 1 
  }
}else{
  samtools_location <- options[samtools_location + 1]
  
  
  samtools_location <- paste(samtools_location, "/", sep = "")
  
}

if(length(bcftools_location) > 1 | identical(check_for_default, bcftools_location)){
  print("Metapop requires the installation location of bcftools so that it may be called. Ignore this is bcftools is in the path.")
  flush.console()
  bcftools_location <- ""
  check_bcftools <- system(paste(bcftools_location, "bcftools --version", sep = ""), intern = T)[1]
  if(grep("bcftools 1.", check_bcftools)){
    print("bcftools found in environment. Proceeding")
  }else{
    quitNow <- 1 
  }
}else{
  bcftools_location <- options[bcftools_location + 1]
  
  bcftools_location <- paste(bcftools_location, "/", sep = "")
}

if(quitNow == 1){
  quit(save = "no")
}


if(length(library_location) > 1 | identical(check_for_default, library_location)){
  print("Metapop will use the default location for R libraries during this session. This will likely fail in HPC environments.")
  flush.console()
  library_location <- .libPaths()
}else{
  library_location <- options[library_location + 1]
}


#libraries
suppressMessages(suppressWarnings(library(doParallel, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(data.table, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(Rsamtools, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(Biostrings, lib.loc = library_location)))

check_samtools <- system(paste(samtools_location, "samtools --version", sep =""), intern = T)[1]
check_bcftools <- system(paste(bcftools_location, "bcftools --version", sep = ""), intern = T)[1]

#add in the bit
if(sum(grepl("1.1[0-9]+", strsplit(check_samtools, split = " ")[[1]][2]))>0){
  
  print("Newest Samtools")
}else{
  if(as.numeric(strsplit(check_samtools, split = " ")[[1]][2]) < 1.19){
    print("Samtools version out of date. Update to samtools version > 1.19 for Metapop to run.") 
  }
}

if(sum(grepl("1.1[0-9]+", strsplit(check_bcftools, split = " ")[[1]][2]))>0){
  print("Newest BCFtools")
}else{
  if(as.numeric(strsplit(check_bcftools, split = " ")[[1]][2]) < 1.19){
    print("BCFtools version out of date. Update to bcftools version > 1.2 for Metapop to run.") 
  }
}

#check for options or assume default
ANI <- which(grepl("-id", options))
min_read_length <- which(grepl("-min", options))
cov <- which(grepl("-cov", options))
depth <- which(grepl("-dep", options))
trunc <- which(grepl("-trunc", options))
threads <- which(grepl("-threads", options))
assem_file <- which(grepl("-assem", options))
genes_file <- which(grepl("-genes", options))
counts <- which(grepl("-ct", options))
global <- which(grepl("-global", options))

force_genes <- which(grepl("-force_genes", options))

setwd(directory_name)
print(paste("Working in:", getwd()))

#Block for getting assembled contigs
if(identical(assem_file, check_for_default)){
  print("MetaPop requires you to specify the contigs to which reads were aligned.")
  print("Please provide an absolute path to this file.")
  print("Exiting.")
  quit(save = "no")
  
}else{
  assem_file <- as.character(options[assem_file + 1])
  
  if(!file.exists(assem_file)){
    print(paste("Contigs file", assem_file, "not found. Exiting."))
    quit(save = "no")
  }
  
  if(!dir.exists("Metapop_Assembly")){
    system("mkdir Metapop_Assembly")
  }
  
  }


#
if(identical(genes_file, check_for_default)){
  print("No genes FASTA file specified. MetaPop will attempt to generate this file")
  
}else{
  genes_file <- as.character(options[genes_file + 1])
  
  if(!identical(force_genes, check_for_default)){
    print("Generating prodigal genes file with the specified name.")
  }else{
  
  if(!file.exists(genes_file)){
    print("No genes file found. MetaPop's microdiversity analyses cannot complete without this file. Because this was specified, MetaPop will exit.")
    print("Further microdiversity analyses cannot properly complete without this file. Is the path to this file correct, and does it exist?")
    print("You also may not use the -genes argument, and MetaPop will attempt to generate it for you.")
    quit(save = "no")
  }
  
  genes <- readDNAStringSet(genes_file)
  format_correct <- sum(grepl("#", names(genes)[1])> 0)
  
  if(format_correct){
    print("Genes are correctly formatted. Proceeding")
  }else{
    print("Genes are not formatted correctly. MetaPop requires them to be in prodigal FASTA format (prodigal -d). Exiting")
    quit(save = "no")
  }
  
  }
  
  
}

if(identical(check_for_default, prodigal_location) & (identical(genes_file, check_for_default) || !identical(force_genes, check_for_default))){
  print("Metapop will attempt to generate a genes file using prodigal if it doesn't already exist. This message only appears if -prodigal not specified. Metapop will assume this means prodigal exists in PATH.")
  flush.console()
  prodigal_location <- "prodigal"
}else{
  prodigal_location <- options[prodigal_location + 1]
  print(paste("Prodigal location specified as:", prodigal_location))
  flush.console()
}

if((identical(genes_file, check_for_default)) | !identical(force_genes, check_for_default)){
  print("Attempting to make genes file...")
  
  if(identical(genes_file, check_for_default)){
    genes_file <- "default_genes.fna"
  }
  
  try({
    system(paste0(prodigal_location, " -p meta -q -i ", assem_file, " -o temp_prod_out.txt -d Metapop_Assembly/", genes_file))
  })
  
  genes_file <- paste0("Metapop_Assembly/", genes_file)
}

if(file.exists("temp_prod_out.txt")){
  system("rm temp_prod_out.txt")
}

#Initial Folder Creation Start

if(!dir.exists("Metapop_BAMs")){
  system("mkdir Metapop_BAMs")
  system("mv *.bam Metapop_BAMs")
  
  original_aligned_reads <- list.files(path = "Metapop_BAMs", full.names = T)
  
}else{
  original_aligned_reads <- list.files(path = "Metapop_BAMs", full.names = T)
}

#Initial folder creation end


if(identical(global, check_for_default)){
  print("Calculating %ID over aligned region. Use flag -global to change this.")
  do_global_id <- F
}else{
  print("Calculating %ID over whole read, including unaligned region.")
  do_global_id <- T
}

if(identical(ANI, check_for_default)){
  print("Defaulting pct. ID to 95%")
  ANI <- 95
}else{
  ANI <- as.numeric(options[ANI + 1])
}

if(identical(min_read_length, check_for_default)){
  print("Defaulting minimum read length to 30 bp")
  min_read_length <- 30
}else{
  min_read_length <- as.numeric(options[min_read_length + 1])
}

if(identical(cov, check_for_default)){
  print("Defaulting horizontal coverage (percent of bases covered at least once in a genome to be considered) to 70%")
  cov <- 70
}else{
  cov <- as.numeric(options[cov + 1])
}

if(identical(trunc, check_for_default)){
  print("Defaulting depth truncation level to 10. Truncation level removes the top and bottom quantiles of depth of coverage over covered positions. Covers the middle 80% by default.")
  trunc <- 10
}else{
  trunc <- as.numeric(options[trunc + 1])
}

if(identical(depth, check_for_default)){
  print("Defaulting minimum average truncated depth (average depth of coverage over only covered positions between trunc and 1-trunc quantile, per genome) to 10")
  depth <- 10
}else{
  depth <- as.numeric(options[depth + 1])
}

if(identical(threads, check_for_default)){
  print("Defaulting minimum threads to 4. If there are fewer than 4 threads available, the lowest number will be used.")
  threads <- 4
}else{
  threads <- as.numeric(options[threads + 1])
}


###############################################################
############OPTIONS SECTION DONE, PROCESSING BEGINS############
###############################################################

if(!dir.exists("metapop_run_settings")){
  dir.create("metapop_run_settings")
}

run_stats <- data.table(parameter = c("Directory", "Samtools Location", "BCFTools Location", "Library Location", "Assembly", "Genes", "ID Cutoff", "Min. Read Length", "Coverage", "Depth", "Truncation", "Pct. ID Context", "Threads"), 
                        setting = c(directory_name, samtools_location, bcftools_location, library_location[1], assem_file, genes_file[1], ANI,min_read_length,cov,depth,trunc, ifelse(do_global_id, "Global", "Local") ,threads))

run_summaries <- data.table(stage = c("Preprocessing", "SNP Calling", "SNP Linkage", "Microdiversity", "Macrodiversity"), begin = c(as.character(Sys.time()), rep(NA, 4)), end = c(as.character(Sys.time()), rep(NA, 4)))

fwrite(run_summaries, "metapop_run_settings/run_dates.tsv", sep = "\t")

fwrite(run_stats, "metapop_run_settings/run_settings.tsv", sep = "\t")


#Acquire base names
base_names <- unlist(strsplit(original_aligned_reads, split = "Metapop_BAMs/"))[c(F,T)]
base_names <- substr(base_names, 1, nchar(base_names)-4)

if(!dir.exists("metapop_preprocessed_reads")){
  system("mkdir metapop_preprocessed_reads") 
}

if(!dir.exists("metapop_temp")){
  system("mkdir metapop_temp")
}

#The format for this might be a problem from assembly to assembly; I'll need to specify this assumes the contig is the first bit and that any detail is space-sep.
fastaLengths <- readDNAStringSet(assem_file)
fastaLengths <- data.table(contig = names(fastaLengths), num_bases = lengths(fastaLengths))

if(grepl(" ", fastaLengths$contig[1])){
 print("Your contigs contain a space, and this will result in faulty output from prodigal and likely everything else later.")
 print("If you do not understand this error, stop and replace spaces in the FASTA file contig names with underscores and realign your reads to this new file.")
 print("If you do understand this error, proceed with caution.")
}

fastaLengths$contig <- unlist(lapply(fastaLengths$contig, function(x){return(strsplit(x, split=" ")[[1]][1])}))

#Function that gets the middle chunk of depths. This chunk is the set of data falling between the trunc and 1-trunc quantiles, but without removing more than trunc % of the data from either end
#EXAMPLE: sequence {1,1,1,1,2,3,4,5,6,7,8,9,9,9,50} with trunc 10 will have lower quantile = 1, upper quantile = 9, and refuses to remove more than floor(14/10) = 1 value from either end.
#EXAMPLE: lower bound and upper bound indicate the data should be {2,3,4,5,6,7,8}, but this would remove more than allowed (one value). 
#EXAMPLE: Thus, {1,1,1,2,3,4,5,6,7,8,9,9,9} is the resulting set, and truncated avg depth = 5.41 for the whole set, instead of the pre-trunc avg of 8.28
get_TAD <- function(dataset, length, t = trunc){
  
  extras <- rep(0, length-length(dataset))
  
  if(length(extras) > 0){
    dataset <- c(extras, dataset)
  }
  
  dataset <- sort(dataset)
  quants <- quantile(dataset, c((t/100), ((100-t)/100)))
  dataLen <- length(dataset)
  
  #Don't remove more than trunc percent of the data from either end, anyway
  max_remove <- floor(length(dataset)/t)
  
  
  #Get the first and last position which are appropriately greater/less than the desired values to remove
  #The greater/less than must be present as occasionally the quantiles will be either decimals or not actually present in the data, depending on precisely what the data look like and chosen quantiles
  bottom_cap <- tail(which(dataset <= quants[1]), 1) + 1
  top_cap <- head(which(dataset >= quants[2]), 1) - 1
  
  #If the lower truncation would remove too much, remove only the allowed amount, else remove the lower truncation
  if(bottom_cap > max_remove){
    remove_lower <- max_remove
  }else{
    remove_lower <- bottom_cap
  }
  
  #If the upper truncation would remove too much, remove only the allowed amount, else remove the upper truncation.
  if(top_cap < (dataLen-max_remove)){
    remove_upper <- (dataLen-max_remove)
  }else{
    remove_upper <- top_cap
  }
  
  #Get the truncated depth between the lower and upper allowable bounds
  return(mean(dataset[remove_lower:remove_upper]))
  
}

#We'll use this for the final output of the following loop
if(!dir.exists("metapop_cov_and_depth")){
  system("mkdir metapop_cov_and_depth")
}

if(!dir.exists("metapop_depth_by_pos")){
  system("mkdir metapop_depth_by_pos")
}

#Using samtools, sorts raw reads, then fills the MD field (usually already done by bowtie2, but it's fast enough and is more robust to guarantee it here), indexes
cl <- makeCluster(min(threads, detectCores(), length(base_names)))
registerDoParallel(cl)

#This loop does a lot. It does the following in the following order:

#(1): Sorts original reads, ensures that pct ID can be identified, and cleans reads for those with % identity > ANI AND read length > min_read_length
#(2): uses samtools to get depths by contig and position from the intermediate BAMs made in (2) and saves these
#(3): uses awk to find the number of unique positions covered at least once, per contig, via the output of (3) and reads these into R
#(4): finds the horizontal coverage of each contig by comparing the number of positions covered on each contig found in (4) to the original references pulled into R before the loop as 'fastaLengths'
#(4b): If there are any contigs with high enough coverage, writes all contigs and their num. pos. cov and pct. cov
#(4c): removes files with no high-coverage contigs from further analysis; files from (5b) are never written if this is so
#(5): retains the contigs with horizontal coverage % > cov and reads the depth per position and corresponding contig for these sufficiently high coverage contigs into R
#(5b): deletes the by-position depth file to free up disk space
#(6): uses R to group the depths by contig and calculate their truncated average depth between quantiles trunc and 1-trunc; details in the function descr. above
#(6b): writes the summarized truncated depth file; these will only have contigs with enough horizontal coverage, and so will have fewer contigs than (5b) in most cases
#(7): retains the names of the contigs with truncated average depth >= depth and returns this vector of contigs

silent <- foreach(i=1:length(base_names), .packages = c("Rsamtools", "data.table")) %dopar% {
  
  system(paste(samtools_location, "samtools sort ", original_aligned_reads[i], " > metapop_preprocessed_reads/", base_names[i], ".sorted.bam", sep = ""), wait = T)
  system(paste(samtools_location, "samtools calmd -b metapop_preprocessed_reads/", base_names[i], ".sorted.bam ", assem_file, " > metapop_preprocessed_reads/", base_names[i], ".ANI_Prep.bam", sep = ""), wait = T)
  
  system(paste("rm metapop_preprocessed_reads/", base_names[i], ".sorted.bam", sep = ""), wait = T)
  
  if(do_global_id){
  #global alignment %ID
  system(paste(samtools_location, "samtools view -h metapop_preprocessed_reads/", base_names[i], ".ANI_Prep.bam | grep -E '@|MD:Z:' | awk '{ if( substr ($0, 0, 1) == \"@\" ){print $0} else{ match($0, /MD:Z:[0-9A-Z\\^]*/,m ); split(m[0],v,/[\\^:]/); nmatch = split(v[3],vmatch, /[^0-9]/); cmatch=0; for(i=1; i<=nmatch; i++) cmatch+=vmatch[i]; if( length($10)>=", min_read_length, " && (cmatch*100)/ (length($10))>=",ANI, " ){print($0)};}}' | " ,samtools_location, "samtools view -bS - > ","metapop_temp/", base_names[i], "_ANI_filt.bam", sep = ""))
  }else{
  #Local alignment %ID
  system(paste(samtools_location, "samtools view -h metapop_preprocessed_reads/", base_names[i], ".ANI_Prep.bam | grep -E '@|MD:Z:' | awk '{ if( substr ($0, 0, 1) == \"@\" ){print $0} else{ match($0, /MD:Z:[0-9A-Z\\^]*/,m ); split(m[0],v,/[\\^:]/); nmatch = split(v[3],vmatch, /[^0-9]/); cmatch=0; for(i=1; i<=nmatch; i++) cmatch+=vmatch[i]; if( length($10)>=", min_read_length, " && (cmatch*100/(cmatch+nmatch-1))>=",ANI, " ){print($0)};}}' | " ,samtools_location, "samtools view -bS - > ","metapop_temp/", base_names[i], "_ANI_filt.bam", sep = ""))
  }
    
  system(paste("rm metapop_preprocessed_reads/", base_names[i], ".ANI_Prep.bam", sep = ""), wait = T)
  
  system(paste(samtools_location, "samtools index metapop_temp/", base_names[i], "_ANI_filt.bam", sep = ""), wait = T)
  
  #Translate the bams to essentially an mpileup format, then to number of positions covered. Only read contigs
  system(paste(samtools_location, "samtools depth metapop_temp/", base_names[i], "_ANI_filt.bam > metapop_depth_by_pos/", base_names[i], ".depth_by_pos.tsv",  sep = ""), wait = T)
  
  #convert these into an R-friendly format
  depthSets <- fread(paste("metapop_depth_by_pos/", base_names[i], ".depth_by_pos.tsv", sep = ""), select = c(1,3))
  
  
  if(nrow(depthSets) == 0){
    system(paste0("rm ", "metapop_depth_by_pos/", base_names[i], ".depth_by_pos.tsv"))
    system(paste0("rm ", "metapop_temp/", base_names[i], "_ANI_filt.bam"))
    system(paste0("rm ", "metapop_temp/", base_names[i], "_ANI_filt.bam.bai"))
    return(NA)
  }
  
  colnames(depthSets) = c("contig", "depths")
  
  cov_values <- depthSets[, length(depths), by = contig]
  colnames(cov_values) = c("contig", "numCov")
  
  cov_values[, percent_cov := (numCov/(fastaLengths$num_bases[match(contig, fastaLengths$contig)]))*100]
  
  #Get the tad N values using data.table and select the contigs that pass, refilter the samtools stuff.
  depthSets <- depthSets[, get_TAD(depths, fastaLengths$num_bases[match(contig, fastaLengths$contig)]) , by = contig]
  colnames(depthSets)[2] = "trunc_depth"
  
  cov_values$depth <- depthSets$trunc_depth[match(cov_values$contig, depthSets$contig)]
  
  #Output depth summary from cov. passing contigs
  fwrite(cov_values, paste("metapop_cov_and_depth/", base_names[i], ".cov_and_depth_by_contig.tsv", sep=""), sep = "\t", col.names = F)
  
  #Only passing horiz covs
  depthSets <- depthSets[depthSets$contig %in% cov_values$contig[cov_values$percent_cov >= cov],]
  
  #Reads are already cleaned. now we clean contigs
  contigsToKeep <- unique(depthSets$contig[depthSets$trunc_depth >= depth])
  rm(depthSets)
  
  if(length(contigsToKeep) >= 0){
    
    filterClean <- FilterRules(list(covered=function(y){
      y$rname %in% contigsToKeep
    }))
    
    paramClean = ScanBamParam(what=c("rname"))
    
    #Filter the bams for high enough ANI values
    filterBam(paste0("metapop_temp/", base_names[i], "_ANI_filt.bam"), destination = paste("metapop_preprocessed_reads/", base_names[i], "_preprocessed.bam", sep = ""), filter = filterClean, param = paramClean)
    
  }
  
  system(paste("rm ", paste0("metapop_temp/", base_names[i], "_ANI_filt.bam"), sep = ""))
  system(paste("rm ", paste0("metapop_temp/", base_names[i], "_ANI_filt.bam.bai"), sep = ""))
  
}

stopCluster(cl)

system("rm -r metapop_temp")

run_summaries$end[run_summaries$stage=="Preprocessing"] <- as.character(Sys.time())

fwrite(run_summaries, "metapop_run_settings/run_dates.tsv", sep = "\t")

print(paste("End Preprocessing:", Sys.time()))







