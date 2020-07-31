options <- commandArgs()

options <- options[!grepl("--args", options)]

check_for_default <- integer(0)

#check for help
help <- which(grepl("-help", options))

if(!identical(help, check_for_default)){
  print("Welcome to MetaPop.")
  print("")
  print("usage: MetaPop.R [OPTS] -dir -assem -ct")
  print("")
  print("Program options:")
  print("-preprocess_only : Flag indicating to filter reads for %ID, length, and depth of coverage and stop.")
  print("-micro_only : Flag indicating to perform only macrodiversity calculations. Assumes preprocess has been done.")
  print("-macro_only : Flag indicating to perform only microdiversity calculations. Assumes preprocess has been done.")
  print("-viz_only : Flag indicating to only produce visualizations. Assumes preprocess has been done.")
  
  print("-no_micro Flag : indicating to skip microdiversity and only perform preprocess and macrodiversity")
  print("-no_macro Flag : indicating to skip macrodiversity and only perform preprocess and microdiversity")
  print("-no_viz Flag : indicating to not attempt to visualize results.")
  print("")
  
  print("Mandatory arguments:")
  print("-dir : Directory containing mapped read files in BAM format")
  print("-assem : Absolute path to assembled contigs")
  print("-ct : Absolute path to counnts normalization file (only mandatory if performing macrodiversity calculations)")
  print("")
  print("Preprocessing Arguments:")
  print("-sam : absolute path to samtools - optional if samtools is in PATH")
  print("-bcf : absolute path to bcftools - optional if samtools is in PATH")
  print("-prodigal : absolute path to prodigal gene prediction tool - optional if prodigal is in PATH or if -genes specified with existing prodigal FASTA format gene calls")
  print("-genes : absolute path to prodigal FASTA format gene file for assembled contiigs. May be left absent if you want metapop to generate this file.")
  print("-id INT : reads below this percent identity (mismatch/alignment length) are removed. Use -global to calculate as (mismatch/read length). Default 95")
  print("-min INT : reads shorter than this are removed. Default 30.")
  print("-cov INT : contigs with breadth of coverage (#bases covered/contig length) less than this are removed from microdiversity. Default 70")
  print("-dep INT : contigs with truncated average depth of coverage (mean of the 10th - 90th percentile depths of coverage) less than this are removed from microdiversity. Default 10.")
  print("-trunc INT : sets the percentiles at which depths of coverage will be truncated for depth. Default 10.")
  print("-threads INT : Number of threads to parallelize processes for. Default 4.")
  print("-global FLAG : Flag indicating to calculate percent identity using read length instead of alignment length.")
  print("-force_genes FLAG: Flag indicating that the file specified with -genes doesn't exist, and requests a file with this name be generated using prodigal.")
  print("")
  print("Variant Calling Arguments:")
  print("-first STRING : prefix of sample to be used as the SNP reference point for all other samples. A prefix is the set of characters before .bam extension in a sample, e.g. a_file.bam has prefix a_file.")
  print("-obs INT : Minimum number of observations of a variant allele at a particular locus for it to be called a SNP. Default 4")
  print("-rep INT : Minimum percent of the population a variant allele at a particular locus must represent for it to be called a SNP. Default 1")
  print("-var_qual INT : Minimum PHRED score for a base to be used in initial variant calling. Default 20.")
  print("")
  print("Microdiversity Arguments:")
  print("-subsamp INT : SNP loci will be subsampled proportionallydown to this depth of coverage for microdiversity calculations. Default 10.")
  
  
  print("Macrodiversity Arguments:")
  print("-complete_bact FLAG : Assumes that each contig supplied is a complete bacterial genome. Lowers the threshold of detection from 70% coverage to 20% coverage.")
  print("-min_det INT : Minimum percent of bases covered for a contig to be considered detected. Default 70; will change if -cov is set to a different value of -complete_bact is used.")
  print("-min_bp INT : Minimum number of positions required to be covered for a contig to be considered detected. Default 5000.")
  
  
  print("Visualization arguments:")
  print("-all FLAG : Metapop will print all contigs from microdiversity results. This will likely take a long time. Default prints top 3 genomes by highest % of genes under positive selection in each sample.")
  print("-snp_scale [local, global, both] : Metapop will print microdiversity results using SNPs assessed at the local (per sample) or global (across all samples) levels, or both. Defaults to local.")
  
  
  quit(save = "no")
}

directory_name <- which(grepl("-dir", options))

if(length(directory_name) > 1 | identical(check_for_default, directory_name)){
  print("Metapop requires the name of the directory where the BAM format aligned reads exist. This must be a single directory. Exiting")
  flush.console()
  quit(save="no")
}else{
  directory_name <- options[directory_name + 1]
}
setwd(directory_name)
print(paste("Working in:", getwd()))

#Modules:
print("Loading MetaPop Modules. This may take a minute.")
#Preproc
metapop_preprocess <- function(options){
  
  check_for_default <- integer(0)
  
  library_location <- which(grepl("-lib", options))
  samtools_location <- which(grepl("-sam", options))
  bcftools_location <- which(grepl("-bcf", options))
  prodigal_location <- which(grepl("-prodigal", options))
  
  print(paste("Start Preprocessing:", Sys.time()))
  
  check_for_default <- integer(0)
  
  quitNow <- 0
  
  
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
    library_location <- .libPaths()
  }else{
    library_location <- options[library_location + 1]
  }
  
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
  global <- which(grepl("-global", options))
  
  force_genes <- which(grepl("-force_genes", options))
  
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
  
  
  
  
  
  
  
}

#Call SNP
metapop_call_snps <- function(options){
  
  check_for_default <- integer(0)
  
  
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
    print("Defaulting vcf min phred call quality to 30")
    var_quality <- 30
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
  
  
  run_parameters <- read.csv("metapop_run_settings/run_settings.tsv", sep = "\t", stringsAsFactors = F)
  
  library_location <- as.character(run_parameters$setting[run_parameters$parameter=="Library Location"])
  
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
  clusterExport(cl, c("library_location", "useless", "refReplace", "piles"), envir=environment())
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
  
  clusterExport(cl, varlist = c("genes", "match_name", "gene_base_replace", "library_location"), envir=environment())
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
  
  clusterExport(cl, varlist = c("assembled_contigs", "match_name", "base_replace", "library_location"), envir=environment())
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
}

#Mine Reads
metapop_mine_reads <- function(options){
  
  check_for_default <- integer(0)
  
  run_parameters <- read.csv("metapop_run_settings/run_settings.tsv", sep = "\t", stringsAsFactors = F)
  
  library_location <- as.character(run_parameters$setting[run_parameters$parameter=="Library Location"])
  
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
  clusterExport(cl=cl, varlist = c("library_location", "unique_codons", "unique_unique_sources", "construct"), envir = environment())
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
}

#Microdiversity
metapop_microdiv <- function(options){
  
  
  check_for_default <- integer(0)
  

  
  run_parameters <- read.csv("metapop_run_settings/run_settings.tsv", sep = "\t")
  run_parameters$parameter <- as.character(run_parameters$parameter)
  run_parameters$setting <- as.character(run_parameters$setting)
  
  library_location <- run_parameters$setting[run_parameters$parameter == "Library Location"]
  threads <- as.numeric(run_parameters$setting[run_parameters$parameter == "Threads"])
  
  original_assembly <- as.character(run_parameters$setting[run_parameters$parameter=="Assembly"])
  geneFile <- as.character(run_parameters$setting[run_parameters$parameter=="Genes"])
  
  cov_cutoff <- as.numeric(run_parameters$setting[run_parameters$parameter=="Coverage"])
  depth_cutoff <- as.numeric(run_parameters$setting[run_parameters$parameter=="Depth"])
  
  sub_samp <- which(grepl("-subsamp", options))
  
  if(identical(sub_samp, check_for_default)){
    print("Defaulting max sub sample size to 10.")
    sub_samp <- 10
  }else{
    sub_samp <- as.numeric(options[sub_samp + 1])
  }
  
  if(identical(which(run_parameters$parameter == "Subsample Size"), check_for_default)){
    run_parameters <- rbind(run_parameters, data.frame(parameter = "Subsample Size", setting = sub_samp))
  }else{
    run_parameters$setting[which(run_parameters$parameter == "Subsample Size")] <- sub_samp
  }
  fwrite(run_parameters, "metapop_run_settings/run_settings.tsv", sep = "\t")
  
  
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
  
  run_summaries <- fread("metapop_run_settings/run_dates.tsv", sep = "\t")
  
  run_summaries$begin[run_summaries$stage=="Microdiversity"] <- as.character(Sys.time())
  
  fwrite(run_summaries, "metapop_run_settings/run_dates.tsv", sep = "\t")
  
  genes <- fread("metapop_codon_bias/codon_bias_sequences.tsv", sep = "\t", header = F, select = c(3,4))
  colnames(genes) = c("V1", "V2")
  
  #Codon bias
  if(dir.exists("metapop_codon_bias") & !file.exists("metapop_codon_bias/gene_IQR_and_mean.tsv")){
    
    
    
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
  
  #subsample and microdiversity
  # Occasionally the subsample will take 1 more than the sub_samp value indicates it should. 
  # This happens because the particular counts and scaling factor sample_prop round up more than once
  
  {
    # The normal data
    microdiv_data <- fread("metapop_called_snps/called_snps.tsv", sep = "\t", header = T)
    
    setkeyv(microdiv_data, c("source", "contig"))
    
    non_genic_data <- fread("metapop_called_snps/non_genic_snps.tsv", sep = "\t", header = T)
    
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
    assembled_contigs <- readDNAStringSet(original_assembly)
    assembled_contigs <- data.table(contig = names(assembled_contigs), length = nchar(assembled_contigs))
    
    assembled_contigs$contig <- unlist(lapply(assembled_contigs$contig, function(y){return(strsplit(y, split = " ")[[1]][1])}))
    
    gene_assembly <- fread("metapop_codon_bias/codon_bias_sequences.tsv", sep = "\t")
    colnames(gene_assembly) = c("seq", "OC", "gene", "contig")
    gene_assembly$length <- nchar(gene_assembly$seq)
    
    #Likely needs to be separate
    depth_file_names <- list.files(path = "metapop_cov_and_depth", full.names = T)
    
    #returns a 0 length data.table if there are no passing contigs
    depth_info <- lapply(depth_file_names, function(x){
      tmp <- fread(x, sep = "\t")
      tmp$source <- substr(x, 23, nchar(x)-28)
      tmp <- tmp[V3 >= cov_cutoff & V4 >= depth_cutoff,]
      return(tmp)
    })
    
    names(depth_info) = substr(depth_file_names, 23, nchar(depth_file_names)-28)
    
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
      
      #Fixed up to linked. Fantastic, time to do THIS...
      #If linked SNPs was done, corrects the microdiv_data for the directly observed read results
      if(dir.exists("metapop_linked_snps")){
        print("Correcting subsamples for linked snps")
        
        mined_reads <- fread("metapop_linked_snps/linked_snp_results.tsv", sep = "\t")
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
        print("We recommend running MetaPop_Mine_Reads.R before the microdiversity module to correct SNP data over reads.")
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
      
      
      if(!dir.exists("metapop_microdiversity")){
        system("mkdir metapop_microdiversity")
      }
      
      fwrite(microdiv_data, "metapop_microdiversity/global_raw_microdiversity_data_snp_loci_only.tsv", sep = "\t")
      
      fwrite(gene_level_microdiv, "metapop_microdiversity/global_gene_microdiversity.tsv", sep = "\t")
      fwrite(contig_level_microdiv, "metapop_microdiversity/global_contig_microdiversity.tsv", sep = "\t")
      
      fwrite(codon_pos_count, "metapop_microdiversity/global_codon_position_summary.tsv", sep = "\t")
      
    }
    
    if(sum(codon_pos_count$third_pos) < sum(codon_pos_count$second_pos) | sum(codon_pos_count$third_pos) < sum(codon_pos_count$first_pos)){
      print("More SNPs were observed in either the first or second positions of codons than in the third position. This is suspicious, and results should be treated cautiously.")
      flush.console()
    }
    
    
    #Local level stuff
    vcfs <- list.files(path = "metapop_variants", full.names = T)
    
    cl <- makeCluster(min(threads, detectCores(), length(vcfs)))
    registerDoParallel(cl)
    clusterExport(cl, c("library_location", "vcfs"), envir=environment())
    clusterEvalQ(cl=cl, expr = library(data.table, lib.loc = library_location))
    
    # Read and process mpileup files in parallel
    vcf_sums <- foreach(i=vcfs) %dopar%{
      tmp <- fread(i, sep = "\t", col.names = c("contig", "pos", "ref", "alt"))
      tmp[, source := substr(i, 18, nchar(i)-13)]
      tmp[, key := paste0(contig,"_", pos, source)]
      
      return(tmp)
    }
    
    stopCluster(cl)
    
    vcf_sums <- rbindlist(vcf_sums)
    
    #This should just be microdviersity_data in the real one, no fread needed
    #microdiversity_data <- fread("metapop_microdiversity/raw_microdiversity_data_snp_loci_only.tsv")
    
    
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
      if(dir.exists("metapop_linked_snps")){
        print("Correcting local subsamples for linked snps")
        
        microdiv_data_local$match_key <- paste0(microdiv_data_local$source, microdiv_data_local$contig_gene, microdiv_data_local$codon)
        mined_reads$match_key <- paste0(mined_reads$source, mined_reads$contig_gene, mined_reads$codon)
        
        mined_reads <- mined_reads[mined_reads$match_key %in% microdiv_data_local$match_key,]
        
        selector <- microdiv_data_local$match_key %in% mined_reads$match_key
        
        #sets counts on each codon for microdiv data to zero to prevent double counting; sets codons to blank for correcting
        
        microdiv_data_local[, obsN := ifelse(selector, 0, obsN)]
        microdiv_data_local[, obsS := ifelse(selector, 0, obsS)]
        
        microdiv_data_local[, affected_by_linked_snps := F]
        microdiv_data_local[, affected_by_linked_snps := ifelse(selector, T, F)]
        
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
      
      fwrite(microdiv_data_local, "metapop_microdiversity/local_raw_microdiversity_data_snp_loci_only.tsv", sep = "\t")
      
      fwrite(gene_level_microdiv, "metapop_microdiversity/local_gene_microdiversity.tsv", sep = "\t")
      fwrite(contig_level_microdiv, "metapop_microdiversity/local_contig_microdiversity.tsv", sep = "\t")
      
      fwrite(codon_pos_count, "metapop_microdiversity/local_codon_position_summary.tsv", sep = "\t")
      
    }
    
  }
  
  #FST
  
  {
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
    
    fwrite(fst_results, "metapop_microdiversity/fixation_index.tsv", sep = "\t")
    
  }
  
  run_summaries$end[run_summaries$stage=="Microdiversity"] <- as.character(Sys.time())
  
  fwrite(run_summaries, "metapop_run_settings/run_dates.tsv", sep = "\t")
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}

#Macrodiversity
metapop_macrodiv <- function(options){
  
  check_for_default <- integer(0)
  
  
  run_parameters <- read.csv("metapop_run_settings/run_settings.tsv", sep = "\t")
  run_parameters$parameter <- as.character(run_parameters$parameter)
  run_parameters$setting <- as.character(run_parameters$setting)
  
  
  library_location <- run_parameters$setting[run_parameters$parameter == "Library Location"]
  threads <- as.numeric(run_parameters$setting[run_parameters$parameter == "Threads"])
  
  
  run_summaries <- fread("metapop_run_settings/run_dates.tsv", sep = "\t")
  
  run_summaries$begin[run_summaries$stage=="Macrodiversity"] <- as.character(Sys.time())
  
  fwrite(run_summaries, "metapop_run_settings/run_dates.tsv", sep = "\t")
  
  counts <- which(grepl("-ct", options))
  
  if(identical(counts, check_for_default)){
    print("Macrodiversity requires a file with read or base pair counts per sample. Please specify this file using -ct")
    quit(save = "no")
  }else{
    counts <- options[counts + 1]
  }
  
  whole_mag <- which(grepl("-complete_bact", options))
  
  mag_cutoff <- which(grepl("-min_det", options))
  
  if(identical(whole_mag, check_for_default)){
    print("Macrodiversity will assume that each contig is a species. Set the -complete_bact to change this")
    whole_mag <- F
  }else{
    print("Macrodiversity will assume that each contig is a complete, closed bacterial genome.")
    whole_mag <- T
  }
  
  if(whole_mag){
    if(identical(mag_cutoff, check_for_default)){
      print("MAGs will be considered detected if 20% of the genome is covered. Set this with -min_det for a different value.")
      mag_cutoff <- 20
    }else{
      mag_cutoff <- as.numeric(options[mag_cutoff + 1])
    }
  }
  
  
  
  if(!whole_mag){
    min_contig_length <- which(grepl("-min_bp", options))
    
    #Looks complicated, but this code checks for user input on min_bp_detect and uses that if found. 
    #Otherwise checks for a previously entered value in the run params and uses that instead. Updates this value if it's found and the user entered something
    #Finally, if the user didn't enter a value, and the params file doesn't have it either, defaults to 5000 and adds it to the params file.
    #Note - this parameter updating code will be used sparingly. Most of the time, values should ABSOLUTELY NOT be changed from an original preprocessing run, as later steps mostly depend on these specific values.
    if(sum(run_parameters$parameter=="BP for Detect") == 1 & identical(min_contig_length, check_for_default)){
      min_contig_length <- as.numeric(run_parameters$setting[run_parameters$parameter=="BP for Detect"])
      contig_len <- data.frame(parameter = "BP for Detect", setting=min_contig_length)
    }else{
      
      if(identical(min_contig_length, check_for_default)){
        print("Setting number of bp for contigs to be considered detected to 5 Kbp. Use -min_bp [INT # of bp] to set a different value.")
        min_contig_length <- 5e3
      }else{
        min_contig_length <- as.numeric(options[min_contig_length + 1])
      }
      
      if(sum(run_parameters$parameter=="BP for Detect") == 1){
        run_parameters$setting[run_parameters$parameter=="BP for Detect"] <- min_contig_length
      }else{
        contig_len <- data.frame(parameter = "BP for Detect", setting=min_contig_length)
      }
      
      contig_len <- data.frame(parameter = "BP for Detect", setting=min_contig_length)
      
    }
    
    coverage_cutoff <- as.numeric(run_parameters$setting[run_parameters$parameter=="Coverage"])
    depth_cutoff <- as.numeric(run_parameters$setting[run_parameters$parameter=="Depth"])  
  }else{
    coverage_cutoff <- mag_cutoff
    depth_cutoff <- 0
    
    contig_len <- data.frame(parameter = c("Whole Bacterial Genomes", "Bacterial Detection Threshold"), setting=c(TRUE, mag_cutoff))
    run_parameters <- rbind(run_parameters, contig_len)
  }
  
  if(!dir.exists("metapop_macrodiversity")){
    system("mkdir metapop_macrodiversity")
  }
  
  norm_file <- data.frame(parameter = c("Read or BP Counts File"), setting=c(counts))
  run_parameters <- rbind(run_parameters, contig_len, norm_file)
  fwrite(run_parameters, "metapop_run_settings/run_settings.tsv", sep = "\t")
  
  cov_depth <- list.files(full.names = T, path = "metapop_cov_and_depth")
  
  norm_file <- counts
  norm_file <- fread(norm_file, sep = "\t")
  
  original_assembly <- as.character(run_parameters$setting[run_parameters$parameter=="Assembly"])
  geneFile <- as.character(run_parameters$setting[run_parameters$parameter=="Genes"])
  
  original_assemblies <- readDNAStringSet(original_assembly)
  
  #Hopefully nothing breaks here
  names(original_assemblies) <- unlist(lapply(names(original_assemblies), function(x){
    return(strsplit(x, split = " ")[[1]][1])
  }))
  
  #Long format is generally better for operations like this. It's not a matrix yet, but you can ignore that
  abundance_matrix <- CJ(names(original_assemblies), unlist(lapply(cov_depth, function(x){
    
    strsplit(strsplit(x, split = "metapop_cov_and_depth/")[[1]][2], split = ".cov_and_depth_by_contig.tsv")[[1]][1]
    
  })))
  
  #Initialize abundance counts at zero.
  abundance_matrix$abundance = 0
  
  #Figures out the normalization factors on each sample.
  norm_file$factor <- norm_file$V2/max(norm_file$V2)
  
  #This both creates and normalizes the abundance table, with the caveats that only things with >= min. coverage or >= min bp covered are considered with depth >= 1
  for(i in cov_depth){
    
    match_name <- strsplit(strsplit(i, split = "metapop_cov_and_depth/")[[1]][2], split = ".cov_and_depth_by_contig.tsv")[[1]][1]
    
    tmp <- fread(i, sep = "\t")
    
    tmp$V4[!(tmp$V3 >= coverage_cutoff | (tmp$V2*(tmp$V3/100))>=min_contig_length)] <- 0
    
    #Needs the raw table too
    #tmp$V4 <- tmp$V4/norm_file$factor[match(match_name, norm_file$V1)]
    
    #This looks odd, but gets the rows which match the sample, then matches the contigs from the depths and adds in the values to the correct rows.
    abundance_matrix[V2==match_name,][match(tmp$V1, V1)]$abundance <- tmp$V4
    
  }
  
  
  
  #Transform the abundance matrix into an actual matrix format. Rows are contigs, columns are samples.
  raw_abundance_matrix <- dcast(abundance_matrix, V1~V2, value.var = "abundance")
  
  
  #Goes ahead and outputs the table. This is what a person looking to send this table would want to send.
  fwrite(raw_abundance_matrix, "metapop_macrodiversity/raw_abundances_table.tsv", sep = "\t")
  
  rm(raw_abundance_matrix)
  
  #normalize the abundance matrix
  abundance_matrix[, abundance := (abundance / norm_file$factor[match(V2, norm_file$V1)] ) ]
  
  abundance_matrix <- dcast(abundance_matrix, V1~V2, value.var = "abundance")
  rn = abundance_matrix$V1
  
  fwrite(abundance_matrix, "metapop_macrodiversity/normalized_abundances_table.tsv", sep = "\t")
  
  #Technically the first column was row names. Removes it.
  abundance_matrix[,V1:=NULL]
  
  #Coverts to the R representation of a matrix instead of a data frame.
  abundance_matrix <- as.matrix(abundance_matrix)
  
  #Names the rows with their contig name.
  rownames(abundance_matrix) = rn
  
  if(!dir.exists("metapop_visualizations")){
    dir.create("metapop_visualizations")
  }
  
  #Determine distibution of values in matrix
  data_distribution_matrix <- abundance_matrix
  data_distribution_matrix <- ifelse(data_distribution_matrix==0, NA, data_distribution_matrix)
  top_75_quantile <- quantile(data_distribution_matrix, probs = c(.75), na.rm = TRUE)
  
  #Create heatmap pdf of normalized abundances
  pheatmap_matrix <- abundance_matrix 
  pheatmap_matrix[pheatmap_matrix > top_75_quantile[[1]]] <- top_75_quantile[[1]]
  
  #Calculate alpha-diversity metrics
  Richness<-specnumber(t(abundance_matrix)) #calculates richness
  Shannons_H<- diversity(t(abundance_matrix), index = "shannon") #calculates shannon index
  Simpson<- diversity(t(abundance_matrix), index = "simpson") #calculates simpson index
  InvSimpson<- diversity(t(abundance_matrix), index = "invsimpson") #calculates inverse simpson index
  Peilous_J <- Shannons_H/log(Richness) #calculates peilou's J
  mode(abundance_matrix) <- "integer"
  Fisher <- fisher.alpha(t(abundance_matrix))
  Estimation_data <- estimateR(t(abundance_matrix))
  Estimation_data <- t(as.data.frame(Estimation_data))
  colnames(Estimation_data) <- c("Observed","Chao1","se.Chao1","ACE","se.ACE")
  alpha_diversity <- cbind(Estimation_data, Shannons_H, Simpson, InvSimpson, Fisher, Peilous_J)
  write.table(alpha_diversity, file = "metapop_macrodiversity/Alpha_diveristy_stats.tsv", sep = "\t", row.names = TRUE, col.names=NA)
  
  #Make plots of all alpha-diversity metrics with line showing the mean and median
  #Make the pdfs output singletons and the 4X page for ezch
  
  alpha_diversity<-as.data.frame(alpha_diversity)
  alpha_diversity_number_samples <- nrow(alpha_diversity)
  alpha_diversity_names <- rownames(alpha_diversity)
  alpha_diversity_data <- as.data.frame(cbind(alpha_diversity_names, alpha_diversity))
  
  Richness_plot <- ggplot(data=alpha_diversity_data)+
    geom_point(aes(x=alpha_diversity_data$alpha_diversity_names, y=alpha_diversity_data$Observed), fill="grey", color = "black",shape=21, size = 2)+
    geom_hline(aes(yintercept=median(alpha_diversity_data$Observed, na.rm=TRUE)),linetype="dotted",color="red")+
    geom_hline(aes(yintercept=mean(alpha_diversity_data$Observed, na.rm=TRUE)),color="red")+
    #change x = 4 to max(num_samples)*.25
    annotate("text", x=alpha_diversity_number_samples/4, y=median(alpha_diversity_data$Observed, na.rm=TRUE) ,label = "Median", color = "red", size = 4, vjust = 1)+
    annotate("text", x=alpha_diversity_number_samples*3/4, y=mean(alpha_diversity_data$Observed, na.rm=TRUE) ,label = "Mean", color = "red", size = 4, vjust = 0)+
    xlab("Samples")+
    scale_y_continuous(name ="Species\nRichness")+
    theme_bw()+
    theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, size = 6),axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  
  Chao1_plot <- ggplot(data=alpha_diversity_data)+
    geom_point(aes(x=alpha_diversity_data$alpha_diversity_names, y=alpha_diversity_data$Chao1), fill="grey", color = "black",shape=21, size = 2)+
    geom_hline(aes(yintercept=median(alpha_diversity_data$Chao1, na.rm=TRUE)),linetype="dotted",color="red")+
    geom_hline(aes(yintercept=mean(alpha_diversity_data$Chao1, na.rm=TRUE)),color="red")+
    annotate("text", x=alpha_diversity_number_samples/4, y=median(alpha_diversity_data$Chao1, na.rm=TRUE) ,label = "Median", color = "red", size = 4, vjust = 1)+
    annotate("text", x=alpha_diversity_number_samples*3/4, y=mean(alpha_diversity_data$Chao1, na.rm=TRUE) ,label = "Mean", color = "red", size = 4, vjust = 0)+
    xlab("Samples")+
    scale_y_continuous(name ="Chao1")+
    theme_bw()+
    theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, size = 6),axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  
  
  ACE_plot <- ggplot(data=alpha_diversity_data)+
    geom_point(aes(x=alpha_diversity_data$alpha_diversity_names, y=alpha_diversity_data$ACE), fill="grey", color = "black",shape=21, size = 2)+
    geom_hline(aes(yintercept=median(alpha_diversity_data$ACE, na.rm=TRUE)),linetype="dotted",color="red")+
    geom_hline(aes(yintercept=mean(alpha_diversity_data$ACE, na.rm=TRUE)),color="red")+
    annotate("text", x=alpha_diversity_number_samples/4,y=median(alpha_diversity_data$ACE, na.rm=TRUE) ,label = "Median", color = "red", size = 4, vjust = 1)+
    annotate("text", x=alpha_diversity_number_samples*3/4,y=mean(alpha_diversity_data$ACE, na.rm=TRUE) ,label = "Mean", color = "red", size = 4, vjust = 0)+
    xlab("Samples")+
    scale_y_continuous(name ="ACE")+
    theme_bw()+
    theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, size = 6),axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  
  
  Shannons_H_plot <- ggplot(data=alpha_diversity_data)+
    geom_point(aes(x=alpha_diversity_data$alpha_diversity_names, y=alpha_diversity_data$Shannons_H), fill="grey", color = "black",shape=21, size = 2)+
    geom_hline(aes(yintercept=median(alpha_diversity_data$Shannons_H, na.rm=TRUE)),linetype="dotted",color="red")+
    geom_hline(aes(yintercept=mean(alpha_diversity_data$Shannons_H, na.rm=TRUE)),color="red")+
    annotate("text", x=alpha_diversity_number_samples/4,y=median(alpha_diversity_data$Shannons_H, na.rm=TRUE) ,label = "Median", color = "red", size = 4, vjust = 1)+
    annotate("text", x=alpha_diversity_number_samples*3/4,y=mean(alpha_diversity_data$Shannons_H, na.rm=TRUE) ,label = "Mean", color = "red", size = 4, vjust = 0)+
    xlab("Samples")+
    scale_y_continuous(name ="Shannon's H")+
    theme_bw()+
    theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, size = 6),axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  
  
  Simpson_plot <- ggplot(data=alpha_diversity_data)+
    geom_point(aes(x=alpha_diversity_data$alpha_diversity_names, y=alpha_diversity_data$Simpson), fill="grey", color = "black",shape=21, size = 2)+
    geom_hline(aes(yintercept=median(alpha_diversity_data$Simpson, na.rm=TRUE)),linetype="dotted",color="red")+
    geom_hline(aes(yintercept=mean(alpha_diversity_data$Simpson, na.rm=TRUE)),color="red")+
    annotate("text", x=alpha_diversity_number_samples/4,y=median(alpha_diversity_data$Simpson, na.rm=TRUE) ,label = "Median", color = "red", size = 4, vjust = 1)+
    annotate("text", x=alpha_diversity_number_samples*3/4,y=mean(alpha_diversity_data$Simpson, na.rm=TRUE) ,label = "Mean", color = "red", size = 4, vjust = 0)+
    xlab("Samples")+
    ylab("Simpson")+
    theme_bw()+
    theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, size = 6),axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  
  
  InvSimpson_plot <- ggplot(data=alpha_diversity_data)+
    geom_point(aes(x=alpha_diversity_data$alpha_diversity_names, y=alpha_diversity_data$InvSimpson), fill="grey", color = "black",shape=21, size = 2)+
    geom_hline(aes(yintercept=median(alpha_diversity_data$InvSimpson, na.rm=TRUE)),linetype="dotted",color="red")+
    geom_hline(aes(yintercept=mean(alpha_diversity_data$InvSimpson, na.rm=TRUE)),color="red")+
    annotate("text", x=alpha_diversity_number_samples/4,y=median(alpha_diversity_data$InvSimpson, na.rm=TRUE) ,label = "Median", color = "red", size = 4, vjust = 1)+
    annotate("text", x=alpha_diversity_number_samples*3/4,y=mean(alpha_diversity_data$InvSimpson, na.rm=TRUE) ,label = "Mean", color = "red", size = 4, vjust = 0)+
    xlab("Samples")+
    scale_y_continuous(name ="Inverse Simpson")+
    theme_bw()+
    theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, size = 6),axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  
  Fisher_plot <- ggplot(data=alpha_diversity_data)+
    geom_point(aes(x=alpha_diversity_data$alpha_diversity_names, y=alpha_diversity_data$Fisher), fill="grey", color = "black",shape=21, size = 2)+
    geom_hline(aes(yintercept=median(alpha_diversity_data$Fisher, na.rm=TRUE)),linetype="dotted",color="red")+
    geom_hline(aes(yintercept=mean(alpha_diversity_data$Fisher, na.rm=TRUE)),color="red")+
    annotate("text", x=alpha_diversity_number_samples/4,y=median(alpha_diversity_data$Fisher, na.rm=TRUE) ,label = "Median", color = "red", size = 4, vjust = 1)+
    annotate("text", x=alpha_diversity_number_samples*3/4,y=mean(alpha_diversity_data$Fisher, na.rm=TRUE) ,label = "Mean", color = "red", size = 4, vjust = 0)+
    xlab("Samples")+
    scale_y_continuous(name ="Fisher")+
    theme_bw()+
    theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, size = 6),axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  
  
  Peilous_J_plot <- ggplot(data=alpha_diversity_data)+
    geom_point(aes(x=alpha_diversity_data$alpha_diversity_names, y=alpha_diversity_data$Peilous_J), fill="grey", color = "black",shape=21, size = 2)+
    geom_hline(aes(yintercept=median(alpha_diversity_data$Peilous_J, na.rm=TRUE)),linetype="dotted",color="red")+
    geom_hline(aes(yintercept=mean(alpha_diversity_data$Peilous_J, na.rm=TRUE)),color="red")+
    annotate("text", x=alpha_diversity_number_samples/4,y=median(alpha_diversity_data$Peilous_J, na.rm=TRUE) ,label = "Median", color = "red", size = 4, vjust = 1)+
    annotate("text", x=alpha_diversity_number_samples*3/4,y=mean(alpha_diversity_data$Peilous_J, na.rm=TRUE) ,label = "Mean", color = "red", size = 4, vjust = 0)+
    xlab("Samples")+
    ylab("Peilous J")+
    theme_bw()+
    theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 90, size = 6),axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  
  #Make pdf of alpha diversity scatterplots
  pdf("metapop_visualizations/Alpha_diversity_scatterplots.pdf", width = 16, height = 16/3)
  
  print(Richness_plot)
  print(Shannons_H_plot)
  print(Simpson_plot)
  print(Peilous_J_plot)
  
  print(Chao1_plot)
  print(ACE_plot)
  print(InvSimpson_plot)
  print(Fisher_plot)
  
  
  dev.off()
  
  
  #Calculate beta-diversity metrics
  clr_transformation <- clr(t(abundance_matrix))
  clr_euc <-vegdist(clr_transformation, method="euclidean") 
  clr_euc <- as.matrix(clr_euc)
  write.table(clr_euc, file = "metapop_macrodiversity/Beta_diveristy_clr_transformed_euclidean_distances.tsv", sep = "\t", row.names = TRUE, col.names=NA)
  bray <-vegdist(t(abundance_matrix), method="bray") 
  bray <- as.matrix(bray)
  write.table(bray, file = "metapop_macrodiversity/Beta_diveristy_bray_distances.tsv", sep = "\t", row.names = TRUE, col.names=NA)
  jaccard <-vegdist(t(abundance_matrix), method="jaccard") 
  jaccard <- as.matrix(jaccard)
  write.table(jaccard, file = "metapop_macrodiversity/Beta_diveristy_jaccard_distances.tsv", sep = "\t", row.names = TRUE, col.names=NA)
  
  if(ncol(clr_euc) < 3){
    print("There are not enough samples to produce a meaningful PCA, PCoA, or NMDS analysis. At least 3 are required.")
    print("There will be data outputs for these analyses in metapop_macrodiversity/, but no visualizations.")
    print("As these outputs are distances, there is not much meaning in the data output.")
  }else{
    
    #Make PCA plots of all beta-diversity metrics
    clr_euc.pca<-prcomp(clr_euc, scale = TRUE)
    clr_euc.eig2 <- eigenvals(clr_euc.pca)
    clr_euc.percentage_variance_explained <- clr_euc.eig2 / sum(clr_euc.eig2)
    clr_euc.PC1_percent <- as.numeric(format(round((clr_euc.percentage_variance_explained[[1]]*100), 2), nsmall = 2))
    clr_euc.PC2_percent <- as.numeric(format(round((clr_euc.percentage_variance_explained[[2]]*100), 2), nsmall = 2))
    clr_euc.PC1_percent= sprintf("%.2f%%", clr_euc.PC1_percent)
    clr_euc.PC2_percent= sprintf("%.2f%%", clr_euc.PC2_percent)
    clr_euc.PC1_percent= paste ("PC1 (",clr_euc.PC1_percent,")")
    clr_euc.PC2_percent= paste ("PC2 (",clr_euc.PC2_percent,")")
    clr_euc.components <- as.data.frame(clr_euc.pca$x)
    clr_euc.components.richness <- as.data.frame(cbind(clr_euc.components, alpha_diversity)) ##add alpha diversity values
    
    CLR_EUC_PCA_PLOT <- ggplot(data=clr_euc.components.richness)+
      geom_point(aes(x=clr_euc.components.richness$PC1, y=clr_euc.components.richness$PC2, fill=clr_euc.components.richness$Observed), color = "black",shape=21, size = 4)+
      scale_fill_gradientn(colours = terrain.colors(10))+
      ggtitle("Centered-log Ratio Transformed Euclidean Distances")+
      labs(fill = "Species Richness")+xlab(clr_euc.PC1_percent)+
      ylab(clr_euc.PC2_percent)+
      theme_bw()+
      xlim(c(min(c(clr_euc.components.richness$PC1, clr_euc.components.richness$PC2), na.rm=T), max(c(clr_euc.components.richness$PC1, clr_euc.components.richness$PC2), na.rm = T))) +
      ylim(c(min(c(clr_euc.components.richness$PC1, clr_euc.components.richness$PC2), na.rm=T), max(c(clr_euc.components.richness$PC1, clr_euc.components.richness$PC2), na.rm = T))) +
      coord_fixed() +
      theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    
    bray.pca<-prcomp(bray, scale = TRUE)
    bray.eig2 <- eigenvals(bray.pca)
    bray.percentage_variance_explained <- bray.eig2 / sum(bray.eig2)
    bray.PC1_percent <- as.numeric(format(round((bray.percentage_variance_explained[[1]]*100), 2), nsmall = 2))
    bray.PC2_percent <- as.numeric(format(round((bray.percentage_variance_explained[[2]]*100), 2), nsmall = 2))
    bray.PC1_percent= sprintf("%.2f%%", bray.PC1_percent)
    bray.PC2_percent= sprintf("%.2f%%", bray.PC2_percent)
    bray.PC1_percent= paste ("PC1 (",bray.PC1_percent,")")
    bray.PC2_percent= paste ("PC2 (",bray.PC2_percent,")")
    bray.components <- as.data.frame(bray.pca$x)
    bray.components.richness <- as.data.frame(cbind(bray.components, alpha_diversity)) ##add alpha diversity values
    
    BRAY_PCA_PLOT <- ggplot(data=bray.components.richness)+
      geom_point(aes(x=bray.components.richness$PC1, y=bray.components.richness$PC2, fill=bray.components.richness$Observed), color = "black",shape=21, size = 4)+
      scale_fill_gradientn(colours = terrain.colors(10))+
      ggtitle("Bray-Curtis (Dissimilarity) Distances")+
      labs(fill = "Species Richness")+xlab(bray.PC1_percent)+
      ylab(bray.PC2_percent)+
      theme_bw()+
      xlim(c(min(c(bray.components.richness$PC1, bray.components.richness$PC2), na.rm=T), max(c(bray.components.richness$PC1, bray.components.richness$PC2), na.rm = T))) +
      ylim(c(min(c(bray.components.richness$PC1, bray.components.richness$PC2), na.rm=T), max(c(bray.components.richness$PC1, bray.components.richness$PC2), na.rm = T))) +
      coord_fixed() +
      theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    
    
    jaccard.pca<-prcomp(jaccard, scale = TRUE)
    jaccard.eig2 <- eigenvals(jaccard.pca)
    jaccard.percentage_variance_explained <- jaccard.eig2 / sum(jaccard.eig2)
    jaccard.PC1_percent <- as.numeric(format(round((jaccard.percentage_variance_explained[[1]]*100), 2), nsmall = 2))
    jaccard.PC2_percent <- as.numeric(format(round((jaccard.percentage_variance_explained[[2]]*100), 2), nsmall = 2))
    jaccard.PC1_percent= sprintf("%.2f%%", jaccard.PC1_percent)
    jaccard.PC2_percent= sprintf("%.2f%%", jaccard.PC2_percent)
    jaccard.PC1_percent= paste ("PC1 (",jaccard.PC1_percent,")")
    jaccard.PC2_percent= paste ("PC2 (",jaccard.PC2_percent,")")
    jaccard.components <- as.data.frame(jaccard.pca$x)
    jaccard.components.richness <- as.data.frame(cbind(jaccard.components, alpha_diversity)) ##add alpha diversity values
    
    JACCARD_PCA_PLOT <- ggplot(data=jaccard.components.richness)+
      geom_point(aes(x=jaccard.components.richness$PC1, y=jaccard.components.richness$PC2, fill=jaccard.components.richness$Observed), color = "black",shape=21, size = 4)+
      scale_fill_gradientn(colours = terrain.colors(10))+
      ggtitle("Jaccard (Similarity) Distances")+
      labs(fill = "Species Richness")+
      xlab(jaccard.PC1_percent)+ylab(jaccard.PC2_percent)+
      theme_bw()+
      xlim(c(min(c(jaccard.components.richness$PC1, jaccard.components.richness$PC2), na.rm=T), max(c(jaccard.components.richness$PC1, jaccard.components.richness$PC2), na.rm = T))) +
      ylim(c(min(c(jaccard.components.richness$PC1, jaccard.components.richness$PC2), na.rm=T), max(c(jaccard.components.richness$PC1, jaccard.components.richness$PC2), na.rm = T))) +
      coord_fixed() +
      theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    
    #Make pdf of PCA plots
    pdf("metapop_visualizations/PCA_CLR.EUC_BRAY_JACCARD_plot.pdf", width = 9, height = 9)
    print(CLR_EUC_PCA_PLOT)
    print(BRAY_PCA_PLOT)
    print(JACCARD_PCA_PLOT)
    dev.off()
    
    #Make PCOA plots of all beta-diversity metrics
    clr_euc.pcoa<-capscale(clr_euc~-1)
    clr_euc.eig2 <- eigenvals(clr_euc.pcoa)
    clr_euc.percentage_variance_explained <- clr_euc.eig2 / sum(clr_euc.eig2)
    clr_euc.PC1_percent <- as.numeric(format(round((clr_euc.percentage_variance_explained[[1]]*100), 2), nsmall = 2))
    clr_euc.PC2_percent <- as.numeric(format(round((clr_euc.percentage_variance_explained[[2]]*100), 2), nsmall = 2))
    clr_euc.PC1_percent= sprintf("%.2f%%", clr_euc.PC1_percent)
    clr_euc.PC2_percent= sprintf("%.2f%%", clr_euc.PC2_percent)
    clr_euc.PC1_percent= paste ("PCo1 (",clr_euc.PC1_percent,")")
    clr_euc.PC2_percent= paste ("PCo2 (",clr_euc.PC2_percent,")")
    clr_euc.components <- as.data.frame(scores(clr_euc.pcoa, display=c("sites")))
    clr_euc.components.richness <- as.data.frame(cbind(clr_euc.components, alpha_diversity)) ##add alpha diversity values
    
    CLR_EUC_PCOA_PLOT <- ggplot(data=clr_euc.components.richness)+
      geom_point(aes(x=clr_euc.components.richness$MDS1, y=clr_euc.components.richness$MDS2, fill=clr_euc.components.richness$Observed), color = "black",shape=21, size = 4)+
      scale_fill_gradientn(colours = terrain.colors(10))+
      ggtitle("Centered-log Ratio Transformed Euclidean Distances")+
      labs(fill = "Species Richness")+
      xlab(clr_euc.PC1_percent)+
      ylab(clr_euc.PC2_percent)+
      theme_bw()+
      xlim(c(min(c(clr_euc.components.richness$MDS1, clr_euc.components.richness$MDS2), na.rm=T), max(c(clr_euc.components.richness$MDS1, clr_euc.components.richness$MDS2), na.rm = T))) +
      ylim(c(min(c(clr_euc.components.richness$MDS1, clr_euc.components.richness$MDS2), na.rm=T), max(c(clr_euc.components.richness$MDS1, clr_euc.components.richness$MDS2), na.rm = T))) +
      coord_fixed() +
      theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    
    
    bray.pcoa<-capscale(bray~-1)
    bray.eig2 <- eigenvals(bray.pcoa)
    bray.percentage_variance_explained <- bray.eig2 / sum(bray.eig2)
    bray.PC1_percent <- as.numeric(format(round((bray.percentage_variance_explained[[1]]*100), 2), nsmall = 2))
    bray.PC2_percent <- as.numeric(format(round((bray.percentage_variance_explained[[2]]*100), 2), nsmall = 2))
    bray.PC1_percent= sprintf("%.2f%%", bray.PC1_percent)
    bray.PC2_percent= sprintf("%.2f%%", bray.PC2_percent)
    bray.PC1_percent= paste ("PCo1 (",bray.PC1_percent,")")
    bray.PC2_percent= paste ("PCo2 (",bray.PC2_percent,")")
    bray.components <- as.data.frame(scores(bray.pcoa, display=c("sites")))
    bray.components.richness <- as.data.frame(cbind(bray.components, alpha_diversity)) ##add alpha diversity values
    
    BRAY_PCOA_PLOT <- ggplot(data=bray.components.richness)+
      geom_point(aes(x=bray.components.richness$MDS1, y=bray.components.richness$MDS2, fill=bray.components.richness$Observed), color = "black",shape=21, size = 4)+
      scale_fill_gradientn(colours = terrain.colors(10))+
      ggtitle("Bray-Curtis (Dissimilarity) Distances")+
      labs(fill = "Species Richness")+
      xlab(bray.PC1_percent)+
      ylab(bray.PC2_percent)+
      theme_bw()+
      xlim(c(min(c(bray.components.richness$MDS1, bray.components.richness$MDS2), na.rm=T), max(c(bray.components.richness$MDS1, bray.components.richness$MDS2), na.rm = T))) +
      ylim(c(min(c(bray.components.richness$MDS1, bray.components.richness$MDS2), na.rm=T), max(c(bray.components.richness$MDS1, bray.components.richness$MDS2), na.rm = T))) +
      coord_fixed() +
      theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    
    
    jaccard.pcoa<-capscale(jaccard~-1)
    jaccard.eig2 <- eigenvals(jaccard.pcoa)
    jaccard.percentage_variance_explained <- jaccard.eig2 / sum(jaccard.eig2)
    jaccard.PC1_percent <- as.numeric(format(round((jaccard.percentage_variance_explained[[1]]*100), 2), nsmall = 2))
    jaccard.PC2_percent <- as.numeric(format(round((jaccard.percentage_variance_explained[[2]]*100), 2), nsmall = 2))
    jaccard.PC1_percent= sprintf("%.2f%%", jaccard.PC1_percent)
    jaccard.PC2_percent= sprintf("%.2f%%", jaccard.PC2_percent)
    jaccard.PC1_percent= paste ("PCo1 (",jaccard.PC1_percent,")")
    jaccard.PC2_percent= paste ("PCo2 (",jaccard.PC2_percent,")")
    jaccard.components <- as.data.frame(scores(jaccard.pcoa, display=c("sites")))
    jaccard.components.richness <- as.data.frame(cbind(jaccard.components, alpha_diversity)) ##add alpha diversity values
    
    JACCARD_PCOA_PLOT <- ggplot(data=jaccard.components.richness)+
      geom_point(aes(x=jaccard.components.richness$MDS1, y=jaccard.components.richness$MDS2, fill=jaccard.components.richness$Observed), color = "black",shape=21, size = 4)+
      scale_fill_gradientn(colours = terrain.colors(10))+
      ggtitle("Jaccard (Similarity) Distances")+
      labs(fill = "Species Richness")+
      xlab(jaccard.PC1_percent)+
      ylab(jaccard.PC2_percent)+
      theme_bw()+
      xlim(c(min(c(jaccard.components.richness$MDS1, jaccard.components.richness$MDS2), na.rm=T), max(c(jaccard.components.richness$MDS1, jaccard.components.richness$MDS2), na.rm = T))) +
      ylim(c(min(c(jaccard.components.richness$MDS1, jaccard.components.richness$MDS2), na.rm=T), max(c(jaccard.components.richness$MDS1, jaccard.components.richness$MDS2), na.rm = T))) +
      coord_fixed() +
      theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    
    #Make pdf of PCoA plots
    pdf("metapop_visualizations/PCoA_CLR.EUC_BRAY_JACCARD_plot.pdf", width = 9, height = 9)
    print(CLR_EUC_PCOA_PLOT)
    print(BRAY_PCOA_PLOT)
    print(JACCARD_PCOA_PLOT)
    dev.off()
    
    #Make NMDS plots of all beta-diversity metrics
    clr_euc.nmds<-metaMDS(clr_euc, k=2)
    clr_euc.stress <- clr_euc.nmds$stress
    clr_euc.stress= sprintf("%#.3f", clr_euc.stress)
    clr_euc.stress = paste ("stress = ",clr_euc.stress,"")
    clr_euc.points <- clr_euc.nmds$points 
    clr_euc.points.richness <- as.data.frame(cbind(clr_euc.points, alpha_diversity)) ##add alpha diversity values
    
    CLR_EUC_NMDS_PLOT <- ggplot(data=clr_euc.points.richness)+
      geom_point(aes(x=clr_euc.points.richness$MDS1, y=clr_euc.points.richness$MDS2, fill=clr_euc.components.richness$Observed), color = "black",shape=21, size = 4)+
      annotate(geom = 'text', label = clr_euc.stress, x = -Inf, y = Inf, hjust = 0, vjust = 1)+scale_fill_gradientn(colours = terrain.colors(10))+
      ggtitle("Centered-log Ratio Transformed Euclidean Distances")+
      labs(fill = "Species Richness")+xlab("NMDS1")+
      ylab("NMDS2")+
      theme_bw()+
      theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
      xlim(min(c(clr_euc.points.richness$MDS1, clr_euc.points.richness$MDS2), na.rm = T),  max(c(clr_euc.points.richness$MDS1, clr_euc.points.richness$MDS2), na.rm = T)) +
      ylim(min(c(clr_euc.points.richness$MDS1, clr_euc.points.richness$MDS2), na.rm = T),  max(c(clr_euc.points.richness$MDS1, clr_euc.points.richness$MDS2), na.rm = T)) +
      coord_fixed()
    
    
    
    bray.nmds<-metaMDS(bray, k=2)
    bray.stress <- bray.nmds$stress
    bray.stress= sprintf("%#.3f", bray.stress)
    bray.stress = paste ("stress = ",bray.stress,"")
    bray.points <- bray.nmds$points 
    bray.points.richness <- as.data.frame(cbind(bray.points, alpha_diversity)) ##add alpha diversity values
    
    BRAY_NMDS_PLOT <- ggplot(data=bray.points.richness)+
      geom_point(aes(x=bray.points.richness$MDS1, y=bray.points.richness$MDS2, fill=bray.components.richness$Observed), color = "black",shape=21, size = 4)+
      annotate(geom = 'text', label = bray.stress, x = -Inf, y = Inf, hjust = 0, vjust = 1)+
      scale_fill_gradientn(colours = terrain.colors(10))+
      ggtitle("Bray-Curtis (Dissimilarity) Distances")+
      labs(fill = "Species Richness")+
      xlab("NMDS1")+
      ylab("NMDS2")+
      theme_bw()+
      theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
      xlim(min(c(bray.points.richness$MDS1, bray.points.richness$MDS2), na.rm = T),  max(c(bray.points.richness$MDS1, bray.points.richness$MDS2), na.rm = T)) +
      ylim(min(c(bray.points.richness$MDS1, bray.points.richness$MDS2), na.rm = T),  max(c(bray.points.richness$MDS1, bray.points.richness$MDS2), na.rm = T)) +
      coord_fixed()
    
    
    jaccard.nmds<-metaMDS(jaccard, k=2)
    jaccard.stress <- jaccard.nmds$stress
    jaccard.stress= sprintf("%#.3f", jaccard.stress)
    jaccard.stress = paste ("stress = ",jaccard.stress,"")
    jaccard.points <- jaccard.nmds$points 
    jaccard.points.richness <- as.data.frame(cbind(jaccard.points, alpha_diversity)) ##add alpha diversity values
    
    JACCARD_NMDS_PLOT <- ggplot(data=jaccard.points.richness)+
      geom_point(aes(x=jaccard.points.richness$MDS1, y=jaccard.points.richness$MDS2, fill=jaccard.components.richness$Observed), color = "black",shape=21, size = 4)+
      annotate(geom = 'text', label = jaccard.stress, x = -Inf, y = Inf, hjust = 0, vjust = 1)+
      scale_fill_gradientn(colours = terrain.colors(10))+
      ggtitle("Jaccard (Similarity) Distances")+
      labs(fill = "Species Richness")+
      xlab("NMDS1")+
      ylab("NMDS2")+
      theme_bw()+
      theme(axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
      xlim(min(c(jaccard.points.richness$MDS1, jaccard.points.richness$MDS2), na.rm = T),  max(c(jaccard.points.richness$MDS1, jaccard.points.richness$MDS2), na.rm = T)) +
      ylim(min(c(jaccard.points.richness$MDS1, jaccard.points.richness$MDS2), na.rm = T),  max(c(jaccard.points.richness$MDS1, jaccard.points.richness$MDS2), na.rm = T)) +
      coord_fixed()
    
    
    #Make pdf of NMDS plots
    pdf("metapop_visualizations/NMDS_CLR.EUC_BRAY_JACCARD_plot.pdf", width = 9, height = 9)
    par(pty="s")
    print(CLR_EUC_NMDS_PLOT)
    print(BRAY_NMDS_PLOT)
    print(JACCARD_NMDS_PLOT)
    dev.off()
    
  }
  
  run_summaries$end[run_summaries$stage=="Macrodiversity"] <- as.character(Sys.time())
  
  fwrite(run_summaries, "metapop_run_settings/run_dates.tsv", sep = "\t")
  
  print("Heatmaps are very memory intensive. MetaPop will attempt to print normalized abundance heatmaps, but these may fail.")
  
  
  try({
    
    abundance_matrix_nozero <- abundance_matrix[ rowSums(abundance_matrix)!=0, ]
    
    pdf("metapop_visualizations/normalized_abundances_heatmap.pdf")
    
    pheatmap(abundance_matrix_nozero,
             scale="none",
             cluster_cols = T,
             cluster_rows = T,
             show_rownames = F, 
             show_colnames = T,
             fontsize_col = 5+1/log10(ncol(abundance_matrix)), 
             breaks=c(-0.00001,0,seq(0.0000001,max(abundance_matrix),length.out=50)), 
             border_color=NA, 
             col=c("white",colorRampPalette(brewer.pal(11,"Spectral"))(59)))
    
    
    dev.off()
    
  })
  
  
  
  
  try({
    
    pheatmap_matrix_nozero <- pheatmap_matrix[ rowSums(pheatmap_matrix)!=0, ]
    
    pdf("metapop_visualizations/normalized_abundances_heatmap_75quantile_removed.pdf")
    
    pheatmap(pheatmap_matrix_nozero,
             scale="none",
             cluster_cols = T,
             cluster_rows = T,
             show_rownames = F, 
             show_colnames = T,
             fontsize_col = 5+1/log10(ncol(pheatmap_matrix)), 
             breaks=c(-0.00001,0,seq(0.0000001,top_75_quantile[[1]],length.out=50)), 
             border_color=NA, 
             col=c("white",colorRampPalette(brewer.pal(11,"Spectral"))(59)))
    
    dev.off()
    
  })
  
  
  run_summaries$end[run_summaries$stage=="Macrodiversity"] <- as.character(Sys.time())
  
  fwrite(run_summaries, "metapop_run_settings/run_dates.tsv", sep = "\t")
}

#Micro viz
metapop_microdiv_viz <- function(options){
  
  check_for_default <- integer(0)

  
  plot_all <- which(grepl("-all", options))
  
  snp_scale <- which(grepl("-snp_scale", options))
  
  
  if(identical(check_for_default, plot_all)){
    print("Selecting the 3 contigs with the highest count of genes under selection from each sample, and plotting each appearance of these contigs across all samples.")
    print("Use flag -all to plot all contigs. Note: This may take a very long time.")
    plot_all <- F
  }else{
    plot_all <- T
  }
  
  
  if(identical(check_for_default, snp_scale)){
    print("Defaulting to plotting microdiversity analyses according to local SNP calls. Set -snp_scale local, global, or both to change this.")
    snp_scale <- "local"
  }else{
    snp_scale <- options[snp_scale + 1]
  }
  
  
  
  run_parameters <- read.csv("metapop_run_settings/run_settings.tsv", sep = "\t")
  run_parameters$parameter <- as.character(run_parameters$parameter)
  run_parameters$setting <- as.character(run_parameters$setting)
  
  library_location <- run_parameters$setting[run_parameters$parameter == "Library Location"]
  threads <- as.numeric(run_parameters$setting[run_parameters$parameter == "Threads"])
  original_assembly <- run_parameters$setting[run_parameters$parameter == "Assembly"]
  
  
  # SNP positions in codons code
  
  codon_pos_count <- fread("metapop_microdiversity/codon_position_summary.tsv", sep = "\t")
  
  
  codon_pos_count <- codon_pos_count[, list(sum(first_pos), sum(second_pos), sum(third_pos)), by = source]
  
  codon_pos_count_all <- melt.data.table(codon_pos_count, id.vars = c("source"))
  codon_pos_count_all <- codon_pos_count_all[, sum(value), by = source]
  
  p_all <- ggplot(codon_pos_count_all, aes(x = source, y = V1)) +
    geom_bar(stat = "identity", position = "stack") +
    ylab("Count of SNPs in sample") +
    xlab("Sample of origin") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90), 
          axis.title.x = element_blank(),
          axis.title = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          axis.line = element_line("black"),
          axis.ticks.y = element_line("black"),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14))+
    scale_y_continuous(labels = scales::comma)
  
  
  codon_pos_count[, sum := V1+V2+V3,]
  
  codon_pos_count[, V1 := V1/sum]
  codon_pos_count[, V2 := V2/sum]
  codon_pos_count[, V3 := V3/sum]
  
  codon_pos_count[, sum := NULL]
  
  colnames(codon_pos_count)[2:4] = c("First Pos.", "Second Pos.", "Third Pos.")
  
  codon_pos_count <- melt.data.table(codon_pos_count, id.vars = c("source"))
  
  p <- ggplot(codon_pos_count, aes(x = source, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "stack") +
    ylab("Proportion of SNPs by position in codon") +
    xlab("Sample of origin") +
    theme_minimal() +
    scale_fill_manual("Position\nof SNP in\nCodon", labels = c("Pos. 1", "Pos. 2", "Pos. 3"), values = c("grey50", "#2ca9e1", "#FF0000")) + 
    theme(axis.text.x = element_text(angle = 90), 
          axis.title.x = element_blank(),
          axis.title = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          axis.line = element_line("black"),
          axis.ticks.y = element_line("black"),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14))
  
  
  if(!dir.exists("metapop_visualizations")){
    system("mkdir metapop_visualizations")
  }
  
  pdf("metapop_visualizations/Third_Pos_SNP_summary.pdf", 11, 11)
  
  print(p)
  
  print(p_all)
  
  dev.off()
  
  #FST
  
  fixation_data <- fread("metapop_microdiversity/fixation_index.tsv", sep = "\t")
  fixation_data_filtered <- fixation_data[!is.na(fixation_data$fst),]
  contigs_names <- unique(fixation_data_filtered$contig)
  
  
  namesave = unique(fixation_data_filtered$contig)
  fixation_data <- fixation_data_filtered[, list(list(.SD)), by = contig]$V1
  
  cl <- makeCluster(min(threads, detectCores()))
  clusterExport(cl, varlist = c("library_location", "fixation_data"), envir = environment())
  clusterEvalQ(cl, suppressMessages(suppressWarnings(library(data.table, lib.loc = library_location))))
  
  registerDoParallel(cl)
  
  fixation_data <- foreach(x = 1:length(fixation_data)) %dopar% {
    x <- fixation_data[[x]]
    
    tmp_contigs <- unique(c(x$row_samp, x$col_samp))
    
    x2 <- x[, list(col_samp, row_samp, fst) ]
    colnames(x2)[1:2] = c("row_samp", "col_samp")
    
    missing <- data.table(row_samp = tmp_contigs, col_samp = tmp_contigs, fst = 0)
    
    x <- rbind(x, x2, missing)
    
    organize_rows_and_cols <- sort(unique(c(x$row_samp, x$col_samp)))
    
    x$row_samp <- factor(x$row_samp, levels = organize_rows_and_cols)
    x$col_samp <- factor(x$col_samp, levels = organize_rows_and_cols)
    
    #Silly, but works as a method of ensuring that the plotting output is correct in shape.
    x2 <- dcast(x, col_samp ~ row_samp, value.var ="fst")
    rownames(x2) = x2$col_samp
    x2[, col_samp := NULL]
    
    x2[upper.tri(x2)] <- NA
    
    x2[, col_samp := rownames(x2)]
    
    x <- melt.data.table(x2, na.rm = T, id.vars = c("col_samp"))
    
    colnames(x)[2:3] = c("row_samp", "fst")
    
    
    return(x)
  }
  
  stopCluster(cl)
  
  
  names(fixation_data) = namesave
  
  groups <- (1:length(fixation_data))%/%threads
  unique_groups <- unique(groups)
  
  cl <- makeCluster(min(threads, detectCores()))
  clusterExport(cl, varlist = c("library_location", "groups", "unique_groups"), envir = environment())
  clusterEvalQ(cl, suppressMessages(suppressWarnings(library(data.table, lib.loc = library_location))))
  clusterEvalQ(cl, suppressMessages(suppressWarnings(library(ggplot2, lib.loc = library_location))))
  clusterEvalQ(cl, suppressMessages(suppressWarnings(library(RColorBrewer, lib.loc = library_location))))
  registerDoParallel(cl)
  
  pdf("metapop_visualizations/global_fst_genome_heatmap_plots.pdf", width = 17, height = 11)
  
  for(k in unique_groups){
    
    heatmaps <- foreach(i = fixation_data[groups==k], j = names(fixation_data)[groups==k]) %dopar% {
      
      fst_heatmap <- ggplot(i, aes(y = row_samp, x = col_samp, fill = fst))+
        geom_raster()+
        theme_classic()+
        theme(plot.title = element_text(size = 18, face = "bold"),
              axis.text.x = element_text(angle = 90, hjust = 0.5),
              axis.text.y = element_text(vjust = 0.5),
              plot.margin = unit(c(0,0,0,0), "cm"),
              legend.title = element_text(size = 16),
              legend.text = element_text(size = 14))+ 
        ylab("")+
        xlab("")+
        ggtitle(j)+
        scale_fill_gradient2(low = "#2ca9e1", mid="grey80", high="#ff0000", na.value = "black", midpoint = 0.5, limits= c(0,1))+
        guides(fill = guide_legend(title = "Fst"))
      
      return(fst_heatmap)
      
    }
    
    for(i in heatmaps){
      print(i)
    }
    
  }
  
  stopCluster(cl)
  
  dev.off()
  
  
  #main genomes
  
  geneFile <- as.character(run_parameters$setting[run_parameters$parameter=="Genes"])
  
  
  if(snp_scale=="local" | snp_scale == "both"){
    
    print("Creating local scale SNP plots...")
    
    genes <- readDNAStringSet(geneFile)
    
    s <- strsplit(names(genes), "[# \t]+") # split names by tab/space
    genes <- data.table(matrix(unlist(s), ncol=5, byrow=T))
    
    names(genes)[1:4] = c("contig_gene", "start", "end", "OC")
    
    genes$start <- as.numeric((genes$start))
    genes$end <- as.numeric((genes$end))
    
    genes <- genes[,-5]
    
    genes$parent_contig <- gsub("_\\d+$", "", genes$contig_gene)
    
    gene_microdiv <- fread("metapop_microdiversity/local_gene_microdiversity.tsv", sep = "\t")
    
    gene_microdiv$parent_contig <- gsub("_\\d+$", "", gene_microdiv$contig_gene)
    
    setkeyv(gene_microdiv, c("parent_contig", "source"))
    
    if(!plot_all){
      highest_selected_contigs <- gene_microdiv[, valuable <- sum(pNpS_ratio > 1, na.rm = T), by = key(gene_microdiv)]
      highest_selected_contigs <- highest_selected_contigs[order(highest_selected_contigs$source, V1),]
      retained_contigs <- unique(highest_selected_contigs[, tail(parent_contig, 3), by = source]$V1)
      gene_microdiv <- gene_microdiv[gene_microdiv$parent_contig %in% retained_contigs,]
    }
    
    num_src <- length(unique(gene_microdiv$source))
    
    genes <- genes[genes$parent_contig %in% gene_microdiv$parent_contig,]
    
    genes$OC <- as.numeric(genes$OC)
    reverser <- genes[,list(start = ifelse(OC < 0, end, start), end = ifelse(OC < 0, start, end))]
    genes$start <- reverser$start
    genes$end <- reverser$end
    rm(reverser)
    genes$pnps <- NA
    genes$pi <- NA
    genes$theta <- NA
    genes$tajD <- NA
    genes[, midpt := (start+end)/2]
    
    unique_contigs <- genes[, list(list(.SD)), by = parent_contig]$V1
    names(unique_contigs) = unique(genes$parent_contig)
    rm(genes)
    
    
    NS = unique(gene_microdiv$parent_contig)
    
    gene_microdiv <- gene_microdiv[, list(list(.SD)), by = parent_contig]$V1
    names(gene_microdiv) = NS
    rm(NS)
    
    
    fastaLengths <- readDNAStringSet(original_assembly)
    fastaLengths <- data.table(contig = names(fastaLengths), num_bases = lengths(fastaLengths))
    fastaLengths$contig <- unlist(lapply(fastaLengths$contig, function(x){return(strsplit(x, split=" ")[[1]][1])}))
    
    fastaLengths <- fastaLengths[ contig %in% names(gene_microdiv),]
    
    fastaLengths[, bigness := round(log10(num_bases))]
    
    setkey(fastaLengths, "bigness")
    
    
    unique_contigs <- unique_contigs[fastaLengths$contig]
    
    
    depth_info <- list.files(path = "metapop_depth_by_pos", full.names = T)
    depth_names <- substr(depth_info, 22, nchar(depth_info)-17)
    
    
    cl <- makeCluster(min(threads, detectCores()))
    
    clusterExport(cl, varlist = c("library_location", "depth_info"), envir = environment())
    clusterEvalQ(cl, suppressMessages(suppressWarnings(library(data.table, lib.loc = library_location))))
    
    registerDoParallel(cl)
    
    depth_bins <- foreach(x = depth_info) %dopar% {
      
      tmp <- fread(x, sep = "\t")
      colnames(tmp) = c("contig", "pos", "depth")
      tmp[,pos := ((pos %/% 250))*250,]
      
      setkeyv(tmp, c("contig", "pos"))
      contig_names <- unique(tmp$contig)
      
      tmp <- tmp[, sum(depth)/250, by = key(tmp)]
      
    }
    
    stopCluster(cl)
    
    names(depth_bins) = depth_names
    
    
    size_category <- fastaLengths[, list(list(contig)), by = bigness]
    
    
    
    #Tajima's D plot needs a legend to specify the annotations, but without coloring the points. The easiest way to do this is to make a fake one and save it for later.
    taj_d_legend_plot <- ggplot(data = data.table(dots = c(1,2,3)), aes(x = dots, fill = factor(dots))) +
      geom_bar()+
      scale_fill_manual(name = "Tajima's D\nSelection", values = alpha(c("red","grey50", "lightblue"), 0.35), labels = c("Positive", "Neutral", "Purifying"))+
      theme(legend.text = element_text(size = 14),
            legend.title = element_text(size = 14))
    
    
    taj_d_legend <- get_legend(taj_d_legend_plot)
    
    
    #dev
    #k <- 0
    #i<-1
    
    #print(size_category)
    
    cl <- makeCluster(min(detectCores(), threads, length(unique_contigs)))
    
    clusterExport(cl, varlist = c("library_location", "groups", "unique_groups", "depth_info", "depth_names", "depth_bins"), envir=environment())
    
    clusterEvalQ(cl, library(data.table, lib.loc = library_location))
    clusterEvalQ(cl, library(ggplot2, lib.loc = library_location))
    clusterEvalQ(cl, library(gggenes, lib.loc = library_location))
    clusterEvalQ(cl, library(cowplot, lib.loc = library_location))
    
    
    registerDoParallel(cl)
    
    for(z in size_category$bigness){
      
      current_conts <- unlist(size_category$V1[size_category$bigness == z])
      
      groups <- (1:length(current_conts)) %/% threads
      unique_groups <- unique(groups)
      
      pdf(paste0("metapop_visualizations/local_contigs_bp_cat_",as.character(format(10^z, scientific = F)),"_microdiversity_viz.pdf"), width = 13 + (2^(z+1)), height = 11)
      
      for(k in unique_groups){
        
        all_contigs <- foreach(i = 1:length(unique_contigs[current_conts][groups == k])) %dopar% {
          
          name <- names(unique_contigs[current_conts][groups == k])[i]
          sub <- unique_contigs[[name]]
          sub_data <- gene_microdiv[[name]]
          
          full_len <- fastaLengths[contig == name, num_bases]
          
          #Extremely high pNpS is largely irrelevant to answering the question of selection. This fix keeps scaling consistent, but always displays high pNpS as pos. select.
          sub_data$pNpS_ratio <- ifelse(sub_data$pNpS_ratio > 2, 2, sub_data$pNpS_ratio)
          sub_data$pNpS_ratio <- ifelse(!sub_data$snps_present, 0, sub_data$pNpS_ratio)
          
          a_contig <- lapply(unique(sub_data$source), function(j){
            
            depth_dat <- depth_bins[[j]][depth_bins[[j]]$contig == name,]
            
            missing_bins <- seq(0, full_len, by = 250)
            missing_bins <- missing_bins[!missing_bins %in% depth_dat$pos]
            
            if(length(missing_bins) > 0){
              missing_bins <- data.table(contig = name, pos = missing_bins, V1 = 0)
              depth_dat <- rbind(depth_dat, missing_bins) 
            }
            
            sub$pnps[match(sub_data$contig_gene[sub_data$source == j], sub$contig_gene)] <- sub_data$pNpS_ratio[sub_data$source == j]
            sub$pi[match(sub_data$contig_gene[sub_data$source == j], sub$contig_gene)] <- sub_data$pi[sub_data$source == j]
            sub$theta[match(sub_data$contig_gene[sub_data$source == j], sub$contig_gene)] <- sub_data$theta[sub_data$source == j]
            sub$tajD[match(sub_data$contig_gene[sub_data$source == j], sub$contig_gene)] <- sub_data$taj_D[sub_data$source == j]
            
            pnps_fill <- any(!is.na(sub$pnps))
            pi_fill <- any(!is.na(sub$pi))
            theta_fill <- any(!is.na(sub$theta))
            tajD_fill <- any(!is.na(sub$tajD))
            
            
            if(tajD_fill){
              sub[is.infinite(tajD), tajD := NA]
            }
            
            labs <- c("Purifying (0)", "Neutral (1)", "Positive (>1)")
            
            brk <- c(0, 1, 2)
            
            genes_plot <- ggplot(sub, aes(fill=pnps, xmin = start, xmax = end, y = OC/2)) +
              geom_gene_arrow() +
              xlab("")+
              ylab("Strand") +
              scale_x_continuous(expand = c(0,0), labels = scales::comma)+
              scale_y_continuous(breaks = c(-0.5, 0.5), labels = c("-", "+"), limits = c(-.75, 0.75)) +
              scale_fill_gradient2(low="lightblue", mid = "white", high="red", na.value="black", breaks = brk, labels = labs, name = "pN/pS Ratio", limits = c(0, 2), midpoint = 1)+
              theme_genes() +
              ggtitle("pN/pS Selection by Gene")+
              theme(axis.text.y = element_text(size = 20, vjust = 0.38),
                    axis.title.y = element_text(size = 16),
                    axis.text.x = element_text(size = 14),
                    legend.text = element_text(size = 14),
                    legend.title = element_text(size = 14),
                    title = element_text(size = 14)) + 
              scale_color_manual(values = 'black', labels = 'No pNpS') +
              guides(color = guide_legend(override.aes = list(fill = "black")))
            #Adds in pi, theta, and tajima's D points. Offset by their section of the plot.
            
            gene_legend <- get_legend(genes_plot)
            
            genes_plot <- genes_plot + theme(legend.position = "none",
                                             axis.line.x.bottom = element_line("black"))
            
            depth_by_pos_plot <- ggplot(depth_dat, aes(x = pos, y = V1))+
              geom_step() +
              theme_minimal()+
              ylab("Depth of Coverage")+
              xlab("") +
              scale_x_continuous(expand = c(0,0), labels = scales::comma)+
              ggtitle("Avg. Depth of Coverage over Contig")+
              theme(axis.title.y = element_text(size = 16),
                    axis.text.x = element_text(size = 14),
                    axis.line = element_line("black"),
                    title = element_text(size = 14),
                    axis.ticks = element_line("black"))
            
            thet_dat <- melt.data.table(sub, id.vars = "midpt", measure.vars = c("pi", "theta"))
            
            pi_and_theta_plot <- ggplot(thet_dat, aes(x = midpt, y = value, color = variable))+
              geom_point(size = 2) +
              theme_minimal() +
              xlab("")+
              ylab("Pi and Theta")+
              scale_x_continuous(expand = c(0,0), labels = scales::comma)+
              scale_color_manual("Nucleotide\nDiversity\nMeasure", labels = c("Pi", "Theta"), values = c("#af8dc3", "#7fbf7b")) +
              ggtitle("Pi and Theta Nucleotide Diversity by Gene")+
              theme(axis.title.y = element_text(size = 16),
                    axis.text.x = element_text(size = 14),
                    legend.text = element_text(size = 14),
                    legend.title = element_text(size = 14),
                    axis.line = element_line("black"),
                    title = element_text(size = 14),
                    axis.ticks = element_line("black"))
            
            thet_legend <- get_legend(pi_and_theta_plot)
            
            pi_and_theta_plot <- pi_and_theta_plot + theme(legend.position = "none")
            
            
            #Needs a fill legend for the colors.
            taj_d_plot <- ggplot(sub, aes(x = midpt))+
              geom_point(aes(y = tajD), size = 2) +
              theme_minimal() +
              xlab("")+
              scale_x_continuous(expand = c(0,0), labels = scales::comma)+
              ylab("Tajima's D") +
              ylim(c(min(sub$tajD), max(sub$tajD)))+
              annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = min(-2, max(sub$tajD)), fill = "lightblue", alpha = 0.35)+
              annotate("rect", xmin = -Inf, xmax = Inf, ymin = max(-2, min(sub$tajD)), ymax = min(2, max(sub$tajD)), fill = "grey50", alpha = 0.35)+
              annotate("rect", xmin = -Inf, xmax = Inf, ymin = max(2, min(sub$tajD)), ymax = Inf, fill = "red", alpha = 0.35) +
              ggtitle("Tajima's D Selection by Gene")+
              theme(axis.title.y = element_text(size = 16),
                    axis.text.x = element_text(size = 14),
                    legend.text = element_text(size = 14),
                    legend.title = element_text(size = 14),
                    axis.line = element_line("black"),
                    title = element_text(size = 14),
                    axis.ticks = element_line("black"))
            
            taj_d_plot <- taj_d_plot + theme(legend.position = "none")
            
            
            taj_d_plot <- taj_d_plot+ coord_cartesian(xlim = c(0,max(sub$end)))
            pi_and_theta_plot <- pi_and_theta_plot + coord_cartesian(xlim = c(0,max(sub$end)))
            depth_by_pos_plot <- depth_by_pos_plot+ coord_cartesian(xlim = c(0,max(sub$end)))
            genes_plot <- genes_plot + coord_cartesian(xlim = c(0,max(sub$end)))
            
            genes_plot <- plot_grid( depth_by_pos_plot,genes_plot, pi_and_theta_plot, taj_d_plot, align = "vh", ncol = 1, axis = "x")
            
            
            legends <- plot_grid(NULL, gene_legend, thet_legend, taj_d_legend, ncol = 1,nrow = 4, align = "vh")
            
            legend_area <- 0.1-(0.01 * z)
            
            title <- ggdraw() + 
              draw_label(
                paste("Contig:", name, "Sample:", j),
                fontface = 'bold',
                x = 0,
                hjust = 0
              ) +
              theme(
                # add margin on the left of the drawing canvas,
                # so title is aligned with left edge of first plot
                plot.margin = margin(0, 0, 0, 7)
              )
            
            genes_plot <- plot_grid(genes_plot, legends, ncol = 2, rel_widths = c(1-legend_area, legend_area))
            
            genes_plot <- ggdraw(add_sub(genes_plot, "Contig Position (bp)", vpadding=grid::unit(0,"lines"),y=6, x=0.5, vjust=4.5))
            
            genes_plot <- plot_grid(
              title, genes_plot,
              ncol = 1,
              # rel_heights values control vertical title margins
              rel_heights = c(0.07, 1)
            )
            
            
            return(genes_plot)
            
          })
          
          return(a_contig)
          
        }
        
        silent <- suppressWarnings(lapply(all_contigs, function(f){
          
          suppressWarnings(lapply(f, print))
          
          return(NA)
          
        }))
        
      }
      
      dev.off()  
      
    }
    
  }
  
  if(snp_scale=="global" | snp_scale == "both"){
    
    print("Creating global scale SNP plots...")
    
    
    genes <- readDNAStringSet(geneFile)
    
    s <- strsplit(names(genes), "[# \t]+") # split names by tab/space
    genes <- data.table(matrix(unlist(s), ncol=5, byrow=T))
    
    names(genes)[1:4] = c("contig_gene", "start", "end", "OC")
    
    genes$start <- as.numeric((genes$start))
    genes$end <- as.numeric((genes$end))
    
    genes <- genes[,-5]
    
    genes$parent_contig <- gsub("_\\d+$", "", genes$contig_gene)
    
    gene_microdiv <- fread("metapop_microdiversity/global_gene_microdiversity.tsv", sep = "\t")
    
    gene_microdiv$parent_contig <- gsub("_\\d+$", "", gene_microdiv$contig_gene)
    
    setkeyv(gene_microdiv, c("parent_contig", "source"))
    
    if(!plot_all){
      highest_selected_contigs <- gene_microdiv[, valuable <- sum(pNpS_ratio > 1, na.rm = T), by = key(gene_microdiv)]
      highest_selected_contigs <- highest_selected_contigs[order(highest_selected_contigs$source, V1),]
      retained_contigs <- unique(highest_selected_contigs[, tail(parent_contig, 3), by = source]$V1)
      gene_microdiv <- gene_microdiv[gene_microdiv$parent_contig %in% retained_contigs,]
    }
    
    num_src <- length(unique(gene_microdiv$source))
    
    genes <- genes[genes$parent_contig %in% gene_microdiv$parent_contig,]
    
    genes$OC <- as.numeric(genes$OC)
    reverser <- genes[,list(start = ifelse(OC < 0, end, start), end = ifelse(OC < 0, start, end))]
    genes$start <- reverser$start
    genes$end <- reverser$end
    rm(reverser)
    genes$pnps <- NA
    genes$pi <- NA
    genes$theta <- NA
    genes$tajD <- NA
    genes[, midpt := (start+end)/2]
    
    unique_contigs <- genes[, list(list(.SD)), by = parent_contig]$V1
    names(unique_contigs) = unique(genes$parent_contig)
    rm(genes)
    
    
    NS = unique(gene_microdiv$parent_contig)
    
    gene_microdiv <- gene_microdiv[, list(list(.SD)), by = parent_contig]$V1
    names(gene_microdiv) = NS
    rm(NS)
    
    
    fastaLengths <- readDNAStringSet(original_assembly)
    fastaLengths <- data.table(contig = names(fastaLengths), num_bases = lengths(fastaLengths))
    fastaLengths$contig <- unlist(lapply(fastaLengths$contig, function(x){return(strsplit(x, split=" ")[[1]][1])}))
    
    fastaLengths <- fastaLengths[ contig %in% names(gene_microdiv),]
    
    fastaLengths[, bigness := round(log10(num_bases))]
    
    setkey(fastaLengths, "bigness")
    
    
    unique_contigs <- unique_contigs[fastaLengths$contig]
    
    
    depth_info <- list.files(path = "metapop_depth_by_pos", full.names = T)
    depth_names <- substr(depth_info, 22, nchar(depth_info)-17)
    
    
    cl <- makeCluster(min(threads, detectCores()))
    
    clusterExport(cl, varlist = c("library_location", "depth_info"), envir = environment())
    clusterEvalQ(cl, suppressMessages(suppressWarnings(library(data.table, lib.loc = library_location))))
    
    registerDoParallel(cl)
    
    depth_bins <- foreach(x = depth_info) %dopar% {
      
      tmp <- fread(x, sep = "\t")
      colnames(tmp) = c("contig", "pos", "depth")
      tmp[,pos := ((pos %/% 250))*250,]
      
      setkeyv(tmp, c("contig", "pos"))
      contig_names <- unique(tmp$contig)
      
      tmp <- tmp[, sum(depth)/250, by = key(tmp)]
      
    }
    
    stopCluster(cl)
    
    names(depth_bins) = depth_names
    
    
    size_category <- fastaLengths[, list(list(contig)), by = bigness]
    
    
    
    #Tajima's D plot needs a legend to specify the annotations, but without coloring the points. The easiest way to do this is to make a fake one and save it for later.
    taj_d_legend_plot <- ggplot(data = data.table(dots = c(1,2,3)), aes(x = dots, fill = factor(dots))) +
      geom_bar()+
      scale_fill_manual(name = "Tajima's D\nSelection", values = alpha(c("red","grey50", "lightblue"), 0.35), labels = c("Positive", "Neutral", "Purifying"))+
      theme(legend.text = element_text(size = 14),
            legend.title = element_text(size = 14))
    
    
    taj_d_legend <- get_legend(taj_d_legend_plot)
    
    
    
    cl <- makeCluster(min(detectCores(), threads, length(unique_contigs)))
    
    clusterExport(cl, varlist = c("library_location", "groups", "unique_groups", "depth_info", "depth_names", "depth_bins"), envir=environment())
    
    clusterEvalQ(cl, library(data.table, lib.loc = library_location))
    clusterEvalQ(cl, library(ggplot2, lib.loc = library_location))
    clusterEvalQ(cl, library(gggenes, lib.loc = library_location))
    clusterEvalQ(cl, library(cowplot, lib.loc = library_location))
    
    
    registerDoParallel(cl)
    
    for(z in size_category$bigness){
      
      current_conts <- unlist(size_category$V1[size_category$bigness == z])
      
      groups <- (1:length(current_conts)) %/% threads
      unique_groups <- unique(groups)
      
      pdf(paste0("metapop_visualizations/global_contigs_bp_cat_",as.character(format(10^z, scientific = F)),"_microdiversity_viz.pdf"), width = 13 + (2^(z+1)), height = 11)
      
      for(k in unique_groups){
        
        all_contigs <- foreach(i = 1:length(unique_contigs[current_conts][groups == k])) %dopar% {
          
          name <- names(unique_contigs[current_conts][groups == k])[i]
          sub <- unique_contigs[[name]]
          sub_data <- gene_microdiv[[name]]
          
          full_len <- fastaLengths[contig == name, num_bases]
          
          #Extremely high pNpS is largely irrelevant to answering the question of selection. This fix keeps scaling consistent, but always displays high pNpS as pos. select.
          sub_data$pNpS_ratio <- ifelse(sub_data$pNpS_ratio > 2, 2, sub_data$pNpS_ratio)
          sub_data$pNpS_ratio <- ifelse(!sub_data$snps_present, 0, sub_data$pNpS_ratio)
          
          a_contig <- lapply(unique(sub_data$source), function(j){
            
            depth_dat <- depth_bins[[j]][depth_bins[[j]]$contig == name,]
            
            missing_bins <- seq(0, full_len, by = 250)
            missing_bins <- missing_bins[!missing_bins %in% depth_dat$pos]
            
            if(length(missing_bins) > 0){
              missing_bins <- data.table(contig = name, pos = missing_bins, V1 = 0)
              depth_dat <- rbind(depth_dat, missing_bins) 
            }
            
            sub$pnps[match(sub_data$contig_gene[sub_data$source == j], sub$contig_gene)] <- sub_data$pNpS_ratio[sub_data$source == j]
            sub$pi[match(sub_data$contig_gene[sub_data$source == j], sub$contig_gene)] <- sub_data$pi[sub_data$source == j]
            sub$theta[match(sub_data$contig_gene[sub_data$source == j], sub$contig_gene)] <- sub_data$theta[sub_data$source == j]
            sub$tajD[match(sub_data$contig_gene[sub_data$source == j], sub$contig_gene)] <- sub_data$taj_D[sub_data$source == j]
            
            pnps_fill <- any(!is.na(sub$pnps))
            pi_fill <- any(!is.na(sub$pi))
            theta_fill <- any(!is.na(sub$theta))
            tajD_fill <- any(!is.na(sub$tajD))
            
            
            if(tajD_fill){
              sub[is.infinite(tajD), tajD := NA]
            }
            
            labs <- c("Purifying (0)", "Neutral (1)", "Positive (>1)")
            
            brk <- c(0, 1, 2)
            
            genes_plot <- ggplot(sub, aes(fill=pnps, xmin = start, xmax = end, y = OC/2)) +
              geom_gene_arrow() +
              xlab("")+
              ylab("Strand") +
              scale_x_continuous(expand = c(0,0), labels = scales::comma)+
              scale_y_continuous(breaks = c(-0.5, 0.5), labels = c("-", "+"), limits = c(-.75, 0.75)) +
              scale_fill_gradient2(low="lightblue", mid = "white", high="red", na.value="black", breaks = brk, labels = labs, name = "pN/pS Ratio", limits = c(0, 2), midpoint = 1)+
              theme_genes() +
              ggtitle("pN/pS Selection by Gene")+
              theme(axis.text.y = element_text(size = 20, vjust = 0.38),
                    axis.title.y = element_text(size = 16),
                    axis.text.x = element_text(size = 14),
                    legend.text = element_text(size = 14),
                    legend.title = element_text(size = 14),
                    title = element_text(size = 14)) + 
              scale_color_manual(values = 'black', labels = 'No pNpS') +
              guides(color = guide_legend(override.aes = list(fill = "black")))
            #Adds in pi, theta, and tajima's D points. Offset by their section of the plot.
            
            gene_legend <- get_legend(genes_plot)
            
            genes_plot <- genes_plot + theme(legend.position = "none",
                                             axis.line.x.bottom = element_line("black"))
            
            depth_by_pos_plot <- ggplot(depth_dat, aes(x = pos, y = V1))+
              geom_step() +
              theme_minimal()+
              ylab("Depth of Coverage")+
              xlab("") +
              scale_x_continuous(expand = c(0,0), labels = scales::comma)+
              ggtitle("Avg. Depth of Coverage over Contig")+
              theme(axis.title.y = element_text(size = 16),
                    axis.text.x = element_text(size = 14),
                    axis.line = element_line("black"),
                    title = element_text(size = 14),
                    axis.ticks = element_line("black"))
            
            thet_dat <- melt.data.table(sub, id.vars = "midpt", measure.vars = c("pi", "theta"))
            
            pi_and_theta_plot <- ggplot(thet_dat, aes(x = midpt, y = value, color = variable))+
              geom_point(size = 2) +
              theme_minimal() +
              xlab("")+
              ylab("Pi and Theta")+
              scale_x_continuous(expand = c(0,0), labels = scales::comma)+
              scale_color_manual("Nucleotide\nDiversity\nMeasure", labels = c("Pi", "Theta"), values = c("#af8dc3", "#7fbf7b")) +
              ggtitle("Pi and Theta Nucleotide Diversity by Gene")+
              theme(axis.title.y = element_text(size = 16),
                    axis.text.x = element_text(size = 14),
                    legend.text = element_text(size = 14),
                    legend.title = element_text(size = 14),
                    axis.line = element_line("black"),
                    title = element_text(size = 14),
                    axis.ticks = element_line("black"))
            
            thet_legend <- get_legend(pi_and_theta_plot)
            
            pi_and_theta_plot <- pi_and_theta_plot + theme(legend.position = "none")
            
            
            #Needs a fill legend for the colors.
            taj_d_plot <- ggplot(sub, aes(x = midpt))+
              geom_point(aes(y = tajD), size = 2) +
              theme_minimal() +
              xlab("")+
              scale_x_continuous(expand = c(0,0), labels = scales::comma)+
              ylab("Tajima's D") +
              ylim(c(min(sub$tajD), max(sub$tajD)))+
              annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = min(-2, max(sub$tajD)), fill = "lightblue", alpha = 0.35)+
              annotate("rect", xmin = -Inf, xmax = Inf, ymin = max(-2, min(sub$tajD)), ymax = min(2, max(sub$tajD)), fill = "grey50", alpha = 0.35)+
              annotate("rect", xmin = -Inf, xmax = Inf, ymin = max(2, min(sub$tajD)), ymax = Inf, fill = "red", alpha = 0.35) +
              ggtitle("Tajima's D Selection by Gene")+
              theme(axis.title.y = element_text(size = 16),
                    axis.text.x = element_text(size = 14),
                    legend.text = element_text(size = 14),
                    legend.title = element_text(size = 14),
                    axis.line = element_line("black"),
                    title = element_text(size = 14),
                    axis.ticks = element_line("black"))
            
            taj_d_plot <- taj_d_plot + theme(legend.position = "none")
            
            
            taj_d_plot <- taj_d_plot+ coord_cartesian(xlim = c(0,max(sub$end)))
            pi_and_theta_plot <- pi_and_theta_plot + coord_cartesian(xlim = c(0,max(sub$end)))
            depth_by_pos_plot <- depth_by_pos_plot+ coord_cartesian(xlim = c(0,max(sub$end)))
            genes_plot <- genes_plot + coord_cartesian(xlim = c(0,max(sub$end)))
            
            genes_plot <- plot_grid( depth_by_pos_plot,genes_plot, pi_and_theta_plot, taj_d_plot, align = "vh", ncol = 1, axis = "x")
            
            
            legends <- plot_grid(NULL, gene_legend, thet_legend, taj_d_legend, ncol = 1,nrow = 4, align = "vh")
            
            legend_area <- 0.1-(0.01 * z)
            
            title <- ggdraw() + 
              draw_label(
                paste("Contig:", name, "Sample:", j),
                fontface = 'bold',
                x = 0,
                hjust = 0
              ) +
              theme(
                # add margin on the left of the drawing canvas,
                # so title is aligned with left edge of first plot
                plot.margin = margin(0, 0, 0, 7)
              )
            
            genes_plot <- plot_grid(genes_plot, legends, ncol = 2, rel_widths = c(1-legend_area, legend_area))
            
            genes_plot <- ggdraw(add_sub(genes_plot, "Contig Position (bp)", vpadding=grid::unit(0,"lines"),y=6, x=0.5, vjust=4.5))
            
            genes_plot <- plot_grid(
              title, genes_plot,
              ncol = 1,
              # rel_heights values control vertical title margins
              rel_heights = c(0.07, 1)
            )
            
            
            return(genes_plot)
            
          })
          
          return(a_contig)
          
        }
        
        silent <- suppressWarnings(lapply(all_contigs, function(f){
          
          suppressWarnings(lapply(f, print))
          
          return(NA)
          
        }))
        
      }
      
      dev.off()  
      
    }
    
  }
  
}
#CB Viz
metapop_cbviz <- function(options){
  
  run_parameters <- read.csv("metapop_run_settings/run_settings.tsv", sep = "\t")
  run_parameters$parameter <- as.character(run_parameters$parameter)
  run_parameters$setting <- as.character(run_parameters$setting)
  
  library_location <- run_parameters$setting[run_parameters$parameter == "Library Location"]
  threads <- as.numeric(run_parameters$setting[run_parameters$parameter == "Threads"])
  genes_file <- run_parameters$setting[run_parameters$parameter == "Genes"]
  
  
  
  passing_contigs <- list.files(path = "metapop_cov_and_depth/", pattern = ".tsv", full.names = T)
  
  min_cov <- as.numeric(run_parameters$setting[run_parameters$parameter=="Coverage"])
  min_dep <- as.numeric(run_parameters$setting[run_parameters$parameter=="Depth"])
  
  cl <- makeCluster(min(detectCores(), threads, length(passing_contigs)))
  
  clusterExport(cl, varlist=c("passing_contigs", "min_cov", "min_dep", "library_location"), envir=environment())
  clusterEvalQ(cl, suppressMessages(suppressWarnings(library(data.table, lib.loc = library_location))))
  
  registerDoParallel(cl)
  
  passing_contigs <- unique(foreach(i=passing_contigs, .combine=c) %dopar% {
    
    tmp <- fread(i, sep="\t")
    
    return(tmp$V1[tmp$V3>=min_cov & tmp$V4 >= min_dep])
    
  })
  
  stopCluster(cl)
  
  #Already per contig
  codon_bias_iqr <- fread(list.files(full.names = T, path = "metapop_codon_bias/", pattern="gene_IQR_and_mean.tsv"), sep = "\t")
  codon_bias_iqr <- codon_bias_iqr[codon_bias_iqr$parent_contig %in% passing_contigs,]
  
  genes <- readDNAStringSet(genes_file)
  
  s <- strsplit(names(genes), "[# \t]+") # split names by tab/space
  genes <- data.table(matrix(unlist(s), ncol=5, byrow=T))
  
  names(genes)[1:4] = c("contig_gene", "start", "end", "OC")
  
  genes$start <- as.numeric((genes$start))
  genes$end <- as.numeric((genes$end))
  genes <- genes[,-5]
  
  # Figure out what contig they come from, mostly for cleaning purposes
  genes$parent_contig <- gsub("_\\d+$", "", genes$contig_gene)  
  
  genes <- genes[genes$parent_contig %in% passing_contigs,]
  
  codon_bias_genes <- fread(list.files(full.names = T, path = "metapop_codon_bias/", pattern="gene_euclidean_distances.tsv"), sep = "\t")
  codon_bias_genes <- codon_bias_genes[codon_bias_genes$parent_contig %in% passing_contigs,]
  
  codon_bias_genes$start <- genes$start[match(codon_bias_genes$gene, genes$contig_gene)]
  codon_bias_genes$end <- genes$end[match(codon_bias_genes$gene, genes$contig_gene)]
  codon_bias_genes$strand <- genes$OC[match(codon_bias_genes$gene, genes$contig_gene)]
  
  rm(genes)
  
  #presplit into contigs
  namesave = unique(codon_bias_genes$parent_contig)
  codon_bias_genes <- codon_bias_genes[, list(list(.SD)), by = parent_contig]$V1
  names(codon_bias_genes) = namesave
  
  codon_bias_genes <- lapply(codon_bias_genes, function(x){
    
    l <- min(x$euc_dist)
    u <- max(x$euc_dist)
    
    x$relative_dist <- 2+(((x$euc_dist - l)/(u-l))*2)
    
    return(x)
  })
  
  if(!all(passing_contigs %in% names(codon_bias_genes))){
    print("There's some contigs in the samples that weren't found in the genes file.")
    print("There may be cases where contigs had no predicted genes, or where the codon bias of a set of genes could not be calculated completely.")
    print("Only contigs with predicted genes can/will be plotted.")
    passing_contigs <- passing_contigs[passing_contigs %in% names(codon_bias_genes)]
  }
  
  
  #If I normalize the gene distance on 0-1, then I can make a consistent width plot
  
  color_legend_plot <- ggplot(data = data.table(dots = c(1,2)), aes(x = dots, fill = factor(dots))) +
    geom_bar()+
    scale_fill_manual(name = "Gene Codon Bias", values = alpha(c("#2ca9e1", "#FF0000"), 1), labels = c("Typical Codon Use", "Abnormal Codon Use"))+
    theme(legend.text = element_text(size = 14),
          legend.title = element_text(size = 14))
  
  color_leg <- get_legend(color_legend_plot)
  
  codon_bias_ringplot <- function(contig, gene_profile, thresholds, leg = color_leg){
    #iqr <- round(thresholds[1], 2)
    #mu <- round(thresholds[2], 2)
    #outlier_bound <- round(thresholds[3], 2)
    
    if(all(thresholds == 0)){
      return(NA)
    }
    
    p <- ggplot(gene_profile, aes(fill=factor(outlier_status), xmin=2, xmax=relative_dist, ymin=start, ymax=end))+
      annotate("rect", xmin = 1.992, xmax = 4, ymin=0, ymax = max(gene_profile$end), fill = "grey65", color = "black") +
      geom_rect() +
      coord_polar(theta="y") +
      xlim(c(0, 4)) +
      ylim(c(0, max(gene_profile$end*4/3))) +
      theme(panel.background = element_blank(), 
            axis.text.y = element_blank(), 
            axis.title.y = element_blank(), 
            axis.ticks.y = element_blank(), 
            axis.text.x = element_blank(), 
            axis.title.x = element_blank(), 
            axis.ticks.x = element_blank(),
            axis.line = element_blank())+
      guides(fill = F) +
      ggtitle(paste0(contig, "\nCodon Usage Bias Distances")) +
      scale_fill_manual(values = c("#2ca9e1", "#FF0000")) +
      annotate("text", y = 0, x = 0, hjust = 0.5, label = paste("Contig Length\n", max(gene_profile$end), " bp", sep = ""), vjust = 0.5)+
      annotate("text", x = 3.1, y = max(gene_profile$end)*1.065, label = "Gene\nEuclidean\nDistance", vjust = 0, hjust = 0.5) +
      
      #tick labels
      annotate("text", y = max(gene_profile$end)*1.014, x = 4, label = paste0(round(max(gene_profile$euc_dist), 3), ""), vjust = 1, hjust = 0, angle = 90)+
      annotate("text", y = max(gene_profile$end)*1.030, x = 2, label = paste0(round(min(gene_profile$euc_dist), 3), ""), vjust = 0, hjust = 0, angle = 90)+
      
      #baseline
      #annotate("segment", x = 2, xend = 4, y = max(gene_profile$end)*1.025, yend = max(gene_profile$end)*1.014) +
      
      #ticks
      annotate("segment", x = 1.992, xend = 1.996, y = max(gene_profile$end), yend = max(gene_profile$end)*1.016) +
      annotate("segment", x = 2.66, xend = 2.66, y = max(gene_profile$end), yend = max(gene_profile$end)*1.012) +
      annotate("segment", x = 3.33, xend = 3.33, y = max(gene_profile$end), yend = max(gene_profile$end)*1.010) +
      annotate("segment", x = 4, xend = 4, y = max(gene_profile$end), yend = max(gene_profile$end)*1.008)
    
    
    p <- plot_grid(NULL, p, leg, NULL, ncol = 4, rel_widths = c(.3, .8, .15, .2))
    
    return(p)
    
  }
  
  groups <- (1:length(passing_contigs))%/%(min(threads, detectCores()))
  unique_groups <- unique(groups)
  
  if(!dir.exists("metapop_visualizations")){
    dir.create("metapop_visualizations")
  }
  
  cl <- makeCluster(min(threads, detectCores()))
  clusterExport(cl, varlist = c("library_location", "groups", "unique_groups", "color_leg"), envir = environment())
  clusterEvalQ(cl, expr = suppressMessages(suppressWarnings(library(data.table, lib.loc = library_location))))
  clusterEvalQ(cl, expr = suppressMessages(suppressWarnings(library(ggplot2, lib.loc = library_location))))
  clusterEvalQ(cl, expr = suppressMessages(suppressWarnings(library(cowplot, lib.loc = library_location))))
  
  pdf("metapop_visualizations/codon_bias_plots.pdf", height = 9, width = 12)
  
  registerDoParallel(cl)
  
  for(k in unique_groups){
    
    CB_plots <- foreach(i = passing_contigs[groups == k]) %dopar% {
      
      codon_bias_ringplot(contig = i, gene_profile = codon_bias_genes[[which(names(codon_bias_genes) == i)]], thresholds =  as.numeric(codon_bias_iqr[which(codon_bias_iqr$parent_contig == i),2:4], leg = color_leg))
      
    }
    
    CB_plots <- CB_plots[!is.na(CB_plots)]
    
    #prevents superfluous outputs
    if(length(CB_plots) > 0){
      for(i in CB_plots){
        print(i)
      }
    }
    
  }
  
  stopCluster(cl)
  
  dev.off()
  
}

#Preproc Viz
metapop_preproc_sum <- function(options){
  
  check_for_default <- integer(0)
  
  run_parameters <- read.csv("metapop_run_settings/run_settings.tsv", sep = "\t")
  run_parameters$parameter <- as.character(run_parameters$parameter)
  run_parameters$setting <- as.character(run_parameters$setting)
  
  library_location <- run_parameters$setting[run_parameters$parameter == "Library Location"]
  threads <- as.numeric(run_parameters$setting[run_parameters$parameter == "Threads"])
  
  
  coverage_cutoff <- as.numeric(run_parameters$setting[run_parameters$parameter == "Coverage"])
  depth_cutoff <- as.numeric(run_parameters$setting[run_parameters$parameter == "Depth"])
  
  samtools_loc <- run_parameters$setting[run_parameters$parameter=="Samtools Location"]
  
  original_files <- list.files(full.names = T, path = "Metapop_BAMs")
  
  cov_depth <- list.files(path="metapop_cov_and_depth", full.names = T, pattern="cov")
  c_d_names <- substring(cov_depth, 23, nchar(cov_depth)-28)
  
  output_bams <- list.files(path = "metapop_preprocessed_reads/", full.names = T)
  output_bams <- output_bams[!grepl(".bai", output_bams)]
  
  cl <- makeCluster(min(threads, detectCores()))
  
  registerDoParallel(cl)
  
  orig_read_counts <- foreach(i = original_files, .combine = c) %dopar% {
    
    as.numeric(system(paste0(samtools_loc, "samtools view -c ", i), intern = T))
    
  }
  
  stopCluster(cl)
  
  overall_summaries <- data.table(sample = unlist(lapply(original_files, function(x){
    strsplit(x, split = "Metapop_BAMs/")[[1]][2]
  })), read_counts = orig_read_counts, finish_counts = 0)
  
  overall_summaries$sample = substr(overall_summaries$sample, 1, nchar(overall_summaries$sample)-4)
  
  cl <- makeCluster(min(threads, detectCores()))
  
  registerDoParallel(cl)
  
  finish_read_counts <- foreach(i = output_bams, .combine = c) %dopar% {
    
    as.numeric(system(paste0(samtools_loc, "samtools view -c ", i), intern = T))
    
  }
  
  stopCluster(cl)
  
  overall_summaries$finish_counts[match(c_d_names, overall_summaries$sample)] <- finish_read_counts
  
  overall_summaries[, removed := read_counts - finish_counts]
  overall_summaries[, read_counts := NULL]
  
  overall_summaries <- melt.data.table(overall_summaries, id.vars = c("sample"))
  overall_summaries <- overall_summaries[order(overall_summaries$variable, decreasing = T),]
  
  read_filtering <- ggplot(overall_summaries, aes(x = sample, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    ylab("Read Count") +
    xlab("Sample Origin") +
    scale_fill_manual("Read\nFiltering", values = c("#2ca9e1", "grey65"), labels = c("Read Kept", "Read Removed")) +
    coord_flip() +
    theme_minimal()
  
  
  
  read_filtering <- read_filtering + theme(axis.title = element_text(size = 14),
                                           axis.text = element_text(size = 14))
  
  title <- ggdraw() + 
    draw_label(
      paste0("Preprocessing Overview\nRead Filtering"),
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )
  
  read_filtering <- plot_grid(title, read_filtering, nrow = 2, rel_heights = c(0.1, 0.9))
  
  cl <- makeCluster(min(threads, detectCores()))
  
  registerDoParallel(cl)
  
  coverage_and_depth <- foreach(i = 1:length(cov_depth), .packages = c("data.table")) %dopar% {
    tmp <- fread(cov_depth[i], sep = "\t")
    tmp[, source := c_d_names[i]]
  }
  
  stopCluster(cl)
  
  
  #Making these legends later is a pain, doing it this way is easy.
  color_legend_plot <- ggplot(data = data.table(dots = c(1,2,3)), aes(x = dots, fill = factor(dots))) +
    geom_bar()+
    scale_fill_manual(name = "Contig\nAssessment", values = alpha(c("grey65","#2ca9e1", "#FF0000"), 1), labels = c("Low Coverage", "Low Depth", "Passing Contig"))+
    theme(legend.text = element_text(size = 14),
          legend.title = element_text(size = 14))
  
  
  color_leg <- get_legend(color_legend_plot)
  
  #We have the data to plot hists of coverage/depth passing, and can show both reads passing and contigs passing.
  
  
  create_summary_viz <- function(depth_data, file_name, file_size, cov=coverage_cutoff, depth=depth_cutoff, leg = color_leg){
    
    group.colors <- c('Insufficient Coverage' = "grey65", 'Insufficient Depth' = "#2CA9E1", 'Passing Contig' = "#FF0000")
    
    depth_data$pass_level <- ifelse(depth_data$V3 >= cov, ifelse(depth_data$V4 >= depth, "Passing Contig", "Insufficient Depth"), "Insufficient Coverage")
    depth_data <- depth_data[order(depth_data$pass_level),]
    
    donut_TF <- depth_data
    donut_TF$coord <- seq(0, 1, length.out = nrow(donut_TF))
    
    #Quirk of geom_rect; if the data is left like so, it will plot a single rectangle for each row of data around the donut. 
    #This will result in visual artifacts, slow plotting speed, and massively oversized PDFs.
    
    donut_TF$ymax = cumsum(donut_TF$coord)
    donut_TF$ymin = c(0, head(donut_TF$ymax, n=-1))
    
    maxes <- donut_TF[, max(ymax, na.rm = T), by = pass_level]
    
    donut_TF <- donut_TF[match(c("Insufficient Coverage", "Insufficient Depth", "Passing Contig")[c("Insufficient Coverage", "Insufficient Depth", "Passing Contig") %in% donut_TF$pass_level], donut_TF$pass_level),]
    
    donut_TF$ymax[match(donut_TF$pass_level, maxes$pass_level)] <- maxes$V1
    
    donut_TF[,ymax:=-ymax]
    donut_TF[,ymin:=-ymin]
    
    title <- ggdraw() + 
      draw_label(
        paste0(file_name, "\nDepth and Coverage of Contigs"),
        fontface = 'bold',
        x = 0,
        hjust = 0
      ) +
      theme(
        # add margin on the left of the drawing canvas,
        # so title is aligned with left edge of first plot
        plot.margin = margin(0, 0, 0, 7)
      )
    
    #Produces a contig-level summary of passage by 
    overview_donut <- ggplot(donut_TF, aes(fill=pass_level, xmin=4.5, xmax=6, ymin=ymin, ymax=ymax))+
      geom_rect()+
      coord_polar(theta="y") +
      xlim(c(0, 6))+
      theme(panel.grid=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(), panel.background = element_blank()) +
      #annotate("text", x = 0, y = 0, label = file_name, size = 3) +
      annotate("text", x = 0.70, y = 0, label = paste("Original Contigs:", nrow(depth_data)),size = 3) +
      annotate("text", x = 0.20, y = 0, label = paste("Contigs with Passing Coverage and TAD:", sum(depth_data$pass_level == "Passing Contig")),size = 3) +
      labs(title="")+
      guides(fill=FALSE) +
      scale_fill_manual(values = group.colors)
    
    overview_donut <- plot_grid(NULL, overview_donut, leg, NULL, ncol = 4, rel_widths = c(.5, .85, .15, .3))
    
    overview_donut <- plot_grid(title, overview_donut, nrow = 2, rel_heights = c(0.1, 0.9))
    
    contig_breakdown <- ggplot(depth_data, aes(x = V3, y = V4, fill = pass_level)) +
      geom_point(alpha = 0.6, shape = 21, color = "grey80", size = 4) +
      geom_vline(xintercept = cov)+
      geom_hline(yintercept = depth)+
      theme_minimal() +
      ylab("Truncated Average Sequencing Depth")+
      xlab("Percent of Genome Covered")+
      scale_fill_manual(values = group.colors)+
      annotate("text", x = cov, y = Inf, label = paste("Minimum Coverage:", cov, " "), vjust=1, hjust = 1) +
      annotate("text", x = Inf, y = depth*(19/20), label = paste("Minimum TAD:", depth), vjust = 1, hjust = 1) +
      expand_limits(x = 0, y = 0)
    
    contig_breakdown <- contig_breakdown + theme(legend.position = "none",
                                                 #axis.line = element_line("grey80"),
                                                 axis.title = element_text(size = 14),
                                                 axis.text = element_text(size = 14),
                                                 panel.grid = element_line(colour = "grey80"))
    
    
    depth_boxplots <- ggplot(depth_data, aes(y = V4, fill = pass_level))+
      geom_boxplot(outlier.color = "grey80", outlier.shape = 21, outlier.size = 4) +
      scale_fill_manual(values = group.colors) +
      ylab("Truncated Average Depth")+
      xlab("")+
      theme_minimal()+
      expand_limits(y = 0)
    
    depth_boxplots <- depth_boxplots + theme(legend.position = "none", 
                                             axis.text = element_blank(), 
                                             #axis.title = element_blank(),
                                             axis.line.y = element_line("grey80"),
                                             axis.ticks = element_blank(),
                                             panel.grid = element_line(colour = "grey80"))
    
    cov_boxplots <- ggplot(depth_data, aes(y = V3, fill = pass_level))+
      geom_boxplot(outlier.color = "grey80", outlier.shape = 21, outlier.size = 4) +
      scale_fill_manual(values = group.colors) +
      xlab("")+
      ylab("Percent of Genome Covered") +
      expand_limits(y = 0) +
      theme_minimal()+
      coord_flip()
    
    cov_boxplots <- cov_boxplots + theme(legend.position = "none", 
                                         axis.text = element_blank(), 
                                         #axis.title = element_blank(),
                                         axis.line.x = element_line("grey80"),
                                         axis.ticks = element_blank(),
                                         panel.grid = element_line(colour = "grey80"))
    
    
    
    c_d_plot <- plot_grid(cov_boxplots, NULL, contig_breakdown, depth_boxplots, ncol = 2, align = "hv", 
                          rel_widths = c(3, 1), rel_heights = c(1, 3), axis = "xy")
    
    legend_plot <- plot_grid(NULL, leg, NULL, NULL, ncol = 2, 
                             rel_widths = c(3, 1), rel_heights = c(1, 4))
    
    c_d_plot <- ggdraw(c_d_plot) + draw_plot(legend_plot)
    
    c_d_plot <- plot_grid(
      title, c_d_plot,
      ncol = 1,
      # rel_heights values control vertical title margins
      rel_heights = c(0.1, 0.9)
    )
    
    return(list(overview_donut, c_d_plot))
    
  }
  
  if(!dir.exists("metapop_visualizations")){
    system("mkdir metapop_visualizations")
  }
  
  pdf("metapop_visualizations/preprocessing_summaries.pdf", width = 16, height = 9)
  
  print(read_filtering)
  
  for(i in 1:length(coverage_and_depth)){
    tmp <- create_summary_viz(coverage_and_depth[[i]], c_d_names[i], file.size(cov_depth[i]))
    print(tmp[[1]])
    print(tmp[[2]])
  }
  
  dev.off()
}
print("Modules loaded.")

#Options:

#select components
dont_micro <- which(grepl("-no_micro", options))
dont_macro <- which(grepl("-no_macro", options))
dont_viz <- which(grepl("-no_viz", options))
preproc_only <- which(grepl("-preprocess_only", options))

macro_only <- which(grepl("-macro_only", options))
micro_only <- which(grepl("-micro_only", options))
viz_only <- which(grepl("-viz_only", options))

library_location <- which(grepl("-lib", options))

if(identical(library_location, check_for_default)){
  print("Metapop will use the default location for R libraries during this session. This will likely fail in HPC environments.")
  flush.console()
  library_location <- .libPaths()
}else{
  library_location <- as.character(options[library_location + 1])
}

if(identical(dont_macro, check_for_default)){
  dont_macro <- F
}else{
  dont_macro <- T
}

if(identical(dont_micro, check_for_default)){
  dont_micro <- F
}else{
  dont_micro <- T
}

if(identical(dont_viz, check_for_default)){
  dont_viz <- F
}else{
  dont_viz <- T
}

if(identical(preproc_only, check_for_default)){
  preproc_only <- F
}else{
  preproc_only <- T
}

if(identical(micro_only, check_for_default)){
  micro_only <- F
}else{
  micro_only <- T
}
if(identical(macro_only, check_for_default)){
  macro_only <- F
}else{
  macro_only <- T
}
if(identical(viz_only, check_for_default)){
  viz_only <- F
}else{
  viz_only <- T
}

print("Loading Libraries...")

#Library_Prep
{
  print(paste("Loading libraries from:",library_location))
  
  suppressMessages(suppressWarnings(library(doParallel, lib.loc = library_location)))
  suppressMessages(suppressWarnings(library(data.table, lib.loc = library_location)))
  suppressMessages(suppressWarnings(library(stringr, lib.loc = library_location)))
  suppressMessages(suppressWarnings(library(Biostrings, lib.loc = library_location)))
  suppressMessages(suppressWarnings(library(Rsamtools, lib.loc = library_location)))
  suppressMessages(suppressWarnings(library(ggplot2, lib.loc = library_location)))
  suppressMessages(suppressWarnings(library(cowplot, lib.loc = library_location)))
  
  suppressMessages(suppressWarnings(library(vegan, lib.loc = library_location)))
  suppressMessages(suppressWarnings(library(compositions, lib.loc = library_location)))
  suppressMessages(suppressWarnings(library(pheatmap, lib.loc = library_location)))
  suppressMessages(suppressWarnings(library(RColorBrewer, lib.loc = library_location)))
  
  suppressMessages(suppressWarnings(library(bit64, lib.loc = library_location)))
  suppressMessages(suppressWarnings(library(gggenes, lib.loc = library_location)))
  
}

print("Libraries loaded.")

if(preproc_only){
  print("Only preprocessing...")
  metapop_preprocess(options)
  quit(save = "no")
}else{
  
  if(micro_only){
    metapop_microdiv(options)
    metapop_microdiv_viz(options)
    quit(save = "no")
  }
  
  if(macro_only){
    metapop_macrodiv(options)
    quit(save = "no")
    
  }
  
  if(viz_only){
    try(metapop_preproc_sum(options))
    try(metapop_microdiv_viz(options))
    try(metapop_cbviz(options))
    try(metapop_macrodiv(options))
    quit(save = "no")
  }
  
  print("Beginning MetaPop... ")
  metapop_preprocess(options)
  
  if(!dont_micro){
    print("Beginning microdiversity")
    
    metapop_call_snps(options)
    
    metapop_mine_reads(options)
    
    metapop_microdiv(options)
  }
  
  if(!dont_macro){
    print("Beginning macrodiversity")
    
    metapop_macrodiv(options)
  }
  
  
}

if(!dont_viz){
  print("Beginning visualizations...")
  
  metapop_preproc_sum(options)
  
  if(!dont_micro){
    metapop_microdiv_viz(options)
    metapop_cbviz(options)
  } 
  
}







