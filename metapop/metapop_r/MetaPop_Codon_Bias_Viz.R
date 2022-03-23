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


setwd(directory_name)

suppressMessages(suppressWarnings(library(doParallel, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(data.table, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(ggplot2, lib.loc = library_location)))
#suppressMessages(suppressWarnings(library(Biostrings, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(cowplot, lib.loc = library_location)))



passing_contigs <- list.files(path = "MetaPop/03.Breadth_and_Depth", pattern = ".tsv", full.names = T)


cl <- makeCluster(min(detectCores(), threads, length(passing_contigs)))

clusterExport(cl, varlist=c("passing_contigs", "min_cov", "min_dep", "library_location"))
clusterEvalQ(cl, suppressMessages(suppressWarnings(library(data.table, lib.loc = library_location))))

registerDoParallel(cl)

passing_contigs <- unique(foreach(i=passing_contigs, .combine=c) %dopar% {
  
  tmp <- fread(i, sep="\t", header = F)
  
  return(tmp$V1[tmp$V3>=min_cov & tmp$V4 >= min_dep])
  
})

stopCluster(cl)

#Already per contig
codon_bias_iqr <- fread(list.files(full.names = T, path = "MetaPop/08.Codon_Bias", pattern="gene_IQR_and_mean.tsv"), sep = "\t")
codon_bias_iqr <- codon_bias_iqr[codon_bias_iqr$parent_contig %in% passing_contigs,]

parse_genes = function(file){
  
  genes <- readLines(file)
  
  headers = substr(genes, 1, 1) == ">"
  
  genes = genes[headers]
  genes = substr(genes, 2, nchar(genes))
  
  s <- strsplit(genes, "[# \t]+") # split names by tab/space
  genes <- data.table(matrix(unlist(s), ncol=5, byrow=T))
  
  names(genes)[1:4] = c("contig_gene", "start", "end", "OC")
  
  genes$start <- as.numeric((genes$start))
  genes$end <- as.numeric((genes$end))
  
  genes <- genes[,-5]
  
  genes$parent_contig <- gsub("_\\d+$", "", genes$contig_gene)
  
  return(genes)
  
}


genes = parse_genes(ref_genes)


genes <- genes[genes$parent_contig %in% passing_contigs,]

codon_bias_genes <- fread(list.files(full.names = T, path = "MetaPop/08.Codon_Bias", pattern="gene_euclidean_distances.tsv"), sep = "\t")
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


cl <- makeCluster(min(threads, detectCores()))
clusterExport(cl, varlist = c("library_location", "groups", "unique_groups", "color_leg"), envir = environment())
clusterEvalQ(cl, expr = suppressMessages(suppressWarnings(library(data.table, lib.loc = library_location))))
clusterEvalQ(cl, expr = suppressMessages(suppressWarnings(library(ggplot2, lib.loc = library_location))))
clusterEvalQ(cl, expr = suppressMessages(suppressWarnings(library(cowplot, lib.loc = library_location))))

pdf("MetaPop/12.Visualizations/codon_bias_plots.pdf", height = 9, width = 12)

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

