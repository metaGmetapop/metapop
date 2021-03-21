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

plot_all <- options[9]

if(plot_all == "1"){
  plot_all <- T
}else{
  plot_all <- F
}

snp_scale <- options[10]

setwd(directory_name)

suppressMessages(suppressWarnings(library(Biostrings, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(data.table, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(ggplot2, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(gggenes, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(cowplot, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(doParallel, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(RColorBrewer, lib.loc = library_location)))

# SNP positions in codons code

codon_pos_count <- fread("MetaPop/10.Microdiversity/global_codon_position_summary.tsv", sep = "\t")

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




pdf("MetaPop/12.Visualizations/Third_Pos_SNP_summary.pdf", 11, 11)

print(p)

print(p_all)

dev.off()

#FST
if(file.exists("MetaPop/10.Microdiversity/fixation_index.tsv")){
  
fixation_data <- fread("MetaPop/10.Microdiversity/fixation_index.tsv", sep = "\t")
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

pdf("MetaPop/12.Visualizations/fst_genome_heatmap_plots.pdf", width = 17, height = 11)

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

}
#main genomes


if(snp_scale=="local" | snp_scale == "both"){
  
  print("Creating local scale SNP plots...")

  genes <- readDNAStringSet(ref_genes)
  
  s <- strsplit(names(genes), "[# \t]+") # split names by tab/space
  genes <- data.table(matrix(unlist(s), ncol=5, byrow=T))
  
  names(genes)[1:4] = c("contig_gene", "start", "end", "OC")
  
  genes$start <- as.numeric((genes$start))
  genes$end <- as.numeric((genes$end))
  
  genes <- genes[,-5]
  
  genes$parent_contig <- gsub("_\\d+$", "", genes$contig_gene)
  
gene_microdiv <- fread("MetaPop/10.Microdiversity/local_gene_microdiversity.tsv", sep = "\t")

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


fastaLengths <- readDNAStringSet(ref_fasta)
fastaLengths <- data.table(contig = names(fastaLengths), num_bases = lengths(fastaLengths))
fastaLengths$contig <- unlist(lapply(fastaLengths$contig, function(x){return(strsplit(x, split=" ")[[1]][1])}))

fastaLengths <- fastaLengths[ contig %in% names(gene_microdiv),]

fastaLengths[, bigness := round(log10(num_bases))]

setkey(fastaLengths, "bigness")


unique_contigs <- unique_contigs[fastaLengths$contig]


depth_info <- list.files(path = "MetaPop/04.Depth_per_Pos", full.names = T)
depth_names <- substr(depth_info, 26, nchar(depth_info)-20)


cl <- makeCluster(min(threads, detectCores()))

clusterExport(cl, varlist = c("library_location", "depth_info"), envir = environment())
clusterEvalQ(cl, suppressMessages(suppressWarnings(library(data.table, lib.loc = library_location))))

registerDoParallel(cl)

depth_bins <- foreach(x = depth_info) %dopar% {
  
  tmp <- fread(x, sep = "\t", header = F)
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

clusterExport(cl, varlist = c("library_location", "depth_info", "depth_names", "depth_bins"))

clusterEvalQ(cl, library(data.table, lib.loc = library_location))
clusterEvalQ(cl, library(ggplot2, lib.loc = library_location))
clusterEvalQ(cl, library(gggenes, lib.loc = library_location))
clusterEvalQ(cl, library(cowplot, lib.loc = library_location))


registerDoParallel(cl)

for(z in size_category$bigness){
  
  current_conts <- unlist(size_category$V1[size_category$bigness == z])
  
  groups <- (1:length(current_conts)) %/% threads
  unique_groups <- unique(groups)
  
pdf(paste0("MetaPop/12.Visualizations/local_contigs_bp_cat_",as.character(format(10^z, scientific = F)),"_microdiversity_viz.pdf"), width = 13 + (2^(z+1)), height = 11)
  
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
  
  
  genes <- readDNAStringSet(ref_genes)
  
  s <- strsplit(names(genes), "[# \t]+") # split names by tab/space
  genes <- data.table(matrix(unlist(s), ncol=5, byrow=T))
  
  names(genes)[1:4] = c("contig_gene", "start", "end", "OC")
  
  genes$start <- as.numeric((genes$start))
  genes$end <- as.numeric((genes$end))
  
  genes <- genes[,-5]
  
  genes$parent_contig <- gsub("_\\d+$", "", genes$contig_gene)
    
  gene_microdiv <- fread("MetaPop/10.Microdiversity/global_gene_microdiversity.tsv", sep = "\t")
  
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
  
  
  fastaLengths <- readDNAStringSet(ref_fasta)
  fastaLengths <- data.table(contig = names(fastaLengths), num_bases = lengths(fastaLengths))
  fastaLengths$contig <- unlist(lapply(fastaLengths$contig, function(x){return(strsplit(x, split=" ")[[1]][1])}))
  
  fastaLengths <- fastaLengths[ contig %in% names(gene_microdiv),]
  
  fastaLengths[, bigness := round(log10(num_bases))]
  
  setkey(fastaLengths, "bigness")
  
  
  unique_contigs <- unique_contigs[fastaLengths$contig]
  
  
  depth_info <- list.files(path = "MetaPop/03.Breadth_and_Depth", full.names = T)
  depth_names <- substr(depth_info, 30, nchar(depth_info)-22)
  
  
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
  
  clusterExport(cl, varlist = c("library_location", "groups", "unique_groups", "depth_info", "depth_names", "depth_bins"))
  
  clusterEvalQ(cl, library(data.table, lib.loc = library_location))
  clusterEvalQ(cl, library(ggplot2, lib.loc = library_location))
  clusterEvalQ(cl, library(gggenes, lib.loc = library_location))
  clusterEvalQ(cl, library(cowplot, lib.loc = library_location))
  
  
  registerDoParallel(cl)
  
  for(z in size_category$bigness){
    
    current_conts <- unlist(size_category$V1[size_category$bigness == z])
    
    groups <- (1:length(current_conts)) %/% threads
    unique_groups <- unique(groups)
    
    pdf(paste0("MetaPop/12.Visualizations/global_contigs_bp_cat_",as.character(format(10^z, scientific = F)),"_microdiversity_viz.pdf"), width = 13 + (2^(z+1)), height = 11)
    
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

