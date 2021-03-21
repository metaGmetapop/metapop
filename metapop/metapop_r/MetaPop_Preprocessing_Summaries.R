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

original_files <- options[8]
original_files <- strsplit(original_files, split=",")[[1]]

original_base_names <- options[9]
original_base_names <- strsplit(original_base_names, split=",")[[1]]

setwd(directory_name)

#For my own dev use

suppressMessages(suppressWarnings(library(doParallel, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(data.table, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(ggplot2, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(cowplot, lib.loc = library_location)))

suppressMessages(suppressWarnings(library(bit64, lib.loc = library_location)))

#Original read counts is annoying... I'll have to pass this in
#original_files <- list.files(full.names = T, path = "Metapop_BAMs")


#TODO update lengths, file paths
cov_depth <- list.files(path="MetaPop/03.Breadth_and_Depth", full.names = T, pattern="breadth_and_depth")
c_d_names <- substring(cov_depth, 30, nchar(cov_depth)-22)

output_bams <- list.files(path = "MetaPop/02.Filtered_Samples", full.names = T)
output_bams <- output_bams[!grepl(".bai", output_bams)]

cl <- makeCluster(min(threads, detectCores()))

registerDoParallel(cl)

orig_read_counts <- foreach(i = original_files, .combine = c) %dopar% {
  
  as.numeric(system(paste0("samtools view -c ", i), intern = T))
  
}

stopCluster(cl)

#TODO these names won't work.

overall_summaries <- data.table(sample = original_base_names, read_counts = orig_read_counts, finish_counts = 0)

#'.BAM' already stripped.
#overall_summaries$sample = substr(overall_summaries$sample, 1, nchar(overall_summaries$sample)-4)

cl <- makeCluster(min(threads, detectCores()))

registerDoParallel(cl)

finish_read_counts <- foreach(i = output_bams, .combine = c) %dopar% {
  
  as.numeric(system(paste0("samtools view -c ", i), intern = T))
  
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


create_summary_viz <- function(depth_data, file_name, file_size, cov=min_cov, depth=min_dep, leg = color_leg){
  
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


pdf("MetaPop/12.Visualizations/preprocessing_summaries.pdf", width = 16, height = 9)

print(read_filtering)

for(i in 1:length(coverage_and_depth)){
  tmp <- create_summary_viz(coverage_and_depth[[i]], c_d_names[i], file.size(cov_depth[i]))
  print(tmp[[1]])
  print(tmp[[2]])
}

dev.off()
