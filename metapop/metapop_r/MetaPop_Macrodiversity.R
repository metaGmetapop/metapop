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

#Module specific 

coverage_cutoff <- as.numeric(options[6])
depth_cutoff <- as.numeric(options[7])

norm_file <- options[8]
whole_mag <- options[9]

if(whole_mag == "1"){
  whole_mag <- T
}else{
  whole_mag <- F
}


mag_cutoff <-as.numeric(options[10])

min_contig_length  <- as.numeric(options[11])

setwd(directory_name)

suppressMessages(suppressWarnings(library(doParallel, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(data.table, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(stringr, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(ggplot2, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(cowplot, lib.loc = library_location)))

suppressMessages(suppressWarnings(library(ggrepel, lib.loc = library_location)))

suppressMessages(suppressWarnings(library(vegan, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(compositions, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(pheatmap, lib.loc = library_location)))
suppressMessages(suppressWarnings(library(RColorBrewer, lib.loc = library_location)))

suppressMessages(suppressWarnings(library(bit64, lib.loc = library_location)))

cov_depth <- list.files(full.names = T, path = "MetaPop/03.Breadth_and_Depth")

#We just need the names, here
original_assemblies <- readLines(ref_fasta)
#Headers only.
original_assemblies <- original_assemblies[substr(original_assemblies, 1, 1) == ">"]

original_assemblies <- substr(original_assemblies, 2, nchar(original_assemblies))

original_assemblies <- unlist(lapply(original_assemblies, function(x){
  
  res  = strsplit(x, split = ' ')[[1]]
  if(length(res) > 1){
    res = res[1]
  }
  
  return(res)
  
}))

norm_header = fread(norm_file, header = T, sep = "\t")
norm_file <- fread(norm_file, header= F, sep = "\t")

if(ncol(norm_file) < 2 | ncol(norm_header) < 2){
  print("There's an error with the normalization file. MetaPop did not recognize the file as being a tab-separated, 2 column file with the first column being file names (minus final extension) and the second column being numbers.")
  print("Macrodiversity cannot proceed.")
  quit()
}

colnames(norm_file) = c("V1", "V2")
colnames(norm_header) = c("V1", "V2")

#If the norm file read without header isn't properly typed but the one read with a header is, then proceed with the header file.
if(!is.numeric(norm_file$V2) & is.numeric(norm_header$V2)){
  norm_file = norm_header
}

#Figures out the normalization factors on each sample.
norm_file$factor <- norm_file$V2/max(norm_file$V2, na.rm = T)

#Long format is generally better for operations like this. It's not a matrix yet, but you can ignore that
abundance_matrix <- CJ(original_assemblies, unlist(lapply(cov_depth, function(x){
  
  strsplit(strsplit(x, split = "MetaPop/03.Breadth_and_Depth/")[[1]][2], split = "_breadth_and_depth.tsv")[[1]][1]
  
})))


colnames(abundance_matrix) = c("V1", "V2")


#Initialize abundance counts at zero.
abundance_matrix$abundance = 0

#This both creates and normalizes the abundance table, with the caveats that only things with >= min. coverage or >= min bp covered are considered with depth >= 1
for(i in cov_depth){
  
  match_name <- strsplit(strsplit(i, split = "MetaPop/03.Breadth_and_Depth/")[[1]][2], split = "_breadth_and_depth.tsv")[[1]][1]
  
  #regex match here
  
  tmp <- fread(i, sep = "\t", header = F)
  
  tmp$V4[!(tmp$V3 >= coverage_cutoff | (tmp$V2*(tmp$V3/100))>=min_contig_length)] <- 0
  
  
  #This looks odd, but gets the rows which match the sample, then matches the contigs from the depths and adds in the values to the correct rows.
  abundance_matrix[V2==match_name,][match(tmp$V1, V1)]$abundance <- tmp$V4
  
}


#remove 

#Transform the abundance matrix into an actual matrix format. Rows are contigs, columns are samples.
raw_abundance_matrix <- dcast(abundance_matrix, V1~V2, value.var = "abundance")

#Goes ahead and outputs the table. This is what a person looking to send this table would want to send.
fwrite(raw_abundance_matrix, "MetaPop/11.Macrodiversity/raw_abundances_table.tsv", sep = "\t")

rm(raw_abundance_matrix)

#print(norm_file$factor[match(abundance_matrix$V2, norm_file$V1)])

#normalize the abundance matrix
abundance_matrix[, abundance := (abundance / norm_file$factor[match(V2, norm_file$V1)] ) ]

abundance_matrix <- dcast(abundance_matrix, V1~V2, value.var = "abundance")
rn = abundance_matrix$V1

#print(abundance_matrix)

fwrite(abundance_matrix, "MetaPop/11.Macrodiversity/normalized_abundances_table.tsv", sep = "\t")

#Technically the first column was row names. Removes it.
abundance_matrix[,V1:=NULL]

#Coverts to the R representation of a matrix instead of a data frame.
abundance_matrix <- as.matrix(abundance_matrix)

#Names the rows with their contig name.
rownames(abundance_matrix) = rn


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
write.table(alpha_diversity, file = "MetaPop/11.Macrodiversity/Alpha_diveristy_stats.tsv", sep = "\t", row.names = TRUE, col.names=NA)

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

try({

#Make pdf of alpha diversity scatterplots
pdf("MetaPop/12.Visualizations/alpha_diversity_scatterplots.pdf", width = 16, height = 16/3)

print(Richness_plot)
print(Shannons_H_plot)
print(Simpson_plot)
print(Peilous_J_plot)

print(Chao1_plot)
print(ACE_plot)
print(InvSimpson_plot)
print(Fisher_plot)

dev.off()

})


#Calculate beta-diversity metrics
clr_transformation <- clr(t(abundance_matrix))
clr_euc <-vegdist(clr_transformation, method="euclidean") 
clr_euc <- as.matrix(clr_euc)
write.table(clr_euc, file = "MetaPop/11.Macrodiversity/Beta_diveristy_clr_transformed_euclidean_distances.tsv", sep = "\t", row.names = TRUE, col.names=NA)
bray <-vegdist(t(abundance_matrix), method="bray") 
bray <- as.matrix(bray)
write.table(bray, file = "MetaPop/11.Macrodiversity/Beta_diveristy_bray_distances.tsv", sep = "\t", row.names = TRUE, col.names=NA)
jaccard <-vegdist(t(abundance_matrix), method="jaccard") 
jaccard <- as.matrix(jaccard)
write.table(jaccard, file = "MetaPop/11.Macrodiversity/Beta_diveristy_jaccard_distances.tsv", sep = "\t", row.names = TRUE, col.names=NA)

if(ncol(clr_euc) < 3){
  print("There are not enough samples to produce a meaningful PCA, PCoA, or NMDS analysis. At least 3 are required.")
  print("There will be data outputs for these analyses in MetaPop/11.Macrodiversity/, but no visualizations.")
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
    geom_label_repel(data=clr_euc.components.richness,aes(x=clr_euc.components.richness$PC1, y=clr_euc.components.richness$PC2, label=rownames(clr_euc.components.richness)))+
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
    geom_label_repel(data=bray.components.richness,aes(x=bray.components.richness$PC1, y=bray.components.richness$PC2, label=rownames(bray.components.richness)))+
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
    geom_label_repel(data=jaccard.components.richness,aes(x=jaccard.components.richness$PC1, y=jaccard.components.richness$PC2, label=rownames(jaccard.components.richness)))+
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
  
  
  try({
  
  #Make pdf of PCA plots
  pdf("MetaPop/12.Visualizations/PCA_CLR.EUC_BRAY_JACCARD_plot.pdf", width = 9, height = 9)
  print(CLR_EUC_PCA_PLOT)
  print(BRAY_PCA_PLOT)
  print(JACCARD_PCA_PLOT)
  dev.off()
  
  })
  
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
    geom_label_repel(data=clr_euc.components.richness,aes(x=clr_euc.components.richness$MDS1, y=clr_euc.components.richness$MDS2, label=rownames(clr_euc.components.richness)))+
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
    geom_label_repel(data=bray.components.richness,aes(x=bray.components.richness$MDS1, y=bray.components.richness$MDS2, label=rownames(bray.components.richness)))+
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
    geom_label_repel(data=jaccard.components.richness,aes(x=jaccard.components.richness$MDS1, y=jaccard.components.richness$MDS2, label=rownames(jaccard.components.richness)))+
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
  
  try({
  
  #Make pdf of PCoA plots
  pdf("MetaPop/12.Visualizations/PCoA_CLR.EUC_BRAY_JACCARD_plot.pdf", width = 9, height = 9)
  print(CLR_EUC_PCOA_PLOT)
  print(BRAY_PCOA_PLOT)
  print(JACCARD_PCOA_PLOT)
  dev.off()
  
  })
  
  #Make NMDS plots of all beta-diversity metrics
  clr_euc.nmds<-metaMDS(clr_euc, k=2)
  clr_euc.stress <- clr_euc.nmds$stress
  clr_euc.stress= sprintf("%#.3f", clr_euc.stress)
  clr_euc.stress = paste ("stress = ",clr_euc.stress,"")
  clr_euc.points <- clr_euc.nmds$points 
  clr_euc.points.richness <- as.data.frame(cbind(clr_euc.points, alpha_diversity)) ##add alpha diversity values
  
  CLR_EUC_NMDS_PLOT <- ggplot(data=clr_euc.points.richness)+
    geom_label_repel(data=clr_euc.components.richness,aes(x=clr_euc.components.richness$MDS1, y=clr_euc.components.richness$MDS2, label=rownames(clr_euc.components.richness)))+
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
    geom_label_repel(data=bray.components.richness,aes(x=bray.components.richness$MDS1, y=bray.components.richness$MDS2, label=rownames(bray.components.richness)))+
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
    geom_label_repel(data=jaccard.components.richness,aes(x=jaccard.components.richness$MDS1, y=jaccard.components.richness$MDS2, label=rownames(jaccard.components.richness)))+
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
  
  try({
  #Make pdf of NMDS plots
  pdf("MetaPop/12.Visualizations/NMDS_CLR.EUC_BRAY_JACCARD_plot.pdf", width = 9, height = 9)
  par(pty="s")
  print(CLR_EUC_NMDS_PLOT)
  print(BRAY_NMDS_PLOT)
  print(JACCARD_NMDS_PLOT)
  dev.off()
  })
  
}

if(length(cov_depth) > 1){
  print("Heatmaps are very memory intensive. MetaPop will attempt to print normalized abundance heatmaps, but these may fail.")
  
  try({
    
    abundance_matrix_nozero <- abundance_matrix[ rowSums(abundance_matrix)!=0, ]
    
    pdf("MetaPop/12.Visualizations/normalized_abundances_heatmap.pdf")
    
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
    
    pdf("MetaPop/12.Visualizations/normalized_abundances_heatmap_75_quantile_removed.pdf")
    
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
  
}else{
  print("Cannot make a heatmap with only one sample.")
}
