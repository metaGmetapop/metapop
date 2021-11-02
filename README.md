# MetaPop #
A pipeline for the macro- and micro-diversity analyses and visualization of metagenomic-derived populations<br/><br/>


Description
---------------

MetaPop is a pipeline designed to facilitate the processing of sets of short read data mapped to reference genomes with the twin aims of calculating sample-level diversity metrics such as abundance, population diversity, and similarity across multiple samples, and assessing within-species diversity through the assessment of nucleotide polymorphisms and amino acid substitutions. To further facilitate understanding, the pipeline also produces graphical summaries of its results.

Although it can be run in a single command, MetaPop is divided into 4 modules: (1) Preprocessing, (2) Microdiversity, (3) Macrodiversity, and (4) Visualization. 

The preprocessing module is mandatory for the other modules, but also independently allows for a user to filter their reads for length, quality of alignment as measured through percent identity, and to calculate the depth and breadth of coverage on genomes within each sample. More details on these processes  are available in section 4 of this document. 

The macrodiversity module will calculate normalized abundances of species across all of the samples in a dataset and will calculate and visualize the results of a wide range of diversity indices. Included indices are available in section 4 of this document.

The microdiversity module identifies single nucleotide polymorphism (SNP) loci and performs quality assessments of these loci to ensure their accuracy. Where variant loci occur on predicted genes, MetaPop assigns the loci to all genes on which they fall, and updates the base calls of both the original reference genomes and gene calls to an all-sample consensus at each variant locus. Next, codons containing multiple SNP loci are assessed for co-occurrence to determine amino acid behavior more accurately, with more detail available in section 4. Finally, the result of these preparatory steps are used to calculate within-genome diversity measures per sample, details of which are available in section 4.

The visualization module contains a component designed to summarize preprocessing and a component designed to visualize microdiversity results. The preprocessing component includes a report of the counts of reads retained and removed in each sample, and a per-sample assessment of the depth and breadth of coverage of each genome within the sample. The microdiversity component contains plots of the depth of coverage over the length of each genome, nucleotide diversity statistics measured on each gene, and interpretations of selective pressure for each genome of each sample. It also produces a summary of codon usage bias observed in each genome to provide a starting point for assessing horizontally transferred genes.

Installation
---------------

MetaPop is a combination of python and R scripts which makes additional calls to a variety of external tools and command line utilities. It is expected to be run from within a unix environment. A unix environment can be accessed through a natively unix system, a virtual machine, or an appropriate platform such as the windows ubuntu terminal.

In order to run MetaPop, a user will require the following software:<br/>
Python 3.6+<br/>
R<br/>
Samtools<br/>
BCFTools<br/>
Prodigal prokaryotic gene prediction<br/>
Pysam<br/>

Additionally, the following R packages are required:<br/>
doParallel<br/>
data.table<br/>
stringr<br/>
ggplot2<br/>
cowplot<br/>
vegan<br/>
compositions<br/>
pheatmap<br/>
RColorBrewer<br/>
bit64<br/>
gggenes<br/>

However, MetaPop has been packaged with conda - to install it, all you must do as a user is install anaconda or miniconda (<https://docs.anaconda.com/anaconda/install/>) on your unix system, and then download the metapop conda package. This will include and install all of the needed software and libraries for you.

`conda install -c bioconda bcftools samtools prodigal pysam bioconductor-rsamtools bioconductor-biostrings`<br/>
`conda install -c conda-forge r-base`<br/>
`conda install -c r r-data.table r-ggplot2 r-rcolorbrewer`<br/>
`conda install -c conda-forge r-doparallel r-cowplot r-bit64 r-gggenes r-stringr r-vegan r-compositions r-pheatmap`<br/>
`conda install -c kgerhardt metapop`<br/>

Usage
---------------

MetaPop’s function is contained within an R script designed to be run from the command line. Usage should always start with Rscript MetaPop.R, with arguments that follow:

`Rscript MetaPop.R -dir [BAM file directory] -assem [path to reference genome file] -ct [library counts normalization file] [OPTIONS]`

Entering this command on the command line will process the mapped reads in the directory specified by -dir, with the file the reads were mapped against specified by -assem, and a file containing sample names and library read counts with -ct.

Arguments:<br/>

Mandatory arguments:<br/>
`-dir` : Directory containing mapped read files in BAM format (files do not need to be sorted and indexed)<br/>
`-assem` : Absolute path to assembled contigs<br/>
`-ct` : Absolute path to counts normalization file in tsv format (only mandatory if performing macrodiversity calculations)<br/>

To convert SAM to BAM files, you can use the following command once you activate your conda environment:<br/>
`samtools view -S -b sample.sam > sample.bam`

ct file example (exclude *.bam* suffix and table headers below) 
BAM File Names  | Number Reads or Bps in Read Library
------------- | -------------
readlibrary1  | 52000000 
readlibrary2 | 43000000
readlibrary3  | 73000000

Program options:<br/>
`-preprocess_only` : Flag indicating to filter reads for %ID, length, and depth of coverage and stop.<br/>
`-micro_only` : Flag indicating to perform only macrodiversity calculations. Assumes preprocess has been done.<br/>
`-macro_only` : Flag indicating to perform only microdiversity calculations. Assumes preprocess has been done.<br/>
`-viz_only` : Flag indicating to only produce visualizations. Assumes preprocess has been done.<br/>
`-no_micro` Flag : indicating to skip microdiversity and only perform preprocess and macrodiversity<br/>
`-no_macro` Flag : indicating to skip macrodiversity and only perform preprocess and microdiversity<br/>
`-no_viz` Flag : indicating to not attempt to visualize results.<br/>
`-threads` INT : Number of threads to parallelize processes for. Default 4.<br/>

Preprocessing Arguments:<br/>
`-sam` : absolute path to samtools - optional if samtools is in PATH<br/>
`-bcf` : absolute path to bcftools - optional if samtools is in PATH<br/>
`-prodigal` : absolute path to prodigal gene prediction tool - optional if prodigal is in PATH or if -genes specified with existing prodigal FASTA format gene calls.<br/>
`-genes` : absolute path to prodigal FASTA format gene file for assembled contiigs. May be left absent if you want metapop to generate this file. Please note if you supply your own gene file, the header for each gene in the fastea format neeeds to be in prodigal format.<br/><br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Prodigal FASTA header format example:<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`>contig1234_1 # 2 # 181 # 1 # ID=1_1;partial=10;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.628`
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`>(contig name_gene number) # (start position in contig) # (end position in contig) # (template (1) or antisense (-1) strand) # (can put NA here)).`<br/><br/>
`-id` INT : reads below this percent identity (mismatch/alignment length) are removed. Use -global to calculate as (mismatch/read length). Default 95.<br/>
`-min` INT : reads with alignments shorter than this are removed. Default 30.<br/>
`-cov` INT : contigs with breadth of coverage (#bases covered/contig length) less than this are removed from microdiversity. Default 70.<br/>
`-dep` INT : contigs with truncated average depth of coverage (mean of the 10th - 90th percentile depths of coverage) less than this are removed from microdiversity. Default 10.<br/>
`-trunc` INT : sets the percentiles at which depths of coverage will be truncated for depth. Default 10.<br/>
`-global` FLAG : Flag indicating to calculate percent identity using read length instead of alignment length.<br/>
`-force_genes` FLAG: Flag indicating that the file specified with -genes doesn't exist, and requests a file with this name be generated using prodigal.<br/>

Variant Calling Arguments:<br/>
`-first` STRING : prefix of sample to be used as the SNP reference point for all other samples. A prefix is the set of characters before .bam extension in a sample, e.g. a_file.bam has prefix a_file.<br/>
`-obs` INT : Minimum number of observations of a variant allele at a particular locus for it to be called a SNP. Default 4<br/>
`-rep` INT : Minimum percent of the population a variant allele at a particular locus must represent for it to be called a SNP. Default 1<br/>
`-var_qual` INT : Minimum PHRED score for a base to be used in initial variant calling. Default 20.<br/>

Microdiversity Arguments:<br/>
`-subsamp` INT : SNP loci will be subsampled proportionallydown to this depth of coverage for microdiversity calculations. Default 10.<br/>

Macrodiversity Arguments:<br/>
`-complete_bact` FLAG : Assumes that each contig supplied is a complete bacterial genome. Lowers the threshold of detection from 70% coverage to 20% coverage.<br/>
`-min_det` INT : Minimum percent of bases covered for a contig to be considered detected. Default 70; will change if -cov is set to a different value of -complete_bact is used.<br/>
`-min_bp` INT : Minimum number of positions required to be covered for a contig to be considered detected. Default 5000.<br/>

Visualization arguments:<br/>
`-all` FLAG : Metapop will print all contigs from microdiversity results. This will likely take a long time. Default prints top 3 genomes by highest % of genes under positive selection in each sample.<br/>
`-snp_scale` [local, global, both] : Metapop will print microdiversity results using SNPs assessed at the local (per sample) or global (across all samples) levels, or both. Defaults to local.<br/>

Function Details
---------------

Several sections of MetaPop’s code perform functions which are difficult to briefly describe. In order to organize this README better, more complete descriptions of complex function behavior and their relevant user arguments are placed here instead of with their briefer descriptions in other sections.

Preprocessing Module:<br/>
The preprocessing component of MetaPop is aimed at filtering datasets for high quality reads, then identifying the genomes which are well-covered using only these reads.

The first filter applied is for alignment length. Reads which themselves are short, or which align only over small portions of their length are more likely to be spurious. MetaPop removes reads with alignments shorter than 30 bp by default, which can be modified with the -min option.

Percent identity is then calculated on the remaining reads as the number of base pairs within a read which match the reference genome divided by the length of the alignment (or the entire length of the read if the -global flag is set). By default, reads mapping at less than 95% identity will be removed, but this value can be changed with the -id argument.

Depth and breadth of coverage per genome are then calculated using only reads passing these initial two filters. Breadth of coverage is simply the percent bases in a genome covered by at least 1 read. Depth, however, is calculated as Truncated Average Depth (TAD) which is more complex.

To calculate TAD, depth is counted for every position in the genome (including zeroes), and then these counts are sorted. The bottom and top 10% of depths are removed (e.g. in a 100 bp genome, the lowest 10 and highest 10 depths would be removed) and the average is taken over the middle 80%. This approach provides for robustness against spurious read recruitment to highly conserved regions of the genome and toward genomes with low breadth of coverage. The quantiles at which reads are truncated (default 10) may be set at different values using the -trunc argument.

Macrodiversity Indices:<br/>
The macrodiversity module calculates the following indices of sample-level diversity statistics from abundance data created during preprocessing:

Alpha diversity:<br/>
-Raw abundance counts<br/>
-Abundance counts normalized across samples of different sizes<br/>
-Chao diversity index<br/>
-ACE index<br/>
-Shannon’s H<br/>
-Simpson’s index<br/>
-Inverse Simpson’s<br/>
-Fisher alpha diversity<br/>
-Peilou’s J<br/>

Beta diversity:<br/>
-Bray-Curtis dissimilarity<br/>
-Center Log Ratio Euclidean distances<br/>
-Jaccard distances<br/>

SNP Linking:<br/>
During the calling of SNPs, MetaPop assigns SNPs to their respective genes, and notes their positions within each gene in terms of the codons to which they belong and their positions in the codons. While the positions within the codons are used to ensure that gene and SNP calls are behaving as expected and producing more SNPs in the 3rd positions of codons when compared the first or second positions, MetaPop also notes when 2 or more SNPs occur on the same codon.

When two or more SNPs are located on the same codon of the same gene, MetaPop directly examines the reads in each sample to determine the precise behavior of that codon in each sample. Because even short reads are likely to fully contain individual codons, MetaPop is able to determine if the SNPs occurring on the same codon are occurring simultaneously on the same reads (indicating that the variant codon is the composite of these two SNP loci), or whether the SNPs are mutually exclusive (indicating two distinct variants), or whether they occasionally co-occur and occasional occur separately.

The amino acids resulting from the codons examined in this process are retained for the later calculation of the ratio of synonymous and non-synonymous mutations. This renders these calculations reflective of the actual behavior of the variant sites, rather than approximations.

Microdiversity Indices:<br/>
The microdiversity module of MetaPop utilizes the SNP calls and linked SNP data generated earlier to quantify nucleotide diversity within each genome and to assess evidence of selective pressures being applied to genes on each genome. 

Microdiversity calculations are divided into local and global scales, with the local scale representing the set of variant loci identified independently within each sample, and the global scale representing the set of variant loci identified over a genome across all samples.

Microdiversity indices include:<br/>
-Pi nucleotide diversity<br/>
-Watterson’s Theta nucleotide diversity<br/>
-Tajima’s D<br/>
-pN/pS ratio of non-synonymous to synonymous mutations<br/>
-Codon usage biases of genes on each genome<br/>
-Fixation index (Fst) between samples where the same genome is observed at sufficient depth of coverage<br/>
