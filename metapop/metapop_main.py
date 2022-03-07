import sys
import re
import bisect
import sqlite3
import argparse
import platform
import importlib
import multiprocessing
import os
import zlib
import time
import tempfile
import subprocess
from array import array
from struct import unpack
from urllib.request import pathname2url
import pysam
from datetime import datetime

import metapop.metapop_helper_functions
import metapop.metapop_filter
import metapop.metapop_snp_call
import metapop.metapop_mine_reads
import metapop.metapop_fst


'''
MetaPop main script. Manages option input, flow of control, and issuing calls to functionality in other scripts. Prints logs of options.
'''
	
def gather_opts():
	parser = argparse.ArgumentParser(description=''' MetaPop ''')

	#file inputs
	parser.add_argument("--input_samples", "-i", default = "", dest = "indir", help = "Directory containing mapped reads in SAM/BAM format.")
	parser.add_argument("--reference", "-r", default = "", dest = "ref_genomes", help = "Reference genomes in FASTA format")
	parser.add_argument("--genes", "-g", dest = "genes", default = "", help = "Genes predicted with Prodigal's -d output. If not specified, MetaPop will create genes for you.")
	
	#R library link. Because this should be conda installed, there shouldn't be an issue here.
	parser.add_argument("--library", "-l", dest = "r_lib", default = "", help = "Path containing R libraries. MetaPop will attempt to use the default location of R libraries, but this will likely fail on a server.")

	#BP/read count normalization file. Auto generated if need be
	parser.add_argument("--norm", "-n", dest = "norm", default = "", help = "Base pair or read count normalization file for population diversity metrics. Default will be a count of the reads in each sample before any preprocessing.")	
	#output
	parser.add_argument("--output", "-o", dest = "outdir", default = ".", help = "Directory to place outputs.")
	
	#Filter opts
	#parser.add_argument("--genomes_are_mags", dest = "is_mag", action = "store_true", help = "The option supplied to --reference is a directory containing 1 or more multiFASTA files, each representing 1 MAG")
	
	parser.add_argument("--global", dest = "isglobal", action = "store_true", help = "Calculate pct ID as matches to reference / read length. Default is matches / alignment length")
	parser.add_argument("--id_min", dest = "min_pct_id", default = 95, help = "Reads below this Pct ID will be removed. Default 95")
	parser.add_argument("--min_len", dest = "min_length", default = 50, help = "Minimum read length. Default 50.")
	
	parser.add_argument("--min_cov", dest = "min_cov", default = 20, help = "Minimum percent of genome covered for inclusion. Default 20.")
	parser.add_argument("--min_dep", dest = "min_dep", default = 10, help = "Minimum truncated average depth (TAD) for inclusion. Defualt 10.")
	parser.add_argument("--trunc", dest = "truncation", default = 10, help = "Truncate the highest and lowest n percent of depths per genome when calculating TAD. Default 10.")

	
	#bcf/vcf opts
	parser.add_argument("--min_qual", dest = "min_qual", default = 20, help = "Minimum phred score base call quality for variant calling. Default 20.")
	parser.add_argument("--min_obs", dest = "min_obs", default = 2, help = "Minimum number of bases of the minor allele needed to be called a SNP. Default 2")
	parser.add_argument("--min_pct", dest = "min_pct", default = 1, help = "Minimum percent of the population the minor allele must represent to be called a SNP. Default 1")
	
	#SNP call opts.
	parser.add_argument("--ref_sample", dest = "reference_sample", default = "", help = "Use this sample as the point of reference for SNP calling. Reference bases will be determined by this sample's consensus base and called positions will be limited to those found in this sample. Supply the name of one of the reference BAM files for this option.")
	
	#mine reads
	parser.add_argument("--subsample_size", dest = "subsamp", default = 10, help = "Microdiversity subsampling size. Each SNP site will have up to this many bases subsampled for the calculation of microdiversity stats.")

	#macrodiv
	parser.add_argument("--whole_genomes", dest = "whole_gen", action = 'store_true', help = "Treat references as whole genomes instead of single contigs. Off by default.")
	parser.add_argument("--genome_detection_cutoff", dest = "gen_cov_min", default = 0, help = "Percent of bases that must be covered for a sequence to be considered detected for macrodiversity analyses. Use this or --minimum_bases_for_detection, not both. Default 0.")
	parser.add_argument("--minimum_bases_for_detection", dest = "min_bp", default = 5000, help = "Count of bases that must be covered for a sequence to be considered detected for macrodiversity analyses. Use this or --genome_detection_cutoff, not both. Default 5000.")
	
	#preproc_viz
	parser.add_argument("--plot_all", dest = "plot_all", action = 'store_true', help = "MetaPop normally plots the top 20 genomes with the most genes under selection in microdiversity visualizations. This flag will make MetaPop plot all genomes. Warning: slow.")
	parser.add_argument("--snp_scale", dest = "snp_scale", default = 'local', help = "Choose one of \'local\', \'global\' or \'both\'. Controls whether SNPs detected for a genome in each sample alone are shown, or SNPs for that genome across all samples.")
	
	#misc opts
	parser.add_argument("--threads", dest = "thds", default = 0, help = "How many threads to use for parallelization.")	
	
	parser.add_argument("--preprocess_only", dest = "preproc_only", action = 'store_true', help = "Perform only preprocessing (Sort reads, filter by Pct. identity and min. length, genomes by depth and breadth of coverage)")	
	parser.add_argument("--no_micro", dest = "no_mic", action = 'store_true', help = "Skip microdiversity analyses (SNP calling, codon usage bias, gene selection pressure)")	
	parser.add_argument("--no_macro", dest = "no_mac", action = 'store_true', help = "Skip macrodiversity analyses (Species richness, sample diversity)")	
	parser.add_argument("--no_viz", dest = "no_viz", action = 'store_true',   help = "Skip producing visual summaries of results.")	
	
	parser.add_argument("--viz_only", dest = "only_viz", action = 'store_true',   help = "Skip all data generation and produce visualizations. Only works if you have previously run all of MetaPop and want to make more visualizations.")	
	parser.add_argument("--skip_preproc", dest = "skip_pre", action = 'store_true',   help = "Skip preprocessing stage of MetaPop (filtering reads by pct. ID, genomes by coverage and depth). Meant for rerunning later stages of MetaPop without redoing preprocessing.")	
	parser.add_argument("--skip_snp_calling", dest = "skip_var", action = 'store_true',   help = "Skip variant calling microdiversity stage of MetaPop. Meant for rerunning microdiversity without redoing variant calling.")	

	parser.add_argument("--just_codon_bias", dest = "justCB", action = 'store_true', help = 'Use this flag alongside the --genes, --output, and --library arguments to calculate codon bias for the genes and then quit.')
	
	return parser
	
def main():

	parser = gather_opts()
	
	if len(sys.argv) == 1:
		parser.print_help(sys.stderr)
		quit()
	
	opts = parser.parse_args()
	
	
	cb_and_quit = opts.justCB
	if cb_and_quit:
		r_scripts_loc = os.path.dirname(sys.modules['metapop'].__file__) + "/metapop_r/"
		output_directory_base = opts.outdir
		if output_directory_base.endswith("/"):
			output_directory_base = output_directory_base[:-1]
			
		if output_directory_base != ".":
			if not os.path.exists(output_directory_base):
				os.mkdir(output_directory_base)
			if not os.path.exists(output_directory_base+"/MetaPop"):
				os.mkdir(output_directory_base+"/MetaPop")
		elif not os.path.exists("MetaPop"):
			os.mkdir("MetaPop")
			
		reference_genes = opts.genes
		if reference_genes == "":
			print("You need to give a genes file to just calculate codon usage bias. This needs to be a file produced by Prodigal gene prediction using the -d output option.")
			return None
		
		rlib = opts.r_lib
		
		command = ["Rscript", r_scripts_loc + "MetaPop_Codon_Bias_Calc_Independent.R", output_directory_base, reference_genes, rlib]
		try:
			subprocess.call(command)
		except:
			print("Codon bias calculation failed.")
		return None
	
	#put the opts here.
	original_bams = opts.indir
	reference_fasta = opts.ref_genomes
	reference_genes = opts.genes
	norm_file = opts.norm
	
	if reference_fasta == "":
		print("MetaPop needs a set of reference genomes to proceed. Please provide it with --reference")
		parser.print_help(sys.stderr)
		quit()
	
	if original_bams == "":
		print("MetaPop needs aligned reads. Please specify a directory containing these reads with --input")
		parser.print_help(sys.stderr)
		quit()
	
	#This is for later development
	#treat_as_mag = opts.is_mag
	treat_as_mag = False
	
	output_directory_base = opts.outdir
	if output_directory_base.endswith("/"):
		output_directory_base = output_directory_base[:-1]
	
	is_global = opts.isglobal
	
	min_id = float(opts.min_pct_id)
	min_len = int(opts.min_length)
	
	min_cov = int(opts.min_cov)
	min_dep = int(opts.min_dep)
	trunc = float(opts.truncation)
	
	if trunc > 49.9:
		print("MetaPop cannot truncate all data! Choose a value for truncation between 0 and 49. Exiting.")
		quit()
	
	#pileup opts
	min_obs = int(opts.min_obs)
	min_pct = float(opts.min_pct)
	min_q = int(opts.min_qual)

	#SNP call opts.
	ref_samp = opts.reference_sample
	
	threads = int(opts.thds)
	
	#R stuff
	rlib = opts.r_lib
	
	#program opts
	preproc_only = opts.preproc_only
	
	no_mic = opts.no_mic
	no_mac = opts.no_mac
	no_viz = opts.no_viz
	
	only_viz = opts.only_viz
	skip_pre = opts.skip_pre
	skip_var_call = opts.skip_var
		
	#mine reads
	sub_sample_size = int(opts.subsamp)

	#macrodiv
	norm_file = opts.norm
	
	whole_genome = opts.whole_gen
	if whole_genome:
		whole_genome = "1"
	else:
		whole_genome = "0"

	genome_cov_cutoff = int(opts.gen_cov_min)
	min_bp_len = int(opts.min_bp)

	#preproc_viz
	plot_all = opts.plot_all
	if plot_all:
		plot_all = "1"
	else:
		plot_all = "0"
	snp_scale = opts.snp_scale
	
	#end opts
	
	time_format = "%d/%m/%Y %H:%M:%S"
	
	#Get the maximum number of cores and runs with it by default.
	if threads == 0:
		threads = len(os.sched_getaffinity(0))

	print("Using:", str(threads), "threads.")
	
	#todo move these to helper functions
	#Setup directories - each steps sets itself up, but this is what main has to do.

	metapop.metapop_helper_functions.main_dir_prep(output_directory_base)
	
	if not no_mac:
		if norm_file == "":
			norm_file = metapop.metapop_helper_functions.produce_default_normalization_file(output_directory_base, original_bams, threads)
		else:
			#Produce non-relative path for accessing the normalization file
			norm_file = os.path.abspath(norm_file)
			
	joined_fastas = metapop.metapop_helper_functions.multi_to_single_fasta(output_directory_base, reference_fasta)
	
	if reference_genes == "":
		reference_genes = metapop.metapop_helper_functions.gene_calls(joined_fastas, output_directory_base)
	else:
		#Produce non-relative path for accessing the genes file
		reference_genes = os.path.abspath(reference_genes)
		
	#groups each input file with its contained contigs
	mag_contig_dict, mag_length_dict = metapop.metapop_helper_functions.create_mag_log(reference_fasta, output_directory_base+"/MetaPop/00.Log_and_Parameters/mags.txt", treat_as_mag)

	
	#fill out logs/params here for the R components
	param = open(output_directory_base+"/MetaPop/00.Log_and_Parameters/run_settings.tsv", "w")

	#Header
	print("parameter", "setting", sep = "\t", file = param)
	#Args
	print("Directory", os.path.normpath(output_directory_base+"/MetaPop/"), sep = "\t", file = param)
	print("Samtools Location", "", sep = "\t", file = param)
	print("BCFTools Location", "", sep = "\t", file = param)
	print("Library Location", "", sep = "\t", file = param)
	print("Assembly", joined_fastas, sep = "\t", file = param)
	print("Genes", reference_genes, sep = "\t", file = param)
	print("Normalization File", norm_file, sep = "\t", file = param)
	#print("Using MAGs", treat_as_mag, sep = "\t", file = param)
	print("ID Cutoff", min_id, sep = "\t", file = param)
	print("Min. Read Length", min_len, sep = "\t", file = param)
	print("Coverage", min_cov, sep = "\t", file = param)
	print("Depth", min_dep, sep = "\t", file = param)
	print("Truncation", trunc, sep = "\t", file = param)
	print("Variant Base Call Cutoff", min_q, sep = "\t", file = param)
	print("Subsample Size", sub_sample_size, sep = "\t", file = param)
	print("BP for Detect", min_bp_len, sep = "\t", file = param)
	print("Macrodiversity with whole genomes", whole_genome, sep = "\t", file = param)
	print("Macrodiversity percent coverage minimum", genome_cov_cutoff, sep = "\t", file = param)
	print("All Genomes Plotted", plot_all, sep = "\t", file = param)
	print("SNP Scale", snp_scale, sep = "\t", file = param)	
	print("Threads", "setting", sep = "\t", file = param)

	param.close()
	
	timer = datetime.now()
	printable_time = timer.strftime(time_format)
	print("MetaPop Started at:", printable_time)
	print("")
	
	#Skip all non-viz blocks.
	if not only_viz:
	
		#Structure of the basic commands for filtering. Passing them as arrays makes it easier to use in the filtering commands.
		filter_command_base = ["input_bam", output_directory_base, min_len, min_id, is_global]
		bdb = [min_cov, min_dep, treat_as_mag, trunc]
		
		if not skip_pre:
			metapop.metapop_filter.filt(original_bams, filter_command_base, bdb, threads, mag_contig_dict, mag_length_dict, joined_fastas)
			
			if preproc_only:
				print("Thanks for using MetaPop.")
				quit()
		
		r_scripts_loc = os.path.dirname(sys.modules['metapop'].__file__) + "/metapop_r/"
		#print("R scripts located at", r_scripts_loc)
		
		
		micro_done = False
		if not no_mic:
			micro_done = True
			
			if not skip_var_call:
				#update to base corrected genomes, genes
				joined_fastas, reference_genes = metapop.metapop_snp_call.call_variant_positions(output_directory_base, joined_fastas, min_obs, min_q, min_pct, threads, ref_samp)
				#do codon bias
				cb_call = ["Rscript", r_scripts_loc + "MetaPop_Codon_Bias_Calc.R", output_directory_base, rlib]
				subprocess.call(cb_call)
			else:
				#The R code moves into the output directory, so we have to remove the directory prefix in the path.
				joined_fastas, reference_genes = joined_fastas[(len(output_directory_base)+1):], reference_genes[(len(output_directory_base)+1):]
			
			metapop.metapop_helper_functions.micro_prep(output_directory_base)

			timer = datetime.now()
			printable_time = timer.strftime(time_format)
			print("Linking SNPs starting at:", printable_time+"...", end = "", flush = True)

			#Python version
			linked_file = metapop.metapop_mine_reads.do_mine_reads(output_directory_base, threads)
			
			fishers_call = ["Rscript", r_scripts_loc + "MetaPop_Fisher_Exact.R", linked_file, rlib]
			
			subprocess.call(fishers_call)
			
			#subprocess.call(mine_reads_call, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
			print("done!")
			
			microdiv_call = ["Rscript", r_scripts_loc + "MetaPop_Microdiversity.R", output_directory_base, str(threads), rlib, joined_fastas, reference_genes, str(min_cov), str(min_dep), str(sub_sample_size)]
			
			timer = datetime.now()
			printable_time = timer.strftime(time_format)
			print("Calculating Microdiversity starting at:", printable_time+"...", flush = True)
			subprocess.call(microdiv_call)
			#subprocess.call(microdiv_call, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
			print("Done!")
			
			#microdiv_file, lengths_file, out_dir, threads
			print("Calculating fixation index FST starting at:", printable_time+"...", end = "", flush = True)
			microdiv_file = os.path.normpath(output_directory_base + "/MetaPop/10.Microdiversity/global_raw_microdiversity_data_snp_loci_only.tsv")
			lengths_file = os.path.normpath(output_directory_base + "/MetaPop/10.Microdiversity/global_contig_microdiversity.tsv")
			print("Done!")
			metapop.metapop_fst.perform_fst(microdiv_file, lengths_file, output_directory_base, threads)

		#shorten the names for R unless micro was done, which does this already
		if not micro_done:
			joined_fastas, reference_genes = joined_fastas[(len(output_directory_base)+1):], reference_genes[(len(output_directory_base)+1):]
			
		if not no_mac:
			metapop.metapop_helper_functions.macro_prep(output_directory_base)
			macrodiv_call = ["Rscript", r_scripts_loc + "MetaPop_Macrodiversity.R", output_directory_base, str(threads), rlib, joined_fastas, reference_genes, str(min_cov), str(min_dep), norm_file, whole_genome, str(genome_cov_cutoff), str(min_bp_len)]
			
			timer = datetime.now()
			printable_time = timer.strftime(time_format)
			print("Calculating and visualizing Macrodiversity starting at:", printable_time+"...", end = "", flush = True)
			subprocess.call(macrodiv_call)
			#subprocess.call(macrodiv_call, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
			print("Done!")

	if not no_viz:
		metapop.metapop_helper_functions.viz_prep(output_directory_base)

		samples = os.listdir(original_bams)
		for i in range(0, len(samples)):
			samples[i] = os.path.normpath(original_bams + "/" + samples[i])
			
		names = metapop.metapop_helper_functions.get_base_names(samples)
			
		#passable string to preproc.
		samples = ",".join(samples)
		
		names = ",".join(names)

		preproc_sums_call = ["Rscript", r_scripts_loc + "MetaPop_Preprocessing_Summaries.R", output_directory_base, str(threads), rlib, joined_fastas, reference_genes, str(min_cov), str(min_dep), original_bams, names]
		microdiv_viz = ["Rscript", r_scripts_loc + "MetaPop_Microdiversity_Visualizations.R", output_directory_base, str(threads), rlib, joined_fastas, reference_genes, str(min_cov), str(min_dep), str(sub_sample_size), plot_all, snp_scale]
		cb_viz = ["Rscript", r_scripts_loc + "MetaPop_Codon_Bias_Viz.R", output_directory_base, str(threads), rlib, joined_fastas, reference_genes, str(min_cov), str(min_dep)]
		
		timer = datetime.now()
		printable_time = timer.strftime(time_format)
		print("Visualizing preprocessing starting at:", printable_time+"...", end = "", flush = True)
		subprocess.call(preproc_sums_call)
		#subprocess.call(preproc_sums_call, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		print("Done!")
		
		
		if not no_mic:
			timer = datetime.now()
			printable_time = timer.strftime(time_format)
			print("Visualizing Microdiversity starting at:", printable_time+"...", end = "", flush = True)
			subprocess.call(microdiv_viz)
			#subprocess.call(microdiv_viz, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
			print("Done!")
			
			timer = datetime.now()
			printable_time = timer.strftime(time_format)
			print("Visualizing Codon Bias starting at:", printable_time+"...", end = "", flush = True)
			subprocess.call(cb_viz)
			#subprocess.call(cb_viz, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
			print("Done!")
	
	if os.path.exists(os.path.normpath(output_directory_base + "/Rplots.pdf")):
		os.remove(os.path.normpath(output_directory_base + "/Rplots.pdf"))
	
	print("")
	print("Thanks for using MetaPop.")
	
	return None
	
	
	
if __name__ == "__main__":
	main()