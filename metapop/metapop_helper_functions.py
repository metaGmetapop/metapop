import sys
import os
import subprocess
import multiprocessing

def dir_to_fastas(directory):
	ref_files = []
	for file in os.listdir(directory):
		if not file.upper().endswith(".FAI") and not file.upper().endswith(".R"):
			ref_files.append(os.path.normpath(directory + "/" + file))
		
	return ref_files
	
def contig_log(reference_files):
	log_file = []
	mag_dict = {}
	length_dict = {}
	
	contig = ""
	
	for file in reference_files:
		
		seqlen = 0
				
		fh = open(file)
		
		for line in fh:
			if line.startswith(">"):
				if seqlen > 0 :
					
					if contig != "MissingNo.":
						log_file.append(contig + "\t" + contig + "\t" + str(seqlen))
						mag_dict[contig] = [contig]
						length_dict[contig] = seqlen
					
					#I ran into a wierd case where a seqID was absent, so there was just a ">" on an otherwise empty line. This handles that case.
					if len(line.strip()) > 1:
						contig = line[1:].strip().split()[0]
					else:
						print("Missing sequence ID found.")
						contig = "MissingNo."
						
					seqlen = 0
				else:
					if len(line.strip()) > 1:
						contig = line[1:].strip().split()[0]
					else:
						print("Missing sequence ID found.")
						contig = "MissingNo."
			else:
				seqlen += len(line.strip())
		
		#Last one doesn't normally happen
		if contig != "MissingNo.":
			log_file.append(contig + "\t" + contig + "\t" + str(seqlen))
			mag_dict[contig] = [contig]
			length_dict[contig] = seqlen
		
		fh.close()
		
	return mag_dict, length_dict, log_file
	
def MAG_log(reference_files):
	log_file = []
	mag_dict = {}
	length_dict = {}
	
	contig = ""
	
	for file in reference_files:
		
		mag_dict[file] = []
		
		seqlen = 0
		
		fh = open(file)
		
		for line in fh:
			if line.startswith(">"):
				if seqlen > 0 :
					if contig != "MissingNo.":
						log_file.append(file + "\t" + contig + "\t" + str(seqlen))
						mag_dict[file].append(contig)
						length_dict[contig] = seqlen
					
					if len(line.strip()) > 1:
						contig = line[1:].strip().split()[0]
					else:
						print("Missing sequence ID found.")
						contig = "MissingNo."
					seqlen = 0
				else:
					if len(line.strip()) > 1:
						contig = line[1:].strip().split()[0]
					else:
						print("Missing sequence ID found.")
						contig = "MissingNo."
			else:
				seqlen += len(line.strip())
		
		#Last one doesn't normally happen
		if contig != "MissingNo.":
			log_file.append(file + "\t" + contig + "\t" + str(seqlen))
			mag_dict[file].append(contig)
			length_dict[contig] = seqlen
			
		#reads input are sorted by contig, so I also want the contigs in each mag sorted for later
		
		fh.close()
		
	return mag_dict, length_dict, log_file
	
def create_mag_log(fasta_dir, output_file, is_mags):
	ref_files = dir_to_fastas(fasta_dir)
	
	if is_mags:
		mag_log, length_log, log_file = MAG_log(ref_files)
	else:
		mag_log, length_log, log_file = contig_log(ref_files)
		
	#log = open(output_file, "w")
	#for line in log_file:
	#	print(line, file = log)
	#log.close()
	
	return mag_log, length_log
	
	
def main_dir_prep(output_directory_base):
	if output_directory_base != ".":
		if not os.path.exists(output_directory_base):
			os.mkdir(output_directory_base)
		if not os.path.exists(output_directory_base+"/MetaPop"):
			os.mkdir(output_directory_base+"/MetaPop")
	elif not os.path.exists("MetaPop"):
		os.mkdir("MetaPop")
	
	if not os.path.exists(output_directory_base+"/MetaPop/00.Log_and_Parameters"):
		os.mkdir(output_directory_base+"/MetaPop/00.Log_and_Parameters")
	if not os.path.exists(output_directory_base+"/MetaPop/01.Genomes_and_Genes"):
		os.mkdir(output_directory_base+"/MetaPop/01.Genomes_and_Genes")


def macro_prep(output_directory_base):
	if not os.path.exists(output_directory_base+"/MetaPop/11.Macrodiversity"):
		os.mkdir(output_directory_base+"/MetaPop/11.Macrodiversity")


def micro_prep(output_directory_base):
	if not os.path.exists(output_directory_base+"/MetaPop/09.Linked_SNPs"):
		os.mkdir(output_directory_base+"/MetaPop/09.Linked_SNPs")
	if not os.path.exists(output_directory_base+"/MetaPop/10.Microdiversity"):
		os.mkdir(output_directory_base+"/MetaPop/10.Microdiversity")



def viz_prep(output_directory_base):
	if not os.path.exists(output_directory_base+"/MetaPop/12.Visualizations"):
		os.mkdir(output_directory_base+"/MetaPop/12.Visualizations")

	
def multi_to_single_fasta(directory, fastas):
	
	if not os.path.exists(directory + "/MetaPop/01.Genomes_and_Genes/all_genomes.fasta"):
		out = open(directory + "/MetaPop/01.Genomes_and_Genes/all_genomes.fasta", "w")
		
		fasta_list = dir_to_fastas(fastas)
		
		#Characters that can form valid FASTA sequences.
		accept_lines = set(["A", "T", "C", "G", "N", "a", "t", "c", "g", "n"])
		
		print("Combining files:")
		for file in fasta_list:
			print(file)
		
		for file in fasta_list:
			fh = open(file)
			
			for line in fh:
				#Deals with empty lines, empty names
				clean_line = line.strip()
				if len(clean_line) > 1:
					if clean_line.startswith(">") or set(clean_line) <= accept_lines:
						print(line, end = "", file = out)
			
			fh.close()
			
		out.close()
		
		index_cmd = ["samtools", "faidx", directory + "/MetaPop/01.Genomes_and_Genes/all_genomes.fasta"]
		subprocess.call(index_cmd)
	else:
		print("Combined genomes file already found.")
		
	return directory + "/MetaPop/01.Genomes_and_Genes/all_genomes.fasta"
		
		
def gene_calls(joint_fasta, output_directory_base):
	ref_base = joint_fasta.split("/")
	ref_base = ref_base[len(ref_base)-1]
	ref_base = os.path.splitext(os.path.basename(os.path.normpath(ref_base)))[0]
				
	prod_cmd = ["prodigal", "-q", "-p", "meta", "-i", joint_fasta, "-d", output_directory_base+"/MetaPop/01.Genomes_and_Genes/"+ref_base+"_genes.fasta", "-o", output_directory_base+"/MetaPop/01.Genomes_and_Genes/temp.txt"]
	#do prodigal first
	if not os.path.exists(output_directory_base+"/MetaPop/01.Genomes_and_Genes/"+ref_base+"_genes.fasta"):
		print("Genes will be predicted with prodigal and placed in Genomes and Genes.")
		subprocess.call(prod_cmd)
		os.remove(output_directory_base+"/MetaPop/01.Genomes_and_Genes/temp.txt")
	else:
		print("Automatically predicted genes found. Skipping Prediction.")
		
	reference_genes = os.path.normpath(output_directory_base+"/MetaPop/01.Genomes_and_Genes/"+ref_base+"_genes.fasta")
	
	return(reference_genes)

def count_reads(file):

	samtools_count = ["samtools", "view", "-c", file]
	link = subprocess.Popen(samtools_count, stdout = subprocess.PIPE)
	count = int(link.stdout.readline())
	return count
		
def produce_default_normalization_file(directory_base, sample_path, threads):

	norm_file_path = os.path.normpath(directory_base + "/MetaPop/00.Log_and_Parameters/normalized_counts.txt")
	
	if not os.path.exists(norm_file_path):

		print("Producing read counts for normalization... ", end = "", flush = True)
		#Get file read counts
		samples = os.listdir(sample_path)
		for i in range(0, len(samples)):
			samples[i] = os.path.normpath(sample_path + "/" + samples[i])
		
		pool = multiprocessing.Pool(min(len(samples), int(threads)))
		counts = pool.map(count_reads, samples)
		pool.close()
		pool.join()
		
		#get file base names
		base_names = []
		for sample in samples:
			base_name = os.path.basename(os.path.normpath(sample))
			if base_name.endswith(".bam") or base_name.endswith(".BAM"):
				base_name = base_name[:-4]
			base_names.append(base_name)
			
		norm_file = open(norm_file_path, "w")
		
		for i in range(0, len(base_names)):
			print(base_names[i], counts[i], sep = "\t", file = norm_file)
		
		norm_file.close()
	
		print("done!")
		
	else:
		print("Automatically generated normalization file already found!")
	
	relative_norm_file_path = "MetaPop/00.Log_and_Parameters/normalized_counts.txt"
	
	return relative_norm_file_path
	
	
def get_base_names(bams):
	cleaned_bams = []
	for file in bams:
		if file.upper().endswith(".BAM") or file.upper().endswith(".SAM"):
			cleaned_bams.append(file)
	
	base_names = []
	
	for sample in cleaned_bams:
		base_name = os.path.basename(os.path.normpath(sample))
		if base_name.endswith(".bam") or base_name.endswith(".BAM"):
			base_name = base_name[:-4]
		
		base_names.append(base_name)
		
	return base_names
	