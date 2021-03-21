import pysam
import sys
import multiprocessing
import datetime
import re
import os
import zlib
from datetime import datetime
import tempfile
import subprocess
from array import array
from struct import unpack
from urllib.request import pathname2url
import math
import bisect

def prep(directory):
	if not os.path.exists(directory + "/MetaPop/05.Variant_Calls"):
		os.mkdir(directory + "/MetaPop/05.Variant_Calls")
	if not os.path.exists(directory + "/MetaPop/06.SNP_Counts"):
		os.mkdir(directory + "/MetaPop/06.SNP_Counts")
	if not os.path.exists(directory + "/MetaPop/07.Cleaned_SNPs"):
		os.mkdir(directory + "/MetaPop/07.Cleaned_SNPs")
	if not os.path.exists(directory + "/MetaPop/08.Codon_Bias"):
		os.mkdir(directory + "/MetaPop/08.Codon_Bias")

def gather_files(directory):
	reads = os.listdir(directory + "/MetaPop/02.Filtered_Samples")
	cleaned_reads = []
	for sample in reads:
		if not sample.endswith(".bai"):
			cleaned_reads.append(sample)
			
	return cleaned_reads
	
def output_file_names(directory, reads):
	output_files = []
	
	for sample in reads:
		output_files.append(directory + "/MetaPop/05.Variant_Calls/" + sample[:-17])
		
	return(output_files)
	
def call_variant_positions(directory, ref_fasta, min_obs, min_qual, min_pct, threads, ref_sample):
	
	prep(directory)
	#base names at this point
	reads = gather_files(directory)
	#uses the base names
	outputs = output_file_names(directory, reads)
	
	reads_with_paths = []
	for sample in reads:
		reads_with_paths.append(directory + "/MetaPop/02.Filtered_Samples/" + sample)
	
	num_samps = len(reads_with_paths)
	
	time_format = "%d/%m/%Y %H:%M:%S"
	
	
	plf = open(directory+"/MetaPop/ploidy.txt", "w")
	print("*", "*", "*", "*", 1, sep = "\t", file = plf)
	plf.close()
	
	vcf_commands = []
	
	for i in range(0, num_samps):
	#for i in range(0, 1):
		vcf_commands.append([reads_with_paths[i], outputs[i], ref_fasta, directory+"/MetaPop/ploidy.txt", min_qual])
	
	pool = multiprocessing.Pool(min(threads, num_samps))
	
	timer = datetime.now()
	printable_time = timer.strftime(time_format)
	print("Calling variants starting at:", printable_time)
	
	called_snps = pool.map(create_vcfs, vcf_commands)
	
	pool.close()
	
	#Cleanup
	os.remove(directory+"/MetaPop/ploidy.txt")
		
	
	global _all_variant_positions
	
	unique_contig_pos = {}
	
	for vcf in called_snps:
		fh = open(vcf)
		for line in fh:
			segs = line.strip().split()
			contig = segs[0]
			pos = int(segs[1])
			if contig not in unique_contig_pos:
				unique_contig_pos[contig] = [pos]
			else:
				if pos not in unique_contig_pos[contig]:
					unique_contig_pos[contig].append(pos)
					
	_all_variant_positions = unique_contig_pos
	
	##########
	variant_commands = []
	
	for i in range(0, num_samps):
		current_command = [reads_with_paths[i], outputs[i], min_qual, min_obs, min_pct, called_snps[i], ref_fasta]
		variant_commands.append(current_command)
	
	timer = datetime.now()
	printable_time = timer.strftime(time_format)
	print("Refining SNP calls starting at:", printable_time)

	pool = multiprocessing.Pool(min(threads, num_samps))
	pool.map(pileup, variant_commands)
	pool.close()

	consensus, snp_inputs = choose_consensus(directory, ref_sample)

	new_genomes, new_genes = associate_genes(directory, consensus, snp_inputs, min_obs, min_pct)
	
	timer = datetime.now()
	printable_time = timer.strftime(time_format)
	print("MetaPop SNP refinement finished at:", printable_time)
	
	return new_genomes, new_genes

def create_vcfs(commands):
	#reads, output_base_names, fasta, ploidy, min_q
	read = commands[0]
	base_name = commands[1]
	fasta = commands[2]
	ploidy = commands[3]
	qualarg = '"QUAL>'+str(commands[4])+'"'
	
	time_format = "%d/%m/%Y %H:%M:%S"
	
	printable_name = os.path.basename(base_name)
	
	timer = datetime.now()
	printable_time = timer.strftime(time_format)
	print("Variant calling for:", printable_name, "started at:", printable_time)
	
	subprocess.call(["bcftools", "mpileup", "-Ob", "-I", "-f", fasta, "-o", base_name+".mpileup.txt", read], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
	subprocess.call(["bcftools", "call", "-mv", "-Ob", "--ploidy-file", ploidy, "-o", base_name+".calls.bcf", base_name+".mpileup.txt"])
	os.remove(base_name+".mpileup.txt")
	subprocess.call("bcftools filter -i"+ qualarg +" -Ob -o "+ base_name+".filtered_calls.bcf "+ base_name+".calls.bcf", shell=True)
	os.remove(base_name+".calls.bcf")
	subprocess.call(["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", "-o", base_name+".final_called.txt", base_name+".filtered_calls.bcf"])
	os.remove(base_name+".filtered_calls.bcf")
	
	timer = datetime.now()
	printable_time = timer.strftime(time_format)
	print("Variant calling for:", printable_name, "finished at:", printable_time)
		
	return base_name+".final_called.txt"
	
def pileup(pileup_command):

	read     = pileup_command[0]
	out      = pileup_command[1]+"_filtered_snps.txt"
	
	basename = pileup_command[1]
	
	out = out.split("/")
	out[len(out)-2] = "06.SNP_Counts"
	out = '/'.join(out)
		
	min_qual = pileup_command[2]
	min_obs  = pileup_command[3]
	#Pct comes in as integer, comparison occurs as a decimal.
	min_pct  = pileup_command[4]/100
	
	vcf = pileup_command[5]
	
	reference_genomes = pileup_command[6]
	
	time_format = "%d/%m/%Y %H:%M:%S"
	
	printable_name = os.path.basename(basename)
	
	timer = datetime.now()
	printable_time = timer.strftime(time_format)
	print("SNP refinement for:", printable_name, "started at:", printable_time)

	output_pile = open(out, "w")
	counter = 0
	
	mpileup_call = ["samtools", "mpileup", "--skip-indels", "-f", reference_genomes, read]
	
	pile_stream = subprocess.Popen(mpileup_call, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
	
	for line in pile_stream.stdout:
		segs = line.strip().decode().split()
		contig = segs[0]
		pos = int(segs[1])
		if contig in _all_variant_positions:
			if pos in _all_variant_positions[contig]:
				#Just make sure.
				refbase = segs[2].upper()
				piled_bases = segs[4]
				
				#mpileup encodes reference base as . or , for main, complement strand matches
				dot_comma_ct = piled_bases.count(",") + piled_bases.count(".")
				a_ct = piled_bases.count("a") + piled_bases.count("A")
				t_ct = piled_bases.count("t") + piled_bases.count("T")
				c_ct = piled_bases.count("c") + piled_bases.count("C")
				g_ct = piled_bases.count("g") + piled_bases.count("G")
				
				if refbase == "A":
					a_ct += dot_comma_ct
				if refbase == "T":
					t_ct += dot_comma_ct
				if refbase == "C":
					c_ct += dot_comma_ct
				if refbase == "G":
					g_ct += dot_comma_ct
				
				print(contig, pos, a_ct, t_ct, c_ct, g_ct, sep = "\t", file = output_pile)
				counter += 1
	
	output_pile.close()
		
	if counter == 0:
		print("No SNPs passed your filters for: " + printable_name + ". There will be no output for this file.")
		os.remove(out)
		
	timer = datetime.now()
	printable_time = timer.strftime(time_format)
	print("SNP refinement for:", printable_name, "finished at:", printable_time)
	
	return None
	
def choose_consensus(directory, reference_file = ""):

	inputs = os.listdir(directory+"/MetaPop/06.SNP_Counts")	
	for i in range(0, len(inputs)):
		inputs[i] = directory+"/MetaPop/06.SNP_Counts/"+inputs[i]
	
	#We can just do the same code for the single input that works. Below will not add if contig + pos missing.
	if reference_file == "":
		print("Reference base at each position will be the consensus of all files.")
	else:
		print("Reference base at each position will be the consensus of:", reference_file)
		base_inputs = os.listdir(directory+"/MetaPop/06.SNP_Counts")	
		base_name = os.path.basename(os.path.normpath(reference_file))[:-4]+"_filtered_snps.txt"
		selected_input = inputs[base_inputs.index(base_name)]
		inputs = [selected_input]
			
	sequence_dict = {}
	
	for file in inputs:
		fh = open(file)
		
		for line in fh:
			segs = line.strip().split()
			contig = segs[0]
			#Position here is 1 based, so we have to correct.
			pos = int(segs[1])
			
			atcg_counts = list(map(int, segs[2:]))
			
			if contig not in sequence_dict:
				sequence_dict[contig] = {}
				sequence_dict[contig][pos] = atcg_counts
			else:
				if pos not in sequence_dict[contig]:
					sequence_dict[contig][pos] = atcg_counts
				else:
					sequence_dict[contig][pos][0] += atcg_counts[0]
					sequence_dict[contig][pos][1] += atcg_counts[1]
					sequence_dict[contig][pos][2] += atcg_counts[2]
					sequence_dict[contig][pos][3] += atcg_counts[3]
		
		fh.close()
	
	bases = ["A", "T", "C", "G"]
	
	consenSNP = open(directory+"/MetaPop/07.Cleaned_SNPs/consensus_bases.txt", "w")
			
	for contig in sequence_dict:
		for pos in sequence_dict[contig]:
			idx = sequence_dict[contig][pos].index(max(sequence_dict[contig][pos]))
			refbase = bases[idx]
			#output a file with the calls
			print(contig, pos, refbase, sequence_dict[contig][pos][0], sequence_dict[contig][pos][1], sequence_dict[contig][pos][2], sequence_dict[contig][pos][3], sep = "\t", file = consenSNP)
			#I reuse the index later, so it's better to leave it like this.
			sequence_dict[contig][pos] = idx
	
	consenSNP.close()
	
	#get all inputs again, in case a ref was given.
	inputs = os.listdir(directory+"/MetaPop/06.SNP_Counts")	
	for i in range(0, len(inputs)):
		inputs[i] = directory+"/MetaPop/06.SNP_Counts/"+inputs[i]
		
	return sequence_dict, inputs
	
def seq_to_codon_bias(seq):
	
	#Split into codons
	result = [seq[i:i+3] for i in range(0, len(seq), 3) ]
	
	codons = ["AAA", "AAT", "AAC", "AAG", 
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
			  "GGA", "GGT", "GGC", "GGG"]
			  
	#setup return vector
	counts = {}
	for key in codons:
		counts[key] = 0
	
	codons = set(codons)
	
	#Remove N-containing codons
	for amino in result:
		if amino in codons:
			counts[amino] += 1
			
	return counts.values()

def localize_to_codon(pos, start, stop, str, rbi):
	#doing the codon math here with mods.
	bases_f = ["A", "T", "C", "G"]
	bases_r = ["T", "A", "G", "C"]
	
	if str == "1":
		codon = ((pos - start)//3)+1
		pos_in_codon = ((pos - start) % 3)+1
		base = bases_f[rbi]
	else:
		codon = ((stop - pos )//3)+1
		pos_in_codon = ((end-pos) % 3) + 1
		base  = bases_r[rbi]
	
	
	pos_in_gene = (3 * (codon)-1) + pos_in_codon - 1
		
	return codon, pos_in_codon, base, pos_in_gene
	
def revcomp(bases):
	if len(bases) > 0:
		match = ["A", "T", "C", "G"]
		switch = ["T", "A", "G", "C"]
		
		spl = list(bases)		
		
		new = []
		
		for base in spl:
			new.append(switch[match.index(base)])
		
		bases = "".join(new)
		
	return bases
		
def associate_genes(directory, consensus, inputs, min_obs, min_pct):
	genes_file = directory + "/MetaPop/01.Genomes_and_Genes/all_genomes_genes.fasta"
	
	bases = ["A", "T", "C", "G"]
	
	#Collect the set of contig-pos combos found as SNP calls over all files
	acceptable_contigs = consensus.keys()
	acceptable_pos_by_contig = {}
	for key in acceptable_contigs:
		acceptable_pos_by_contig[key] = list(consensus[key].keys())
		
	contig = ""
	seq = ""
	contig_gene = ""
	
	gene_dict = {}
	
	#number -> decimal conversion
	real_min_pct = min_pct/100.0
	
	#We're gonna need start/stops for the genes to assign snps to them.
	print("Getting codon usage bias...")
		
	gene = open(genes_file)
		
	for line in gene:
		if line.startswith(">"):
		
			if len(seq) > 0 :
			
				if contig_gene != "MissingNo.":
					
					g = contig_gene.split("_")
					g = g[len(g)-1]
					
					c = contig_gene[1:len(contig_gene)-(len(g)+1)]
															
					if c in acceptable_contigs:
										
						if c in gene_dict:
							gene_dict[c][0].append(strand)
							gene_dict[c][1].append(start)
							gene_dict[c][2].append(end)
						else:
							gene_dict[c] = [[strand], [start], [end]]
				
				if len(line.strip()) > 1:
					segs = line.strip().split(" # ")
					contig_gene = segs[0]
					
					start = int(segs[1])
					end = int(segs[2])
					
					strand = bool(segs[3])
				else:
					print("Missing sequence ID found.")
					contig_gene = "MissingNo."
			
			else:
				#This is essentially just first line of the file
				if len(line.strip()) > 1:
					segs = line.strip().split(" # ")
					contig_gene = segs[0]
					start = int(segs[1])
					end = int(segs[2])

					strand = bool(segs[3])
					
					seq = "not_empty"
				else:
					print("Missing sequence ID found.")
					contig_gene = "MissingNo."
		
	#Last one doesn't normally happen
	if contig_gene != "MissingNo.":
		g = contig_gene.split("_")
		g = g[len(g)-1]
		c = contig_gene[1:len(contig_gene)-(len(g)+1)]
		
		if c in acceptable_contigs:
		
			if c in gene_dict:
				gene_dict[c][0].append(strand)
				gene_dict[c][1].append(start)
				gene_dict[c][2].append(end)
			else:
				gene_dict[c] = [[strand], [start], [end]]
		
	seq = ""
	
	gene.close()
	
	#Not every contig in consensus is guaranteed to be added to gene_dict, 
	acceptable_contigs_genes = gene_dict.keys()
	
	#Take care of links later
	#Links are separated by source, too.
	genic = open(directory + "/MetaPop/07.Cleaned_SNPs/genic_snps_uncorrected.tsv", "w")
	
	non_genic = open(directory + "/MetaPop/07.Cleaned_SNPs/non_genic_snps.tsv", "w")
	print("contig_pos\tcontig\tpos\tref_base\tdepth\ta_ct\tt_ct\tc_ct\tg_ct\tsource\tsnps", file = non_genic)
	
	#for genic snps, we need a dict of replacements for the gene updates later - add the gene, the pos. in gene, and the base for correction here.
	gene_correction = {}
	
	linked_snps_tracker = {}
	
	print("Finalizing SNPs...")
	
	for file in inputs:
		source_name = file.split("/")
		source_name = source_name[len(source_name)-1]
		source_name = source_name[0:(len(source_name)-18)]
			
		acceptable_bases = {}
		
		fh = open(file)
		
		#Count each time a codon appears per gene; gene : count_of_appearances.
		linked_snps_tracker[source_name] = {}
		
		for line in fh:
			segs = line.strip().split()
			contig = segs[0] 
			pos = int(segs[1])
			counts = list(map(int, segs[2:]))
			
			
			#Don't bother going further if the c-p combo isn't in the ref file; this shouldn't ever activate if no single ref file was selected.
			if contig not in acceptable_contigs:
				continue
			if pos not in acceptable_pos_by_contig[contig]:
				continue
			
			'''
			#Only print the SNPs that are in the set for each file
			if contig not in acceptable_bases:
				continue
			if pos not in acceptable_bases[contig]:
				continue
			'''
			
			
			refbaseidx = consensus[contig][pos]
			
			depth = sum(counts)
			
			#We don't want to call any reference bases as SNPs
			ref_ct = counts[refbaseidx]
			
			counts[refbaseidx] = 0
			
			snp_bases = []
			
			#This is always going to be the same, so just declare it rather than having the program do more work.
			for i in [0, 1, 2, 3]:
				if counts[i] > min_obs and counts[i]/depth > real_min_pct:
					snp_bases.append(bases[i])
			
			snps = "".join(snp_bases)
			
			#if len(snps) == 0:
				#continue
			
			#for print, reset the count of the ref base
			counts[refbaseidx] = ref_ct
			
			#If a contig has no genes predicted, the next lines will fail. This statement protects that.
			if contig not in acceptable_contigs_genes:
				continue
			
			#Starts should always be sorted
			#cannot guarantee ends are sorted.
			gene_starts = gene_dict[contig][1]
			gene_stops = gene_dict[contig][2]
			gene_strands = gene_dict[contig][0]
			
			#if the start is > pos, don't need to check
			#If the end is < pos, don't need to check
			
			#this finds the rightmost index of starts I have to search to by finding all of the ones that are smaller than pos
			low_idx = bisect.bisect_left(gene_starts, pos)
			
			was_genic = False
			
			for i in range(0, low_idx):
				if gene_stops[i] > pos:
				
					if contig+"_"+str(i+1) not in linked_snps_tracker[source_name]:
						linked_snps_tracker[source_name][contig+"_"+str(i+1)] = {}
				
					if strand:
						print_strand = "1"
						codon, pos_in_codon, base, pos_in_gene = localize_to_codon(pos, gene_starts[i], gene_stops[i], print_strand, refbaseidx)
						print_base = bases[refbaseidx]
							
					else:
						print_strand = "-1"
						counts = [counts[1], counts[0], counts[3], counts[2]]
						snps = revcomp(snps)
						print_base = revcomp(bases[refbaseidx])
						codon, pos_in_codon, base, pos_in_gene = localize_to_codon(pos, gene_starts[i], gene_stops[i], print_strand, refbaseidx)
						
					if contig+"_"+str(i+1) not in gene_correction:
						gene_correction[contig+"_"+str(i+1)] = [[pos_in_gene], [print_base]]
					else:
						#We don't want to add a base/gene pos pair twice.
						if pos_in_gene not in gene_correction[contig+"_"+str(i+1)][0]:
							gene_correction[contig+"_"+str(i+1)][0].append(pos_in_gene)
							gene_correction[contig+"_"+str(i+1)][1].append(print_base)
					
					
					if codon not in linked_snps_tracker[source_name][contig+"_"+str(i+1)]:
						linked_snps_tracker[source_name][contig+"_"+str(i+1)][codon] = False	
					else:
						linked_snps_tracker[source_name][contig+"_"+str(i+1)][codon] = True						
					
					
					#First pass, we check occurrences of 2 SNPs in one codon on one gene; second pass we correct false -> true where needed.
					#We have to go through the file again to check link = True
					
					if len(snps) > 0:
						print(contig+"_"+str(pos), contig, pos, bases[refbaseidx], depth, *counts, source_name, snps, contig+"_"+str(i+1), gene_starts[i], gene_stops[i], print_strand, codon, pos_in_codon, "False", file = genic, sep = "\t")
					
					was_genic = True
					
			if not was_genic:
				if len(snps) > 0:
					print(contig+"_"+str(pos), contig, pos, bases[refbaseidx], depth, *counts, source_name, snps, file = non_genic, sep = "\t")
				
			
		
		fh.close()
	
	genic.close()
	non_genic.close()
	
	genic_raw = open(directory + "/MetaPop/07.Cleaned_SNPs/genic_snps_uncorrected.tsv", "r")
	genic_links = open(directory + "/MetaPop/07.Cleaned_SNPs/genic_snps.tsv", "w")
	
	print("contig_pos\tcontig\tpos\tref_base\tdepth\ta_ct\tt_ct\tc_ct\tg_ct\tsource\tsnps\tcontig_gene\tstart\tend\tOC\tcodon\tpos_in_codon\tlink", file = genic_links)
	
	for line in genic_raw:
		segs = line.strip().split("\t")
		
		#print(segs)
		
		source = segs[9]
		gene = segs[11]
		codon = int(segs[15])
		
		#linked_snps_tracker[source][gene][codon] = False/True
		
		if linked_snps_tracker[source][gene][codon]:
			segs[17] = True
			print(*segs, sep = "\t", file = genic_links)
		else :
			print(line.strip(), file = genic_links)
	
	genic_raw.close()
	genic_links.close()
	
	#cleanup
	os.remove(directory + "/MetaPop/07.Cleaned_SNPs/genic_snps_uncorrected.tsv")
	
	#codon bias goes here, I guess?
	codon_seqs = open(directory + "/MetaPop/08.Codon_Bias/codon_bias_sequences.tsv", "w")
	codon_bias = open(directory + "/MetaPop/08.Codon_Bias/codon_counts_by_gene.tsv", "w")
	updated_genes = open(directory + "/MetaPop/01.Genomes_and_Genes/base_corrected_genes.fa", "w")
	
	gene = open(genes_file)
	
	#Codon bias stuff here - upon print, I need to still update the positions in each gene with the approp. bases; also write corrected genes here, might as well
	
	full_ID = ""
	
	acceptable_genes = gene_correction.keys()
	
	print("Updating genes with consensus bases...")
	
	for line in gene:
		if line.startswith(">"):
			
			if len(seq) > 0 :
			
				if contig_gene in acceptable_genes:
				
					if contig_gene != "MissingNo.":
						seq = list(seq)
						
						for i in range(0, len(gene_correction[contig_gene][0])):
							#The array's indices are 0-indexed, the gene pos are 1-indexed.
							pos_in_gene = gene_correction[contig_gene][0][i]-1
							base_replace = gene_correction[contig_gene][1][i]
							seq[pos_in_gene] = base_replace
						
						seq = "".join(seq)
						counts = seq_to_codon_bias(seq)


						
						print(seq, strand, contig_gene, parent_contig, sep = "\t", file = codon_seqs)

						print(*counts, sep = "\t", file = codon_bias)
						
						print(full_ID, file = updated_genes)
						#update seq
						print(seq, file = updated_genes)
				
				else:
					#just print the original if it doesn't need updates.
					print(seq, strand, contig_gene, parent_contig, sep = "\t", file = codon_seqs)

					counts = seq_to_codon_bias(seq)
					print(*counts, sep = "\t", file = codon_bias)
					print(full_ID, file = updated_genes)
					print(seq, file = updated_genes)
				
				if len(line.strip()) > 1:
					segs = line.strip().split(" # ")

					contig_gene = segs[0][1:]
					
					parent_contig = contig_gene.split("_")
					parent_contig = parent_contig[:(len(parent_contig)-1)]
					parent_contig = "_".join(parent_contig)

					strand = int(bool(segs[3]))

					full_ID = line.strip()
					
				else:
					print("Missing sequence ID found.")
					contig_gene = "MissingNo."
				
				seq = ""
			
			else:
				if len(line.strip()) > 1:
					segs = line.strip().split(" # ")
					contig_gene = segs[0][1:]

					parent_contig = contig_gene.split("_")
					parent_contig = parent_contig[:(len(parent_contig)-1)]
					parent_contig = "_".join(parent_contig)

					strand = int(bool(segs[3]))

					full_ID = line.strip()
				else:
					print("Missing sequence ID found.")
					contig_gene = "MissingNo."
		else:
			seq += line.strip()
	
	#Last one doesn't normally happen
	if contig_gene != "MissingNo.":
	
		if contig_gene in acceptable_genes:
			seq = seq.split("")
					
			for i in range(0, len(gene_correction[contig_gene][0])):
				#The array's indices are 0-indexed, the gene pos are 1-indexed.
				pos_in_gene = gene_correction[contig_gene][0][i]-1
				base_replace = gene_correction[contig_gene][1][i]
				seq[pos_in_gene] = base_replace
						
			seq = "".join(seq)
			counts = seq_to_codon_bias(seq)
			print(seq, strand, contig_gene, parent_contig, sep = "\t", file = codon_seqs)
			print(*counts, sep = "\t", file = codon_bias)
			print(full_ID, file = updated_genes)
			print(seq, file = updated_genes)
		
		else:
			print(seq, strand, contig_gene, parent_contig, sep = "\t", file = codon_seqs)
			counts = seq_to_codon_bias(seq)
			print(*counts, sep = "\t", file = codon_bias)
			print(full_ID, file = updated_genes)
			print(seq, file = updated_genes)
	
	seq = ""
	
	gene.close()
	codon_seqs.close()
	codon_bias.close()
	updated_genes.close()
	
	OG_genomes = open(directory + "/MetaPop/01.Genomes_and_Genes/all_genomes.fasta", "r")
	updated_genomes = open(directory + "/MetaPop/01.Genomes_and_Genes/base_corrected_genomes.fa", "w")
	
	acceptable_contigs = consensus.keys()
	
	print("Updating genomes with consensus bases...")
	
	for line in OG_genomes:
		if line.startswith(">"):
			
			if len(seq) > 0 :
			
				if contig in acceptable_contigs:
				
					if contig != "MissingNo.":
						seq = list(seq)
												
						for i in consensus[contig]:
							#i is a position in the genome, but the same 0 and 1 -indexing applies here
							seq[i-1] = bases[consensus[contig][i]]
						
						seq = "".join(seq)				
						
						print(">"+contig, file = updated_genomes)
						#update seq
						print(seq, file = updated_genomes)
				
				else:						
					print(">"+contig, file = updated_genomes)
					#update seq
					print(seq, file = updated_genomes)
				
				if len(line.strip()) > 1:
					segs = line.strip().split()
					contig = segs[0][1:]
					
					
				else:
					print("Missing sequence ID found.")
					contig = "MissingNo."
				
				seq = ""
			
			else:
				if len(line.strip()) > 1:
					segs = line.strip().split()
					contig = segs[0][1:]
				else:
					print("Missing sequence ID found.")
					contig = "MissingNo."
		else:
			seq += line.strip()
	
	#Last one doesn't normally happen
	if contig != "MissingNo.":
	
		if contig in acceptable_genes:
			seq = seq.split("")
					
			for i in consensus[contig]:
				seq[i-1] = bases[consensus[contig][i]]
						
			seq = "".join(seq)
			
			print(">"+contig, file = updated_genomes)
			#update seq
			print(seq, file = updated_genomes)
			
		else:
			print(">"+contig, file = updated_genomes)
			#update seq
			print(seq, file = updated_genomes)
	
	seq = ""
	
	updated_genomes.close()
	OG_genomes.close()
		
	#R scripts are the only things that are going to use these, so dirs can be relative.
	return "MetaPop/01.Genomes_and_Genes/base_corrected_genomes.fa", "MetaPop/01.Genomes_and_Genes/base_corrected_genes.fa"
	


