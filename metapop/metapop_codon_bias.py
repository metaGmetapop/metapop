#This script works. Keep as a backup.
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
import numpy as np
from collections import defaultdict


def codon_converter():
	cols = ["AAA", "AAT", "AAC", "AAG",
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
			
	values = []
	for i in range(0, 64):
		values.append(i)
		
	cd_dict = dict(zip(cols, values))
	
	return cd_dict

def aa_converter():
	
	aa_set = {}
	
	aa_set["I"]  = ["ATT", "ATC", "ATA"]
	aa_set["L"]  = ["CTT", "CTC", "CTA", "CTG", "TTA", "TTG"]
	aa_set["V"]  = ["GTT", "GTC", "GTA", "GTG"]
	aa_set["FE"] = ["TTT", "TTC"]
	aa_set["C"]  = ["TGT", "TGC"]
	aa_set["A"]  = ["GCT", "GCC", "GCA", "GCG"]
	aa_set["G"]  = ["GGT", "GGC", "GGA", "GGG"]
	aa_set["P"]  = ["CCT", "CCC", "CCA", "CCG"]
	aa_set["TH"] = ["ACT", "ACC", "ACA", "ACG"]
	aa_set["S"]  = ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"]
	aa_set["Y"]  = ["TAT", "TAC"]
	aa_set["Q"]  = ["CAA", "CAG"]
	aa_set["N"]  = ["AAT", "AAC"]
	aa_set["H"]  = ["CAT", "CAC"]
	aa_set["E"]  = ["GAA", "GAG"]
	aa_set["D"]  = ["GAT", "GAC"]
	aa_set["K"]  = ["AAA", "AAG"]
	aa_set["R"]  = ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"]
	aa_set["ST"] = ["TAA", "TAG", "TGA"]
	
	codons = codon_converter()
	
	aa_to_idx = {}
	for aa in aa_set:
		if aa not in aa_to_idx:
			aa_to_idx[aa] = []
			
		for cd in aa_set[aa]:
			aa_to_idx[aa].append(codons[cd])
			
	return aa_set, aa_to_idx
	
	
	
def calculate_codon_bias_per_genome(gene_file, counts_file):
	'''
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
	'''
	
	genes = []
	genomes = []
	fh = open(gene_file)
	for line in fh:
		segs = line.strip().split('\t')
		gene = segs[2].split()[0]
		genes.append(gene)
		genome = segs[3].split()[0]
		genomes.append(genome)
		
	fh.close()
	
	counts = np.loadtxt(counts_file, dtype = np.int32, delimiter = "\t")
	
	converter, col_index = aa_converter()
	
	#There's only 19 in here because we're including the stop codon and don't need to worry about trptophan or start.
	#This is transposed on purpose right now.
	#aa_counts = np.zeros(shape = (19, counts.shape[0]), dtype = np.int32)
	cur_col = 0
	
	proportions = np.zeros(shape = counts.shape, dtype = np.float32)
	
	#don't yell at me over div by zero.
	np.seterr(divide='ignore')
	
	for aa in col_index:
		colsum = np.sum(counts[:, col_index[aa]], axis = 1)
		#aa_counts[cur_col] = colsum
		proportions[:, col_index[aa]] = np.divide(counts[:, col_index[aa]], colsum[:, None])
		cur_col += 1
		#counts[:, col_index[aa]]
		
	#nans to zeros
	np.nan_to_num(proportions, copy = False)
	#we can ignore tryp. and start because distance of 0.0 = 0
	
	#perform by-genome grouped selections
	split_indices = []
	count = 0
	last_genome = genomes[0]
	for g in genomes:
		if g == last_genome:
			count += 1
		else:
			split_indices.append(count)
			#set new genome
			last_genome = g
			count += 1
			
	#divide into per-genome sets
	per_genome_counts = np.vsplit(counts, split_indices)
	#per_genome_aa_sums = np.vsplit(aa_counts, split_indices)
	per_genome_codon_props = np.vsplit(proportions, split_indices)
	
	averages, iq_ranges = [], []
	
	#euc_dist_out = open("MetaPop/08.Codon_Bias/gene_euclidean_distances.tsv", "w")
	euc_dist_out = open("gene_euclidean_distances.tsv", "w")
	
	for i in range(0, len(per_genome_codon_props)):
	#for i in range(0, 1):
	
		cur_count_set = per_genome_counts[i]
		
		print(cur_count_set.shape)
		
		#This is the per-gene CB usage
		cur_prop_set = per_genome_codon_props[i]
		
		mean_cts = np.sum(cur_count_set, axis = 0)
		mean_props = np.zeros(64, dtype = np.float32)
		for aa in col_index:
			colsum = np.sum(mean_cts[col_index[aa]])
			#aa_counts[cur_col] = colsum
			mean_props[col_index[aa]] = np.divide(mean_cts[col_index[aa]], colsum)
		
		#fix nans if present
		np.nan_to_num(mean_props, copy = False)
		
		#euc. dist from the average gene proportions
		dists = np.square(np.subtract(cur_prop_set, mean_props))
		
		#sum dists per gene, sqrt
		dists = np.sqrt(np.sum(dists, axis = 1))
		
		iqr = np.quantile(dists, [0.25, 0.75])
		iqr = iqr[1] - iqr[0]
		mu = np.mean(dists)
		
		averages.append(mu)
		iq_ranges.append(iqr)
		
		is_outlier = dists > (mu + (1.5 * iqr))
		
		#So we can just print out the stuff here.
		for j in range(0, len(is_outlier)):
			print(genomes.pop(0), genes.pop(0), dists[j], is_outlier[j], sep = "\t", file = euc_dist_out)
		
	euc_dist_out.close()
	
		
		
		
		
	

		
	
#genes <- fread("MetaPop/08.Codon_Bias/codon_bias_sequences.tsv", sep = "\t", header = F, select = c(3,4))
#colnames(genes) = c("V1", "V2")
	
genes = sys.argv[1]
counts = sys.argv[2]


calculate_codon_bias_per_genome(genes, counts)

	
	
	
	
	
	
	
	
	
	
	
	
	