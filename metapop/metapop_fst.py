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


def contig_lengths(file):
	lengths = {}
	header = True
	fh = open(file)
	for line in fh:
		if header:
			header = False
		else:
			segs = line.strip().split('\t')
			contig = segs[0]
			contig = contig.split()[0]
			bp = int(segs[4])
			
			lengths[contig] = bp
	
	fh.close()
	
	return lengths

def format_microdiv(file, lengths):
	microdiv = {}
	contig_data = {}
	
	all_sources = []
	
	fh = open(file)	
	header = True
	for line in fh:
		segs = line.strip().split('\t')
		if header:
			header = False
			colnames = segs
			#print(colnames)
			continue
			
		contig = segs[1]
		contig = contig.split()[0]
		source = segs[9]
		
		all_sources.append(source)
		
		#We don't need this to be an int.
		position = segs[2]
		pi = float(segs[23])
		a, t, c, g = int(segs[18]), int(segs[19]), int(segs[20]), int(segs[21]),
		sub_samp_depth = int(segs[22])
		#Needs pos, cov, pi, and subsamp A, T, C, G
		
		
		if contig not in microdiv:
			microdiv[contig] = {}
			#Add length
		if contig not in contig_data:
			#Genome length
			contig_data[contig] = {}
			contig_data[contig]["length"] = lengths[contig]
			
		if source not in microdiv[contig]:
			#List of the data.
			microdiv[contig][source] = {}
			
		microdiv[contig][source][position] = [a, t, c, g, pi, sub_samp_depth]
		
	fh.close()
	
	all_sources = list(set(all_sources))
	
	snp_cts = {}
	to_remove = []
	
	for ct in microdiv:
		if len(microdiv[ct]) < 2:
			to_remove.append(ct)
		else:
			#We don't need to do this work unless the contig is going to actually be used.
			snps = []
			for source in microdiv[ct]:
				#The snps seen on this contig.
				my_snps = list(microdiv[ct][source].keys())
				#Build the set of total SNPs for all contigs
				snps.extend(my_snps)
			#This ends up being the length of the total set of SNPs across all samples
			#unique positions and add total SNP count
			contig_data[ct]["total_snp_count"] = len(list(set(snps)))
	
	#Clean up the data	
	for ct in to_remove:
		#Ensure that the data is removed nicely
		microdiv[ct] = None
		#Get rid of the row.
		microdiv.pop(ct)
	
	return microdiv, contig_data, all_sources
	

#Return all pairwise FST values for a contig
def fst_pair(one_contig_args):
	data, contig_dat, contig_name = one_contig_args[0], one_contig_args[1], one_contig_args[2]
	
	#total_snps = contig["total_snp_count"]
	#contig_length = contig["length"]
	
	fsts = {}
	sources = list(data.keys())
	for i in range(0, len(sources)-1):
		#Select query sample
		query = sources[i]
		#Select query data from the contig's list
		query_data = data[query]
		#Results holder
		fsts[query] = {}
		#self is always 0, so we don't need it
		#fsts[query][query] = 0
		
		cov1 = contig_dat["length"] - (contig_dat["total_snp_count"] - len(query_data))
		pi1 = 0
		for p in query_data:
			try:
				pi1 += query_data[p][4]
			except:
				pass
			
		pi1 = pi1 / cov1
		
		for j in range(i+1, len(sources)):
			target = sources[j]
			
			#Default value of FST is 1, set it here in case there are no shared SNP positions
			fsts[query][target] = 1
			#Select target data from contig's list
			target_data = data[target]
			
			cov2 = contig_dat["length"] - (contig_dat["total_snp_count"] - len(target_data))
			pi2 = 0
			for p in target_data:
				try:
					pi2 += target_data[p][4]
				except:
					pass
			
			pi2 = pi2 / cov2
			
			#Select shared SNP positions
			shared_positions = list(set(list(query_data.keys())).intersection(list(target_data.keys())))

			
			#FST is calculated over shared positions
			if len(shared_positions) > 0:
				#find the adjusted denominator
				
				shared_cov = contig_dat["length"] - (contig_dat["total_snp_count"] - len(shared_positions))
				#print("shared cov is", shared_cov)
				
				per_pos_fst = 0
				
				for pos in shared_positions:
					#arrays are a, t, c, g subsamples, pi, and sub_samp_depth
					q = query_data[pos]
					t = target_data[pos]
					
					#Subsample depth for query and target multiplied
					cc = q[5]*t[5]
					#Pairwise pi values
					at = (q[0] * t[1])/cc
					ac = (q[0] * t[2])/cc
					ag = (q[0] * t[3])/cc
					ta = (q[1] * t[0])/cc
					tc = (q[1] * t[2])/cc
					tg = (q[1] * t[3])/cc
					ca = (q[2] * t[0])/cc
					ct = (q[2] * t[1])/cc
					cg = (q[2] * t[3])/cc
					ga = (q[3] * t[0])/cc
					gt = (q[3] * t[1])/cc
					gc = (q[3] * t[2])/cc
					
					this_pos_fst = at+ac+ag+ta+tc+tg+ca+ct+cg+ga+gt+gc
					
					per_pos_fst += this_pos_fst
					
				try:
					
					fst = per_pos_fst/shared_cov
					
					fst = 1-(( (pi1+pi2)/2) / fst)
					
					if fst < 0:
						fst = 0
					
				except:
					#Divide by zero happened because SNP positions weren't shared for any pairing
					fst = 1
					
				fsts[query][target] = fst
						
					
	#The final iteration is the last item in source and is a self vs. self, which also always has FST = 0
	#fsts[target] = {target: 0}
	
	return fsts, contig_name
			
def calculate_FST(microdiv, contig_data, all_sources, output, threads):
	success = True
	
	args = []
		
	for ct in microdiv:
		args.append((microdiv[ct], contig_data[ct], ct,))
			
	
	fh = open(output, "w")
	print("row_samp", "col_samp", "contig", "fst", sep = "\t", file = fh)
	pool = multiprocessing.Pool(threads)
	
	#We also need to print the "missing" pairs, so we use the list of sources
	for result in pool.map(fst_pair, args):
		dat = result[0]
		ct_name = result[1]
		for i in range(0, len(all_sources) - 1):
			query = all_sources[i]
			for j in range(i+1, len(all_sources)):
				target = all_sources[j]
				
				if query in dat:
					if target in dat[query]:
						print(query, target, ct_name, dat[query][target], sep = "\t", file = fh)
					else:
						print(query, target, ct_name, "NA", sep = "\t", file = fh)
				else:
					print(query, target, ct_name, "NA", sep = "\t", file = fh)

		#for query in dat:
		#	for target in dat[query]:
		#		print(query, target, ct_name, dat[query][target], sep = "\t", file = fh)
	
	pool.close()
	pool.join()
	
	fh.close()
		
	return success
	
	
#Main function
def perform_fst(microdiv_file, lengths_file, out_dir, threads):
	output_file = os.path.normpath(out_dir + "/MetaPop/10.Microdiversity/fixation_index.tsv")
	
	#Get length of each contig
	ct_lens = contig_lengths(lengths_file)

	#Read in and format data; also get total SNP counts for valid genome+sample combos
	md, ctd, sources = format_microdiv(microdiv_file, ct_lens)
	
	#for ct in md:
	#	print(ct, md[ct])
	
	#We don't need this anymore
	ct_lens = None
	
	#Calculate FST and write to output.
	success = calculate_FST(md, ctd, sources, output_file, threads)
	
	if not success:
		print("FST not successfully calculated.")
	
	return success
	
#Testing.	
#lens = "/mnt/c/Users/Kenji/Desktop/metapy/toy_data/MetaPop/10.Microdiversity/global_contig_microdiversity.tsv"
#file = "/mnt/c/Users/Kenji/Desktop/metapy/toy_data/MetaPop/10.Microdiversity/global_raw_microdiversity_data_snp_loci_only.tsv"

#perform_fst(file, lens, "test_toy", 6)
