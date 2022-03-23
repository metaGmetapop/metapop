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

def read_and_format_linked_candidates(file):
	linked_data = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list)))))
	
	fh = open(file)
	for line in fh:
		if line.endswith("True\n"):
			segs = line.strip().split("\t")
			#contig_pos = segs[0]
			contig = segs[1]
			contig = contig.split()[0]
			pos = int(segs[2])
			ref_base = segs[3]
			source = segs[9]
			snps = segs[10]
			contig_gene = segs[11]
			contig_gene = contig_gene.split()[0]
			#if OC == 1, strand = forward, else strand = reverse
			OC = int(segs[14])
			codon = int(segs[15])
			pos_in_codon = int(segs[16])
			
			linked_data[source][contig][contig_gene][codon][OC].append([pos, ref_base, snps, pos_in_codon])
				
	fh.close()
	
	#Convert defdicts into regular dicts. If a leaf node of a dict is empty, trim upwards until the leaf is no longer empty.
	selections = {}
	
	#We need to build a copy because altering dict size in loops is a no-no.
	cleaned_data = dict()

	linked_data = dict(linked_data)
	for source in linked_data:
		cleaned_data[source] = dict(linked_data[source])
		for contig in linked_data[source]:
			cleaned_data[source][contig] = dict(linked_data[source][contig])
			for gene in linked_data[source][contig]:
				cleaned_data[source][contig][gene] = dict(linked_data[source][contig][gene])
				for codon in linked_data[source][contig][gene]:
					cleaned_data[source][contig][gene][codon] = dict(linked_data[source][contig][gene][codon])
					for strand in linked_data[source][contig][gene][codon]:
						#Only select those items for which a true match was found.
						if len(linked_data[source][contig][gene][codon][strand]) > 1:
							#Format the request here
							set = []
							positions_in_codon = []
							for l in linked_data[source][contig][gene][codon][strand]:
								set.append(l[0])
								positions_in_codon.append(l[3])
							set.sort()
							if source in selections:
								selections[source].append((contig, set, positions_in_codon))
								
							else:
								selections[source] = [(contig, set, positions_in_codon)]								
						else:
							cleaned_data[source][contig][gene][codon].pop(strand)
					if len(cleaned_data[source][contig][gene][codon]) == 0:
						cleaned_data[source][contig][gene].pop(codon)
				if len(cleaned_data[source][contig][gene]) == 0:
					cleaned_data[source][contig].pop(gene)
			if len(cleaned_data[source][contig]) == 0:
				cleaned_data[source].pop(contig)
		if len(cleaned_data[source]) == 0:
			cleaned_data.pop(source)
				
	#Clean up.
	linked_data = None
		
	return cleaned_data, selections

def parse_cig_to_pos(cig, seq, left_side, positions):
	
	#We only accept codons that are found in consecutive matches covering the codon. 
	#If a codon we're looking for isn't all matches, we disregard it.
	ret = None
	pos_in_read = 0
	pos_in_genome = 0
	low = min(positions) - left_side
	high = max(positions) - left_side
	
	for ct_op in cig:
		count, op = int(ct_op[0]), ct_op[1]
		if op == "M" or op == "X" or op == "=":
			if pos_in_genome <= low and pos_in_genome + count >= high:
				#How much further we have to go...
				pad_left = low - pos_in_genome
				pad_right = high - pos_in_genome
				ret = [seq[pad_left + pos_in_read], seq[pad_right+ pos_in_read]]
				#print(cig, positions, ret, pad_left, pad_right, pos_in_read, )
				#And we're done, stop looking.
				break
			#Advance spot in both read and genome - a match does both
			pos_in_genome += count
			pos_in_read += count
			continue
			
		if op == "I" or op == "S":
			#Advance spot in read - insertions don't move along the genome
			pos_in_read += count
			continue
		if op == "D" or op == "N":
			#Advance spot in genome - deletions aren't present in the read, but are present in the genome.
			pos_in_genome += count
			continue
		#H and P ops also technically exist, but we don't want to do anything with them.
		if op == "H" or op == "P":
			continue
		#Something was wrong in the op, we just skip it.
		break
			
	return ret
			
def read_one_range(command):
	gss = command[4]
	command = command[0:4]
	
	reader = subprocess.Popen(command, stdout = subprocess.PIPE)
	
	results = {}
	
	for read in reader.stdout:
		segs = read.strip().decode().split("\t")
		
		genome = segs[2]
		genome = genome.split()[0]
		
		#Leftmost is 1-based.
		leftmost = int(segs[3])
		
		cig = segs[5]
		seq = segs[9]
		
		parsed_cig = []
		for c in re.findall(r"([0-9]+)([a-z]+)", cig, re.I):
			parsed_cig.append(c)
		
		if len(parsed_cig) == 0:
			continue
		
		if len(parsed_cig) == 1:
			#This is all matches; selection of the codon is simple.
			chars = []
			for i in range(0, len(gss)):
				try:
					chars.append(seq[gss[i] - leftmost])
				except:
					pass
					
			chars = ''.join(chars)
			#Skip meaningless.
			if len(chars) == 0:
				continue
				
			if chars not in results:
				results[chars] = 1
			else:
				results[chars] += 1
		else:
			#There's either indels or something else complicating the selection of the codon. This function parses the cigar string to select the correct bases.
			try:
				chars = parse_cig_to_pos(parsed_cig, seq, leftmost, gss)
			except:
				chars = None
			
			if chars is not None:
				#If chars is none, there was a parse fail in the cigar string. We just skip those.
				chars = ''.join(chars)
				if chars not in results:
					results[chars] = 1
				else:
					results[chars] += 1
	
	keys = list(results.keys())
	#Just in case, clean up
	for k in keys:
		if len(k) < 2:
			results.pop(k)
	
	return results
	
def access_read_ranges(selected, threads, directory_base):
	ranges = []
	post_data = []
	for source in selected:
		file = os.path.normpath(directory_base + "/MetaPop/02.Filtered_Samples/" + source + "_preprocessed.bam")
		for s in selected[source]:
			genome = s[0]
			begin, end = min(s[1]), max(s[1])
			ranges.append(("samtools", "view", file, genome + ":" + str(begin)+"-"+str(end), s[1]))
			post_data.append([source, genome, s[1], s[2]])
	
	results = {}
	pool = multiprocessing.Pool(threads)
	res = pool.map(read_one_range, ranges)
	pool.close()
	pool.join()
	
	total = []
	
	for i in range(0, len(post_data)):
		total.append([*post_data[i], res[i]])
	
	return total

#Most of these should never activate.
def clean_format_codons(bases, positions):
	#We do revcomp later, so these are just put together in 1, 2, 3 order currently, regardless of their being forward or reverse strand.
	if positions[0] == 1:
		if positions[1] == 2:
			#1, 2
			bases = ''.join([bases[0], bases[1], '*'])
		if positions[1] == 3:
			#1, 3
			bases = ''.join([bases[0], '*', bases[1]])
	if positions[0] == 2:
		if positions[1] == 1:
			#2, 1 gets reversed
			bases = ''.join([bases[1], bases[0], '*'])
		if positions[1] == 3:
			#2, 3
			bases = ''.join(['*', bases[0], bases[1]])
	if positions[0] == 3:
		if positions[1] == 2:
			#3, 2 gets reversed
			bases = ''.join(['*', bases[1], bases[0]])
		if positions[1] == 1:
			#3, 1 gets reversed
			bases = ''.join([bases[1],'*', bases[0]])
	return bases
	

def process_snp_set(combo, strand):
	#how they appear in data
	as_in_data = []
	real_codons = []
	refs = []
	alts = []
	positions = []
	
	if strand == -1:
		combo.reverse()
		
	for one_pos in combo:
		refs.append(one_pos[1])
		alts.append(list(one_pos[2]))
		positions.append(one_pos[3])
	
	all_refs = []
	all_snps = []
	ref_first = []
	ref_second = []
	if len(positions) == 3:
		#This never changes, so it's literally out of the loop
		three_refs = ''.join(refs)
		for x in alts[0]:
			for y in alts[1]:
				for z in alts[2]:

					#Triplet refs are the only way a triple SNP is reported.
					all_refs.append(three_refs)
					all_refs.append(three_refs)
					all_refs.append(three_refs)
					all_refs.append(three_refs)
					#Doublets of refs
					#all_refs.append(''.join([refs[0], refs[1], '*']))
					#all_refs.append(''.join(['*', refs[1], refs[2]]))
					#all_refs.append(''.join([refs[0], '*', refs[2]]))
					
					#Triplet SNPs
					all_snps.append(''.join([x, y, z]))
					#Doublets of SNPs
					all_snps.append(''.join([x, y, '*']))
					all_snps.append(''.join(['*', y, z]))
					all_snps.append(''.join([x, '*', z]))
					
					#A placeholder for matching any triplet
					ref_first.append("***")
					#reference in first pos, 2 SNP combos of (1, 2), (2, 3), (1, 3)
					ref_first.append(''.join([refs[0], y, '*']))
					ref_first.append(''.join(['*', refs[1], z]))
					ref_first.append(''.join([refs[0], '*', z]))
					
					#A placeholder for matching any triplet
					ref_second.append("***")
					#reference in last pos, 2 SNP combos of (1, 2), (2, 3), (1, 3) 
					ref_second.append(''.join([x, refs[1], '*']))
					ref_second.append(''.join(['*', y, refs[2]]))
					ref_second.append(''.join([x, '*', refs[2]]))
		
	else:
		two_refs = clean_format_codons(refs, positions)
		for x in alts[0]:
			for y in alts[1]:
				#Each pairing has a unique set of SNPs it looks at, but the ref is the same
				all_refs.append(two_refs)
				#Both SNPs
				all_snps.append(clean_format_codons([x, y], positions))
				ref_first.append(clean_format_codons([refs[0], y], positions))
				ref_second.append(clean_format_codons([x, refs[1]], positions))
				
	return all_refs, all_snps, ref_first, ref_second
	
	
def formatted_to_valid_combos(snps, mined_reads, output):
	out = open(output, "w")
	#Header
	print("data_ref", "source", "contig_gene", "OC", "codon", "refs", "snp", "ref_count", "snp_count", "ref_first", "ref_second", sep = "\t", file = out)

	current_results = 0
	refiltered = []
	for source in snps:
		for genome in snps[source]:
			for gene in snps[source][genome]:
				for codon in snps[source][genome][gene]:
					for strand in snps[source][genome][gene][codon]:
						one_row = snps[source][genome][gene][codon][strand]
						mined_data = mined_reads[current_results]
						current_results += 1
						all_refs, all_snps, ref_snp, snp_ref = process_snp_set(one_row, strand)
						#print(all_refs, all_snps, ref_snp, snp_ref)
						
						mined_dict = mined_data[4]
						key_results = list(mined_dict.keys())
						
						match_positions = mined_data[2]
						codon_positions = mined_data[3]
						codon_positions.sort()
						
						new_keys = []
						for key in key_results:
							
							if 1 < len(key) < 3:
								#These are the only ones possible
								if codon_positions == [1, 2]:
									new = key + "*"
								if codon_positions == [2, 3]:
									new = "*" + key
								if codon_positions == [1, 3]:
									new = key[0] + "*" + key[1]
								new_keys.append(new)
							else:
								new_keys.append(key)
								
								
						matchable_dict = {}
						for i in range(0, len(key_results)):
							matchable_dict[new_keys[i]] = mined_dict[key_results[i]]
						
						#split_keys = [list(k) for k in new_keys]
						
						#There should always be the same length of each array, here
						for i in range(0, len(all_refs)):
							
							#2 is pretty easy. Just check.
							if len(codon_positions) == 2:
								
								ref = all_refs[i]
								snp = all_snps[i]
								rs  = ref_snp[i]
								sr  = snp_ref[i]
								
								if ref in matchable_dict:
									ref_ct = matchable_dict[ref]
								else:
									ref_ct = 0
								if snp in matchable_dict:
									snp_ct = matchable_dict[snp]
								else:
									snp_ct = 0
								if rs in matchable_dict:
									ref_first = matchable_dict[rs]
								else:
									ref_first = 0
								if sr in matchable_dict:
									ref_second = matchable_dict[sr]
								else:
									ref_second = 0
							
							#3 pos is so much harder
							else:
								ref = all_refs[i].replace("*", ".")
								snp = all_snps[i].replace("*", ".")
								rs  = ref_snp[i].replace("*", ".")
								sr  = snp_ref[i].replace("*", ".")
							
								ref_regex = re.compile(ref)
								snp_regex = re.compile(snp)
								rs_regex = re.compile(rs)
								sr_regex = re.compile(sr)
								
								ref_hits = list(filter(ref_regex.match, new_keys))
								snp_hits = list(filter(snp_regex.match, new_keys))
								rs_hits = list(filter(rs_regex.match, new_keys))
								sr_hits = list(filter(sr_regex.match, new_keys))
								
								#Clean out the parental hits
								for hit in ref_hits:
									if hit in rs_hits:
										rs_hits.pop(rs_hits.index(hit))
									if hit in sr_hits:
										sr_hits.pop(sr_hits.index(hit))
										
								for hit in snp_hits:
									if hit in rs_hits:
										rs_hits.pop(rs_hits.index(hit))
									if hit in sr_hits:
										sr_hits.pop(sr_hits.index(hit))
								
								ref_ct = 0
								snp_ct = 0
								
								ref_first = 0
								ref_second = 0
								
								for hit in ref_hits:
									ref_ct += matchable_dict[hit]
								for hit in snp_hits:
									snp_ct += matchable_dict[hit]
								for hit in rs_hits:
									ref_first += matchable_dict[hit]
								for hit in sr_hits:
									ref_second += matchable_dict[hit]
								
								#Triplet case is further special; ref first and second will be the same, and we int divide the results
								#and put half into each diagonal (worst case for showing linked SNPs)
								if rs == "...":
									ref_second = ref_first // 2
									ref_first = ref_first - ref_second
									
							#phi = calc_phi(ref_ct, snp_ct, ref_first, ref_second)
							print(current_results, source, gene, strand, codon, ref, snp, ref_ct, snp_ct, ref_first, ref_second, sep = "\t", file = out)
						

	out.close()
	
#Main function
def do_mine_reads(output_directory, threads):
	snps = output_directory + "/MetaPop/07.Cleaned_SNPs/genic_snps.tsv"
	output_file = output_directory + "/MetaPop/09.Linked_SNPs/linked_snp_results.tsv"
	
	formatted_snps, selections_to_read = read_and_format_linked_candidates(snps)
	res = access_read_ranges(selections_to_read, threads, output_directory)
	
	formatted_to_valid_combos(formatted_snps, res, output_file)
	
	return output_file


	















