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


def filter_thread_starter(_cmd, _cld, _im, _imd):
	global contig_mag_dict
	global contig_length_dict
	global is_mags
	global inverted_mag_dict
	contig_mag_dict = _cmd
	contig_length_dict = _cld
	is_mags = _im
	inverted_mag_dict = _imd

#The main function for this section. Ensures the proper directories exist and issues calls to the child functions; manages parallelism per sample
def filt(files, command_base, b_d_base, threads, mag_dict, length_dict, ref_file):

	if not os.path.exists(command_base[1] + "/MetaPop/02.Filtered_Samples"):
		os.mkdir(command_base[1] + "/MetaPop/02.Filtered_Samples")
	if not os.path.exists(command_base[1] + "/MetaPop/03.Breadth_and_Depth"):
		os.mkdir(command_base[1] + "/MetaPop/03.Breadth_and_Depth")
	if not os.path.exists(command_base[1] + "/MetaPop/04.Depth_per_Pos"):
		os.mkdir(command_base[1] + "/MetaPop/04.Depth_per_Pos")

		
	time_format = "%d/%m/%Y %H:%M:%S"	

	bams = os.listdir(files)
	
	cleaned_bams = []
	
	for file in bams:
		if file.upper().endswith(".BAM") or file.upper().endswith(".SAM"):
			cleaned_bams.append(file)
	
	sorts = []
	index = []
	commands = []
	b_and_d = []
	
	md_files = []
	
	checked_bams = []
	
	for sample in cleaned_bams:
		base_name = os.path.basename(os.path.normpath(sample))
		if base_name.endswith(".bam") or base_name.endswith(".BAM"):
			base_name = base_name[:-4]
			
		timer = datetime.now()
		printable_time = timer.strftime(time_format)
		print("Checking:", base_name, "for correctness... ", end = "", flush = True)		

		sanity = ["samtools", "quickcheck", os.path.normpath(files+"/"+sample)]
		sane = subprocess.Popen(sanity, stderr = subprocess.PIPE)
		counter = 0
		for line in sane.stderr:
			counter += 1
		#go to the next file if the file was incorrectly fmted.
		if counter > 0:
			print("the file did not appear to be a SAM or BAM format file. Skipping this file.")
			continue
		
		print("file passed.")
		checked_bams.append(sample)
		
	for sample in checked_bams:
		base_name = os.path.basename(os.path.normpath(sample))
		if base_name.endswith(".bam") or base_name.endswith(".BAM"):
			base_name = base_name[:-4]
			
		timer = datetime.now()
		printable_time = timer.strftime(time_format)
		print("Pre-formatting:", base_name, "started at:", printable_time + "...", end = '', flush = True)
	
		#Sort also converts SAM -> BAM and here makes a copy of the original that I am OK to modify or destroy w/o touching originals.
		sort = ["samtools", "sort", "-@", str(threads), "-o", command_base[1] + "/MetaPop/02.Filtered_Samples/"+base_name+"_initial.bam", files+"/"+sample]
		subprocess.call(sort, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		index = ["samtools", "index", "-@", str(threads), command_base[1] + "/MetaPop/02.Filtered_Samples/"+base_name+"_initial.bam"]
		#Index should naturally produce no stdout or stderr.
		subprocess.call(index)
		
		timer = datetime.now()
		printable_time = timer.strftime(time_format)
		print(" finished at", printable_time)
		
		md_files.append([command_base[1] + "/MetaPop/02.Filtered_Samples/"+base_name+"_initial.bam", base_name, ref_file, command_base[1]])
	
	pool = multiprocessing.Pool(min(threads, len(md_files)))
	md_filled = pool.map(fill_mdz, md_files)
	pool.close()
	
	for sample in md_filled:
		#base name based on mdz file, now
		base_name = os.path.basename(os.path.normpath(sample))
		base_name = base_name[:-14]
		
		timer = datetime.now()
		printable_time = timer.strftime(time_format)
		print("Sorting:", base_name, "started at:", printable_time + "...", end = '', flush = True)
		sort = ["samtools", "sort", "-@", str(threads), "-o", command_base[1] + "/MetaPop/02.Filtered_Samples/"+base_name+"_sorted.bam", sample]
		
		subprocess.call(sort, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		#Cleanup mdz file
		os.remove(sample)
		
		index = ["samtools", "index", "-@", str(threads), command_base[1] + "/MetaPop/02.Filtered_Samples/"+base_name+"_sorted.bam"]
		#Index should naturally produce no stdout or stderr.
		subprocess.call(index)
		
		timer = datetime.now()
		printable_time = timer.strftime(time_format)
		print(" finished at", printable_time)
		
		next_comm = command_base[:]
		next_comm[0] = command_base[1]+ "/MetaPop/02.Filtered_Samples/"+base_name+"_sorted.bam"
		next_comm[1] = command_base[1] + "/MetaPop/02.Filtered_Samples/" + base_name + "_filtered.bam"
		next_comm.append(command_base[1]+ "/MetaPop/")
		next_comm.append(base_name)
		next_comm.append(b_d_base[0])
		next_comm.append(b_d_base[1])
		next_comm.append(b_d_base[3])
		next_comm.append(command_base[1] + "/MetaPop/02.Filtered_Samples/" + base_name + "_preprocessed.bam")
		
		commands.append(next_comm)
		
		index.append([command_base[1] + "/MetaPop/02.Filtered_Samples/"+base_name+"_filtered.bam"])
	
	_inverted_mag_dict = {}
	for file in mag_dict:
		for contig in mag_dict[file]:
			_inverted_mag_dict[contig] = file
	
	#Make these accessible for parallel later.
	#global contig_mag_dict
	#global contig_length_dict
	#global is_mags
	#global inverted_mag_dict
	#global ref_fasta
	
	#inverted_mag_dict = _inverted_mag_dict
	#contig_mag_dict = mag_dict
	#contig_length_dict = length_dict
	is_mags = b_d_base[2]
	ref_fasta = ref_file
	
	procs = min(threads, len(bams))
	#(_cmd, _cld, _im, _imd)
	p = multiprocessing.Pool(procs, initializer=filter_thread_starter, initargs = (mag_dict, length_dict, is_mags, _inverted_mag_dict,))
	
	timer = datetime.now()
	printable_time = timer.strftime(time_format)
	#Sort, filter, depth per pos, breadth and TAD
	print("Filtering genomes. Started at:", printable_time)
	p.map(parse_reads, commands)
	p.close()
	
	timer = datetime.now()
	printable_time = timer.strftime(time_format)
	print("MetaPop preprocessing finished at:", printable_time)
	print(" ")
	
	return None

	#This does %ID and length filter, gets depth per position, and gets breadth and TAD for each genome.
def parse_reads(commands):
	
	in_file = commands[0]
	out_file = commands[1]
	min_length = commands[2] 
	min_pct_id = commands[3] 
	is_global = commands[4]
	dir_base = commands[5]
	base_name = commands[6]
	
	min_cov = commands[7]
	min_depth = commands[8]
	truncation = commands[9]
	
	final_filt = commands[10]
	
	time_format = "%d/%m/%Y %H:%M:%S"
	
	timer = datetime.now()
	printable_time = timer.strftime(time_format)
	print("Filtering reads by length and percent identity for:", base_name, "starting at:", printable_time)
	
	#Just to get the header.
	input_reads = pysam.AlignmentFile(in_file, "rb")
	output_reads = pysam.AlignmentFile(out_file, "wb", template=input_reads)
	
	for read in input_reads:
		parse_entry(read, output_reads, is_global = is_global, min_length = min_length, min_pct_id = min_pct_id)
	
	output_reads.close()
	input_reads.close()
	#Cleanup
	os.remove(in_file)
	os.remove(in_file+".bai")
	
	pysam.index(out_file)
	
	timer = datetime.now()
	printable_time = timer.strftime(time_format)
	print("Getting breadth and depth of coverage for:", base_name, "starting at:", printable_time)
	
	#Depth of coverage per pos
	depth_per_pos_name = dir_base + "04.Depth_per_Pos/" + base_name + "_depth_per_pos.tsv"
	breadth_depth_name = dir_base + "03.Breadth_and_Depth/" + base_name + "_breadth_and_depth.tsv"
	
	current_depths = {}
	count_of_breadth = 0
	last_contig = ""
	
	#contig_tads = {}
	#contig_breadths = {}
	
	depth_per_pos_fh = open(depth_per_pos_name, "w")
	breadth_and_depth_fh = open(breadth_depth_name, "w")
	
	#get a dict with a dict of counts of the appearance of each depth over a genome.
	#Designed to work whether there are multiple genomes or just a single genome per MAG.
	mag_depth_dict = {}
		
	first = True
	
	depth_per_pos = ["samtools", "depth", out_file]
	depths = subprocess.Popen(depth_per_pos, stdout=subprocess.PIPE)
	
	for line in depths.stdout:
		line = line.decode().strip()
		print(line, file = depth_per_pos_fh)
		
		segs = line.split("\t")
		
		contig = segs[0]
		#Select only the readID in case of wrong format.
		contig = contig.split()[0]
		#Pos. in genome is segs[1], but it's not needed here.
		depth = int(segs[2])
		
		#First time through we want to do this.
		#Seems to be skipping the first genome
		if first:
			first = False
			last_contig = contig
		
		
		current_mag = inverted_mag_dict[contig]
		if current_mag not in mag_depth_dict:
			mag_depth_dict[current_mag] = {}
		
		if depth in mag_depth_dict[current_mag]:
			mag_depth_dict[current_mag][depth] += 1
		else:
			mag_depth_dict[current_mag][depth] = 1
				
		#Files are sorted by contig and pos; when we see a new contig, we commit the current one.
		if contig != last_contig:
			#cur_breadth, cur_depth = breadth_and_depth_output(current_depths, contig_length_dict[last_contig], last_contig, truncation, breadth_and_depth_fh)
			breadth_and_depth_output(current_depths, contig_length_dict[last_contig], last_contig, truncation, breadth_and_depth_fh)
			#contig_tads[last_contig] = cur_depth
			#contig_breadths[last_contig] = cur_breadth
			
			#Reset these.
			current_depths = {depth : 1}
			count_of_breadth = 1
			last_contig = contig
		else:
			if depth in current_depths:
				current_depths[depth] += 1
			else:
				current_depths[depth] = 1
			
			count_of_breadth += 1
			last_contig = contig
			
	#Final iteration of loop not done otw.
	#cur_breadth, cur_depth = breadth_and_depth_output(current_depths, contig_length_dict[contig], contig, truncation, breadth_and_depth_fh)
	breadth_and_depth_output(current_depths, contig_length_dict[contig], contig, truncation, breadth_and_depth_fh)
	#contig_tads[contig] = cur_depth
	#contig_breadths[contig] = cur_breadth
	
	depth_per_pos_fh.close()
	breadth_and_depth_fh.close()
	
	
	surviving_genomes = []
	
	#Get the breadth and depth of each 
	for mag in mag_depth_dict:
		mag_len = 0
		for contig in contig_mag_dict[mag]:
			mag_len += contig_length_dict[contig]
		#print("MAG", mag, "is", mag_len, "bp long.")
		breadth, tad = mag_breadth_and_depth(mag_depth_dict[mag], mag_len, truncation)
		
		#If the mag is passing, add its contigs to the list.
		if breadth >= min_cov and tad >= min_depth:
			surviving_genomes.extend(contig_mag_dict[mag])
	
	
	#We don't need to make outputs that have no content, or even look at them.
	if len(surviving_genomes) > 0:
	
		timer = datetime.now()
		printable_time = timer.strftime(time_format)
		print("Filtering by genome for:", base_name, "starting at:", printable_time)
		
		#Previous filtered file
		input_reads = pysam.AlignmentFile(out_file, "rb")
		
		filterbam = pysam.AlignmentFile(final_filt, "wb", template = input_reads)
		
		for read in input_reads:
			ref = input_reads.get_reference_name(read.reference_id)
			if ref in surviving_genomes:
				filterbam.write(read)
				
				
		input_reads.close()
		filterbam.close()
		
		subprocess.call(["samtools", "index", final_filt], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		
	else: 
		print("No genomes survived filtering for:", base_name +". There will be no further output for this file.")
	
	#Cleanup the results of the last step.
	os.remove(out_file+".bai")
	os.remove(out_file)
				
	
	timer = datetime.now()
	printable_time = timer.strftime(time_format)
	print("Completed preprocessing for:", base_name, "at:", printable_time)

	return None
		
def parse_entry(entry, out_handle, is_global, min_length, min_pct_id):

	if not entry.has_tag("MD"):
		return None
	else:
		#Get the match/mismatch info for the read
		mdz_seg = entry.get_tag("MD")
		match_count = re.findall('[0-9]+', mdz_seg)
		
		if is_global:
			#pct ID on global is %match/read legnth
			sum=0
			for num in match_count:
				sum+=int(num)
			#total_count = len(''.join([i for i in mdz_seg if not i.isdigit()])) + sum
			#Total count on global is just read length, not aligned section.
			total_count = entry.query_length
			
			pct_id = (sum/(total_count))*100
			#write if it passes
			if pct_id >= min_pct_id and total_count >= min_length:
				out_handle.write(entry)
		else:
			#pct ID on local is %match/alignment length
			sum=0
			for num in match_count:
				sum+=int(num)
			#Total count here is num matches/alignment length
			total_count = len(''.join([i for i in mdz_seg if not i.isdigit()])) + sum
			pct_id = (sum/(total_count))*100
			#write if it passes
			if pct_id >= min_pct_id and total_count >= min_length:
				out_handle.write(entry)
					
	return None
					
def breadth_and_depth_output(counts_dict, length, name, percentile, handle):
	
	#Update this to use the counts dict, which is depth:count_of_that_depth
	num_pos_covered = 0
	
	for depth in counts_dict:
		num_pos_covered += counts_dict[depth]

	#Gets breadth out of the way, onto getting TAD
	breadth = (num_pos_covered/length) * 100		
	cutoff = math.floor(length/percentile)
	
	if num_pos_covered <= cutoff:
		TAD = 0.0
	else:
		#add zeroes to pad it out and count.
		count_zeroes = (length - num_pos_covered)
		counts_dict[0] = count_zeroes
		
		#print(length, num_pos_covered, end = "")
		
		#Unique observed depths
		observed_depths = list(counts_dict.keys())
		#sort ascending
		observed_depths.sort()
		
		
		cutoff = math.floor(length/percentile)
		
		#Remove lower values
		depths_removed_so_far = 0
		for depth in observed_depths:
			observations_of_current_depth = counts_dict[depth]
			#This occurs once we're in the right interval of depths to stop truncating.
			if depths_removed_so_far + observations_of_current_depth >= cutoff:
				counts_dict[depth] -= (cutoff - depths_removed_so_far)
				break
			else:
				depths_removed_so_far += observations_of_current_depth
				counts_dict[depth] = 0
			
		#remove upper values
		depths_removed_so_far = 0
		#Reverse sort so we go from highest -> lowest depth.
		observed_depths.sort(reverse = True)
		for depth in observed_depths:
			observations_of_current_depth = counts_dict[depth]
			#This occurs once we're in the right interval of depths to stop truncating.
			if depths_removed_so_far + observations_of_current_depth >= cutoff:
				counts_dict[depth] -= (cutoff - depths_removed_so_far)
				break
			else:
				depths_removed_so_far += observations_of_current_depth
				counts_dict[depth] = 0
		
		#get TAD
		TAD = 0
		for depth in counts_dict:
			TAD += (depth * counts_dict[depth])
		
		#Average depth over included positions - since we truncate the first and last (cutoff) positions, this works.
		total_count = length - (2*cutoff)
		TAD = round(TAD/total_count, 6)

	print(name, num_pos_covered, breadth, TAD, sep = "\t", file = handle)
	
	return None
	#return breadth, TAD

def mag_breadth_and_depth(counts_dict, length, percentile):
	#Update this to use the counts dict, which is depth:count_of_that_depth
	num_pos_covered = 0
	
	for depth in counts_dict:
		num_pos_covered += counts_dict[depth]

	#Gets breadth out of the way, onto getting TAD
	breadth = (num_pos_covered/length) * 100		
	cutoff = math.floor(length/percentile)
	
	if num_pos_covered <= cutoff:
		TAD = 0.0
	else:
		#add zeroes to pad it out and count.
		count_zeroes = (length - num_pos_covered)
		counts_dict[0] = count_zeroes
		
		#print(length, num_pos_covered, end = "")
		
		#Unique observed depths
		observed_depths = list(counts_dict.keys())
		#sort ascending
		observed_depths.sort()
		
		cutoff = math.floor(length/percentile)
		
		#Remove lower values
		depths_removed_so_far = 0
		for depth in observed_depths:
			observations_of_current_depth = counts_dict[depth]
			#This occurs once we're in the right interval of depths to stop truncating.
			if depths_removed_so_far + observations_of_current_depth >= cutoff:
				counts_dict[depth] -= (cutoff - depths_removed_so_far)
				break
			else:
				depths_removed_so_far += observations_of_current_depth
				counts_dict[depth] = 0
			
		#remove upper values
		depths_removed_so_far = 0
		#Reverse sort so we go from highest -> lowest depth.
		observed_depths.sort(reverse = True)
		for depth in observed_depths:
			observations_of_current_depth = counts_dict[depth]
			#This occurs once we're in the right interval of depths to stop truncating.
			if depths_removed_so_far + observations_of_current_depth >= cutoff:
				counts_dict[depth] -= (cutoff - depths_removed_so_far)
				break
			else:
				depths_removed_so_far += observations_of_current_depth
				counts_dict[depth] = 0
		
		#get TAD
		TAD = 0
		for depth in counts_dict:
			TAD += (depth * counts_dict[depth])
		
		#Average depth over included positions - since we truncate the first and last (cutoff) positions, this works.
		total_count = length - (2*cutoff)
		TAD = round(TAD/total_count, 6)
	
	return breadth, TAD


def fill_mdz(args):
	file, base_name, ref, dir = args[0], args[1], args[2], args[3]
	time_format = "%d/%m/%Y %H:%M:%S"
	timer = datetime.now()
	printable_time = timer.strftime(time_format)
	print("Preparing percent identity for:", base_name, "starting at", printable_time)
	
	output = os.path.normpath(dir + "/MetaPop/02.Filtered_Samples/"+base_name+"_filled_MD.bam")
	
	#sam_to_bam = ["samtools", "view", "-h", "-u", file]
	#Make sure the file has filled MD:Z: tags - it usually will already, this is the default for most read alignment tools now.
	calmd = ["samtools", "calmd", "-u", "-b", file, ref]
	
	md_to_bam = ["samtools", "view", "-h", "-b", "-", "-o", output]
	
	md = subprocess.Popen(calmd, stdout = subprocess.PIPE)
	compress = subprocess.Popen(md_to_bam, stdin=md.stdout)
	
	md.wait()
	compress.wait()
	
	os.remove(file)
	os.remove(file+".bai")
	
	timer = datetime.now()
	printable_time = timer.strftime(time_format)
	print("Percent identity prepared for:", base_name, "finished at", printable_time)
	
	return output
