#!/usr/bin/env python

"""This preprocessing script is part of the CLAM pipeline.

This subcommand (new v1.1) will prepare the input files for CLAM pipeline. As of the current version (v1.1), it looks for 
reads passing QC, splits the input bam file by sorting them into `unique.sorted.bam` and `multi.sorted.bam`, 
and adding an additional tag "RT" (short for Read Tag) to each alignment based which read tagger function the user supplied.

Note that you can also run `CLAM realigner` directly, which will call `preprocessor` and automatically determine
if `preprocessor` has been called in the output folder. 

If you don't want to run `realigner`, you can also run `peakcaller` directly after `preprocessor`.

Example run:
	```
	CLAM preprocessor -i path/to/input/Aligned.out.bam -o path/to/clam/outdir/ --read-tagger-method median
	```
Author:
	Zijun Zhang <zj.z@ucla.edu>

Tested under python 2.7
"""
# from . import config
# __version__ = config.__version__

import os
import sys
import pysam
import numpy as np
from collections import defaultdict
#from tqdm import tqdm
import logging
import datetime
import bisect
import argparse as ap
import inspect


logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('') # preprocessor.parsering


def alignment_mutation(x, mut_ref, mut_obs):
	"""DOCSTRING
	Need to read reference genome
	NotImplemented
	"""	
	raise NotImplementedError()


def read_tagger_collection(alignment, method='median', **kwargs):
	""" tag a read alignment to a genomic locus
	Args:
	Returns:
	"""
	tagger_func = {
		# center of the read; must dicard junction reads
		'median': lambda x: -1 if 'N' in x.cigarstring else int(np.median(x.positions))+1,
		# 5' start site of the read; truncation in iCLIP/eCLIP
		'start': lambda x: -1 if 'N' in x.cigarstring else x.positions[-1]+1 if x.is_reverse else x.positions[0]+1,
		# extend from 5' site to certain length; need kwargs
		'extend': lambda x: -1 if 'N' in x.cigarstring else x.positions[-1]-kwargs['ext_len'] if x.is_reverse else x.positions[0]+kwargs['ext_len'],
		# mutation tag a specific mutation type
		'mutation': lambda x: alignment_mutation(x, kwargs['mut_ref'], kwargs['mut_obs'])
		}
	try:
		tag=tagger_func[method](alignment)
	except:
		tag=-1
	return tag



def filter_bam_multihits(filename, max_tags, max_hits, out_dir, read_tagger_method, strandness):
	"""Pre-processing function for cleaning up the input bam file.
	Args:
	Returns:
	"""
	# logging the parameter values
	frame = inspect.currentframe()
	args, _, _, values = inspect.getargvalues(frame)
	msg = 'Params:\n'
	for i in args:
		msg += "%s = %s \n"%(i, values[i])
	logger.info(msg)
	read_tagger=lambda x: read_tagger_collection(x, method=read_tagger_method)
	logger.info('filtering input bam')
	
	in_bam = pysam.Samfile(filename,'rb')
	# unique read bam
	ubam_fn = os.path.join(out_dir, 'unique.bam')
	sorted_ubam_fn = os.path.join(out_dir, 'unique.sorted.bam')
	ubam=pysam.Samfile(ubam_fn, 'wb', template=in_bam)
	unique_counter = 0
	
	# multi-read bam
	mbam_fn = os.path.join(out_dir, 'multi.bam')
	sorted_mbam_fn = os.path.join(out_dir, 'multi.sorted.bam')
	mbam=pysam.Samfile(mbam_fn, 'wb', template=in_bam)
	mread_set = set()
	
	if "STAR" in [x['ID'] for x in mbam.header['PG']]:
		logger.info('input bam aligner is STAR')
		tool = "STAR" 
	else: # "bowtie2" in [x['ID'] for x in mbam.header['PG']]:
		logger.info('input bam aligner is bowtie2')
		tool = "bowtie2"
	#else:
	#	logger.info('input bam aligner is not detected, use all align records as uniq-align')
	#	tool = "nul"
#if read.is_paired:
	# splitting unique and multi- reads
	# and add the read taggers we need
	if not \
		(os.path.isfile( os.path.join(out_dir,'unique.sorted.bam') ) and \
		os.path.isfile( os.path.join(out_dir,'multi.sorted.bam')) ):
			
		#for read in tqdm(in_bam):
		counter = 0
		for read in in_bam:
			# poor man's progress bar
			counter += 1
			if not counter % 10**6:
				logger.info('tagged %i alignments'%counter)
			read_tag = read_tagger(read)
			## skip reads with unassigned tagger
			if read_tag==-1:
				continue
			read.tags += [('RT', read_tag)] ## add the tag

			tagged_read = pysam.AlignedSegment()
			tagged_read.query_name = read.query_name
			tagged_read.query_sequence = 'N'
			tagged_read.flag = read.flag
			tagged_read.reference_id = read.reference_id
			tagged_read.reference_start = read_tag - 1  # 0-based leftmost coordinate
			tagged_read.mapping_quality = read.mapping_quality
			tagged_read.cigar = ((0,1),)
			tagged_read.template_length = read.template_length
			tagged_read.query_qualities = pysam.qualitystring_to_array("<")

			tmp1 = read.next_reference_start + 1 # convert 0 to 1-based (unmap: -1 --> 0)
			if tmp1 != 0: tmp2 = in_bam.get_reference_name(read.next_reference_id) # SE
			else: tmp2 = "*"  
			tagged_read.next_reference_id = read.next_reference_id # original treated as SE, and omited
			tagged_read.next_reference_start = read.next_reference_start # original treated as SE, and omited

			tagged_read.tags = read.tags
			read_len = sum([i[1] for i in read.cigar if i[0]==0])
			tagged_read.tags += [('RL', read_len)]

			#keep original read info to tag
			tagged_read.tags += [('os', read.query_sequence)]
			tagged_read.tags += [('or', read.reference_start+1)] # 0-based leftmost coordinate
			tagged_read.tags += [('og', read.cigarstring)]
			tagged_read.tags += [('oq', pysam.qualities_to_qualitystring(read.query_qualities))] 
			

			tagged_read.tags += [('oa', tmp2)]
			tagged_read.tags += [('ob', tmp1)]

			# add strandness check
			if strandness != "none":
				tagged_read.is_reverse = (read.is_reverse) ^ (strandness!="same")
			
			#https://www.seqanswers.com/forum/bioinformatics/bioinformatics-aa/21213-bowtie2-uniquely-mapped-reads
			if ( ( tool == "bowtie2" ) and (read.is_secondary or read.has_tag('XS')) ):
				#try:
				# if read.opt("NH") < max_hits: # already filtered in bowtie2
				mbam.write(tagged_read)
				mread_set.add(read.qname)
				#except KeyError:
				#	#print read
				#	raise Exception('%s: missing NH tag when is_secondary=%s'%(read.qname,read.is_secondary))
			elif ( tool == "STAR" ) and (read.is_secondary or (read.has_tag('NH') and read.opt("NH")>1)):
				# logger.info('input bam aligner is STAR')
				if read.opt("NH") < max_hits:
					mbam.write(tagged_read)
					mread_set.add(read.qname)
			else:
				ubam.write(tagged_read)
				unique_counter += 1
		
		ubam.close()
		mbam.close()
		
		# sorting
		pysam.sort('-m', '4G', '-@', '3', '-T', os.path.dirname(sorted_ubam_fn), '-o', sorted_ubam_fn, ubam_fn)
		os.remove(ubam_fn)
		pysam.sort('-m', '4G', '-@', '3', '-T', os.path.dirname(sorted_mbam_fn), '-o', sorted_mbam_fn, mbam_fn)
		os.remove(mbam_fn)
		pysam.index(sorted_ubam_fn)
		pysam.index(sorted_mbam_fn)
		
		# log the statistics
		multi_counter = len(mread_set)
		# logger.info('testtesttesttesttesttesttesttesttesttesttesttest')
		logger.info('Unique reads = %s;  ' %unique_counter + 'Multi reads = %s (%.2f %%)' %( multi_counter, float(multi_counter)/(multi_counter+unique_counter)*100 ))
	else:
		logger.info('found previously sorted tag-bam. checking if need collapsing.')
	
	# filter redundant tags if turned on
	if max_tags>0:
		logger.info('collapsing unique')
		filter_bam_maxtags(os.path.join(out_dir, 'unique.sorted.collapsed.bam'), os.path.join(out_dir, 'unique.sorted.bam'), max_tags)
		logger.info('collapsing multi')
		filter_bam_maxtags(os.path.join(out_dir, 'multi.sorted.collapsed.bam'), os.path.join(out_dir, 'multi.sorted.bam'), max_tags)
	
	in_bam.close()
	return


def collapse_stack(stack, collapse_dict, max_tags):
	"""DOCSTRING
	Args
	Returns
	"""
	new_alignment_list = []
	new_alignment_dict = defaultdict(list)
	for aln in stack:
		new_alignment_dict[aln.query_sequence].append(aln)
	
	# TODO 2017.10.21: 
	# further collapse `new_alignment_dict`
	# based on degeneracy and/or read tags
	
	for seq in new_alignment_dict:
		this_alignment_qname_list = [x.qname for x in new_alignment_dict[seq] ]
		is_collapsed = [True if x in collapse_dict else False for x in this_alignment_qname_list]
		## if any of the alignment is collapsed before,
		## we require all of them to be collapsed
		if any(is_collapsed):
			assert all(is_collapsed)
			target_alignment_qname = collapse_dict[this_alignment_qname_list[0]][0:max_tags]
			assert len(collapse_dict[this_alignment_qname_list[0]]) <= max_tags
			target_alignment = [new_alignment_dict[seq][this_alignment_qname_list.index(x)] for x in target_alignment_qname]
		else:
			target_alignment = new_alignment_dict[seq][0:max_tags]
			for aln_qname in this_alignment_qname_list:
				collapse_dict[aln_qname] = [x.qname for x in target_alignment]
		new_alignment_list.extend( target_alignment )
	return new_alignment_list, collapse_dict


def filter_bam_maxtags(obam_fn, ibam_fn, max_tags=1):
	"""DOCSTRING
	Args
	Returns
	"""
	assert max_tags>0
	# prepare files
	ibam = pysam.Samfile(ibam_fn, 'rb')
	obam = pysam.Samfile(obam_fn, 'wb', template=ibam)
	# init 
	collapse_dict = defaultdict(list)
	chr_list=[x['SN'] for x in ibam.header['SQ']]
	input_counter = 0
	output_counter = 0
	
	for chr in chr_list:
		# empty stack for each new chromosome
		stack = []
		last_pos = -1
		for read in ibam.fetch(chr):
			input_counter += 1
			if not (input_counter % (5*(10**6)) ):
				logger.info('collapsed %i alignments'%input_counter)
			if read.positions[0] > last_pos:
				new_alignment_list, collapse_dict = collapse_stack(stack, collapse_dict, max_tags)
				output_counter += len(new_alignment_list)
				last_pos = read.positions[0]
				stack = [read]
				for new_alignment in new_alignment_list:
					#new_alignment.query_sequence = '*'
					#new_alignment.query_qualities = '0'
					_ = obam.write(new_alignment)
			else:
				stack.append(read)
		new_alignment_list, collapse_dict = collapse_stack(stack, collapse_dict, max_tags)
		output_counter += len(new_alignment_list)
		last_pos = read.positions[0]
		for new_alignment in new_alignment_list:
			#new_alignment.query_sequence = '*'
			#new_alignment.query_qualities = '0'
			_ = obam.write(new_alignment)
	ibam.close()
	obam.close()
	#os.rename(obam_fn, ibam_fn)
	#pysam.sort(obam_fn)
	pysam.index(obam_fn)
	logger.info('Input = %s; Output = %s; Redundancy = %.2f'%(input_counter,output_counter, 1-float(output_counter)/input_counter))
	return
	
	


def parsering(args):
	"""DOCSTRING
	Args
	Returns
	"""
	try:
		in_bam = args.in_bam
		out_dir  = args.out_dir
		if not (os.path.isdir(out_dir) or os.path.exists(args.out_dir)):
			os.makedirs(out_dir) # recursively
		tag_method = args.tag_method
		max_hits = args.max_hits
		## Note: if specified max_tags, need pre-sorted bam
		max_tags = args.max_tags
		strandness = args.strandness
		
		#logger = logging.getLogger('CLAM.Preprocessor')
		logger.info('start')
		logger.info('run info: %s'%(' '.join(sys.argv)))
		
		filter_bam_multihits(in_bam, max_hits=max_hits, max_tags=max_tags, out_dir=out_dir, 
			read_tagger_method=tag_method,
			strandness=strandness)
		
		logger.info('end')
	except KeyboardInterrupt():
		sys.exit(0)
	return



if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser('preprocessor')

	# ag_prep = subparsers.add_parser("preprocessor", help="CLAM Preprocessor: tag read alignments to specific locus")

	# input/output
	parser.add_argument("-i", "--input", dest="in_bam", type=str, required=True, help="Input bam file")
	parser.add_argument("-o", "--out-dir", dest="out_dir", type=str, required=True, help="Output folder")

	# processing
	parser.add_argument("--read-tagger-method", dest="tag_method", type=str, choices= ('median', 'start'), default='median', help="Read tagger method, 'median' for read center, 'start' for read start site; default: median")
	parser.add_argument("--max-multihits", dest="max_hits", type=int, default=100, help="The maximum hits allowed for multi-mapped reads; default: 100")
	parser.add_argument("--max-tags", dest="max_tags", type=int, default=-1, help="The maximum identical tags at given location; default: -1, no filter")

	# strandness
	parser.add_argument("--strandness", dest="strandness", type=str, default="same", choices=['same', 'opposite', 'none'], help="The expected read strandness with transcription direction: same, opposite, or none(i.e. unstranded); default: same")

	# parser = argparse.ArgumentParser('BigWig library')
	# parser.add_argument('--input-file', '-i', type=str, required=True, help='input BigWig file')
	# parser.add_argument('--library', '-l', type=str, required=True, help='path to libjkweb.so')
	args = parser.parse_args()

	parsering(args)

