#! /usr/bin/env python
import argparse, sys, os, errno
from distutils.spawn import spawn
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
# sys.path.append("/BioII/lulab_b/baopengfei/gitsoft/clipper")
# from distutils.extension import Extension
# peaks = Extension("clipper.src.peaks", sources = ['clipper/src/peaksmodule.cc'])
# from src.peaks import shuffle
# from src.readsToWiggle import readsToWiggle_pysam
import numpy as np
import datetime

'''
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
sns.set()
'''
import pandas as pd
from pandas import DataFrame, Series
from scipy.fftpack import fft
from scipy.signal import convolve
# import numba

import pyBigWig
import pysam
from collections import defaultdict
# from statsmodels.sandbox.stats.multicomp import multipletests
import math, random, re
from multiprocessing import Pool

command_handlers = {}
def command_handler(f):
    command_handlers[f.__name__] = f
    return f


#accelrate for loop in func. : https://blog.csdn.net/qq_43657442/article/details/110000940
#@numba.jit('int64(int32[:], int32[:], float64, float64, float64)')

#  consider read overlapped with range
def count_pileup_heights(tlen, reads):
    """
    Sub-routine for do_permutation(...)
    Counts the distribution of pile-up heights for a given gene/permutation
    """
    loc_heights=[0] * tlen
    for id, pos, read_len, score in reads:
        for i in range(max(0,pos), min(pos+read_len, tlen)): # range(0,10) ~ range(10) ~ 0,1,2...9      0-based
            loc_heights[int(i)]+=score
    return loc_heights

## only consider full read inside range
# def count_pileup_heights(tlen, reads):
#     """
#     Sub-routine for do_permutation(...)
#     Counts the distribution of pile-up heights for a given gene/permutation
#     """
#     loc_heights=[0] * tlen
#     for id, pos, read_len, score in reads:
#         for i in range(pos, min(pos+read_len, tlen)): # range(0,10) ~ range(10) ~ 0,1,2...9      0-based
#             loc_heights[int(i)]+=score
#     return loc_heights

#  consider read overlapped with range
def permutate_heights(tlen, reads):
    """
    Sub-routine for do_permutation(...)
    Randomly allocate the read locations.
    """
    loc_heights=[0] * tlen
    for id, pos, read_len, score in reads:
        # if score<1 and random.random() > score:
        #     continue
        # rand_pos=random.randint(0, max(1, tlen-read_len)) # only consider full read inside range
        rand_pos=random.randint(-read_len, max(1, tlen)) # consider read overlapped with range

        # for i in range(rand_pos, min(rand_pos + read_len, tlen)): # only consider full read inside range
        for i in range(max(0,rand_pos), min(rand_pos + read_len, tlen)): # consider read overlapped with range    
            loc_heights[int(i)]+=1
    # print(datetime.datetime.now(),"loc_heights: ", loc_heights)
    return loc_heights

## only consider full read inside range
# def permutate_heights(tlen, reads):
#     """
#     Sub-routine for do_permutation(...)
#     Randomly allocate the read locations.
#     """
#     loc_heights=[0] * tlen
#     for id, pos, read_len, score in reads:
#         # if score<1 and random.random() > score:
#         #     continue
#         rand_pos=random.randint(0, max(1, tlen-read_len)) # only consider full read inside range
#         # rand_pos=random.randint(-read_len, max(1, tlen)) # consider read overlapped with range

#         for i in range(rand_pos, min(rand_pos + read_len, tlen)): # only consider full read inside range
#         # for i in range(max(0,rand_pos), min(rand_pos + read_len, tlen)): # consider read overlapped with range    
#             loc_heights[int(i)]+=1
#     # print(datetime.datetime.now(),"loc_heights: ", loc_heights)
#     return loc_heights

def do_permutation(transcr, read_transcript, max_iter, pval_cutoff, min_cov, seed):	
    """
    Permutes the reads along a given gene length, sub-routine that get called by get_permutation_fdr(..).
    Returns the locally corrected p-values for each observed height on the given gene.
    """
    # global fail_sig
    # fail_sig = False
    
    random.seed(seed)
    # if you have very few reads on a gene, don't waste time
    if len(read_transcript) <= 5:
        res = []
        res.append(min_cov)
        res.append({})
        print(datetime.datetime.now(),"do_permutation len(read_transcript) <= 5, early return min_cov !")
        return res

    chr, tstart, tend, strand, tid = transcr[0:5]
    tid_length = tend-tstart # +1 ?
    obs_heights_count=count_pileup_heights(tlen=tid_length, reads=read_transcript)
    rand_heights_dist=defaultdict(int)
    rand_sum=0
    # need to account for the 'observed' data, since permutation tests should never report p-value as 0. 3/22/16
    # and also make sure the cov value always included in height_to_pval{}
    for i in obs_heights_count:
        if i==0:
            continue
        else:
            rand_heights_dist[int(i)]+=1
            rand_sum+=1
    for B in range(max_iter):
        new_heights_count=permutate_heights(tlen=tid_length, reads=read_transcript)
        for i in new_heights_count:
            if i==0:
                continue
            else:
                rand_heights_dist[int(i)]+=1
                rand_sum+=1
    height_to_pval={}
    for h in set(rand_heights_dist.keys()):
        if h < 1: # skip 
            continue # zero(<1) truncated ?
        else:
            lefter=0
            for j in range(int(h), max(rand_heights_dist)+1):
                lefter+=rand_heights_dist[j]
            height_to_pval[int(h)]=lefter/float(rand_sum) # py2: /; py3: //
    # def min_sig(vals,p=0.05):
    #     vals_sort = sort(vals)
    
    # print(datetime.datetime.now(),"do_permutation obs_heights_count:", obs_heights_count)
    # print(datetime.datetime.now(),"do_permutation len(obs_heights_count):", len(obs_heights_count))
    # print(datetime.datetime.now(),"do_permutation rand_heights_dist:", rand_heights_dist)
    # print(datetime.datetime.now(),"do_permutation height_to_pval:", height_to_pval)
    
    # too many dist model options reported, we just select cov that just pass perm_fdr
    not_sig = {key:val for key,val in height_to_pval.items() if val > pval_cutoff} #  need select h whose p >0.05, get err if dict is blank
    try:
        h0_tmp = min(not_sig, key=lambda k: not_sig[k]) # return min key
    except:
        # when not_sig = {} null dict (means no cov meets p<pval_cutoff) or other situation ?
        h0_tmp = int(np.percentile(list(rand_heights_dist.keys()), (95)))
        # h0_tmp = max(rand_heights_dist, key=lambda k: rand_heights_dist[k])
        print(datetime.datetime.now(),"do_permutation exception !")
    # pval_list=[]
    # for i in obs_heights_count:
    # 	if i<1: # skip here: rm too many <1 cov sites ?
    # 		pval_list.append(1.0) #continue # why not append 1 as pvalue: pval_list.append(1.0)
    # 	pval_list.append(height_to_pval[i])
    # if len(pval_list)<=1:
    # 	return []
    # qval_list=multipletests(pval_list, method='fdr_bh')[1]
    if h0_tmp < min_cov:
        h0_tmp = min_cov
    h0_tmp = int(h0_tmp)
    # print(datetime.datetime.now(),h0_tmp)
    # print(datetime.datetime.now(),height_to_pval)
    res = []
    res.append(h0_tmp)
    res.append(height_to_pval)
    # res.append(fail_sig)
    # print(datetime.datetime.now(),"do_permutation fail_sig:", fail_sig)
    return res

def find_one_maxima_in_center_local(arr):
    """
    Returns a list of boolean values for an array that mark if a value is a local
    maxima or not True for yes false for no
    Importantly for ranges of local maxima the value in the middle of the range
    is chosen as the minimum value
    """
    # import numpy as np
    # to initalize a new array to all false
    maxima = np.empty(len(arr)) # , dtype='bool'
    maxima.fill(0)
    print("array:",arr)
    max_range_start = 0
    increasing = True
    for i in range(len(arr[:-1])): # 

        # update location of maxima start until
        if arr[i] < arr[i + 1]:
            max_range_start = i + 1
            increasing = True

        if (arr[i] > arr[i + 1]) and increasing is True:
            increasing = False
            # gets the local maxima midpoint
            midpoint = int((max_range_start + i) / 2) # math.floor
            maxima[int(midpoint)] = arr[int(midpoint)] # True

    # catches last case
    if increasing:
        midpoint = int((max_range_start + len(arr) ) / 2) # math.floor, (max_range_start + len(arr) - 1)
        maxima[int(midpoint)] = arr[int(midpoint)] # True
    
    # find one most max cov pos in center
    indexed_maxima = {i:val for i, val in enumerate(maxima) if val}
    maxima_pos = max(indexed_maxima, key=lambda k: indexed_maxima[k]) # return key of max value (first key if multi max exist)
    maxima_cov = indexed_maxima[maxima_pos]
    return maxima_pos, int(maxima_cov)

def read_tid_frag_from_bam(tid, reads, is_unique):
    """
    Use pysam to fetch reads info for a given gene and its loci.
    Returns reads, read weights and its mapped loci.
    is_stranded
    """
    tid_reads=[]
    #gene, chr, strand, start, end=tid
    chr, start, end, strand, gene = tid[0:5]
    # if strand=='-':
    #     is_reverse=True
    # else:
    #     is_reverse=False    
    # reads=[x for x in bam.fetch(chr, int(start), int(end)) if x.is_reverse==is_reverse or not is_stranded]
    reads=[x for x in reads if x.pos>=int(start) and x.positions[-1]<=int(end)] # x.pos, x.positions: 0-based; tid: 1-based 
    for read in reads:
        if is_unique:
            # try:
            #     opt_NH=read.opt('NH')
            #     if opt_NH > 1:
            #         continue
            # except:
            #     pass
            score=1
        # else:
        #     try:
        #         opt_AS=read.opt('AS')
        #         if isinstance(opt_AS, float):
        #             score=opt_AS
        #         else:
        #             continue
        #     except:
        #         continue
        try:
            read_length = read.opt('RL')
        except:
            read_length = read.positions[-1] - read.positions[0] + 1 # new pysam: get_reference_positions
        
        if (read.positions[0]-start>=0 and end-read.positions[-1]>=0):    #to aviod long read:  and read_length<200; to avoid junction reads: (not 'N' in read.cigarstring), remove these to keep in line with bw signal 
            tid_reads.append([read.qname, read.pos-start, read_length, score])
    return tid_reads # [["a",1,2],["b",2,3]...]
def read_tid_frag_from_records(tid, reads):
    """
    Use reads records [read.qname, read.pos-start, read_length, score] info for a given gene and its loci.
    Returns reads, read weights and its mapped loci.
    is_stranded, is_unique
    """
    tid_reads=[]
    #gene, chr, strand, start, end=tid
    chr, start, end, strand, gene = tid[0:5]
    reads=[x for x in reads if (x[1]+x[2])>=int(start) and x[1]<=int(end)] # include all overlapped reads, not just full reads inside, to aviod tx_cov/localmax not included in _height_to_pval, should result to fewer nan pvalue
    for read in reads:
        tid_reads.append([read[0], read[1]-start, read[2], read[3]])
    return tid_reads # [["a",1,2],["b",2,3]...]


def _call_peaks_localmax(tx_id, count_signal, tx_reads, min_peak_length=10, max_peak_length=200, bin_width=10, min_cov=5, permutate_pval=0.05, decay=0.05, local=True, max_iter=100, seed=1234):
    '''
    input: bw signals, 
    return: local or global maxima position
    '''
    # cal all tx reads aver length 
    # read_lengths = 0
    # for read in reads:
    #     try:
    #         read_length = read.opt('RL')
    #     except:
    #         read_length = read.positions[-1] - read.positions[0] + 1
    #     read_lengths += read_length
    # aver_read_len = np.mean(np.array(read_lengths))
    #logger(aver_read_len)

    # shuffle tx reads    
    # convert pysam to a wiggle vector, junction, positional count(coverage), read lengths, all_reads, location    
    # global _height_to_pval
    # global tx_cov
    
    x = np.array(count_signal)
    gene_length = x.shape[0] # gene_length = len(x)
    # (wiggle, jxns, pos_counts,lengths, allreads, read_locations) = readsToWiggle_pysam(tx_reads, 0, gene_length, "+", "center", False)
    # # gene_length = int(lengths)
    # aver_read_len = np.mean(np.array(lengths))
    # #logger
    # lengths = [gene_length - 1 if read >= gene_length else read for read in lengths]
    # tx_cov = get_FDR_cutoff_mean(readlengths=lengths,genelength=gene_length,alpha=0.05)   
    # print(datetime.datetime.now(),"\t","original tx reads record: ", tx_reads)
    tx_global = [tx_id,0,gene_length,"+",tx_id]
    tx_reads_global = tx_reads_local = read_tid_frag_from_records(tx_global, tx_reads) # read_tid_frag_from_records(tx_global, tx_reads) # get converted tx reads 
    print(datetime.datetime.now(),"\t","full tx reads num: ", len(tx_reads_global))
    # print(datetime.datetime.now(),"\t","full tx reads record: ", tx_reads_global)
    #tx_cov = _child_get_permutation_fdr(tx_global, tx_reads_global, pval_cutoff=permutate_pval, min_cov=min_cov, max_iter=max_iter, seed=seed)
    res = do_permutation(transcr=tx_global, read_transcript=tx_reads_global, max_iter=max_iter, pval_cutoff=permutate_pval, min_cov=min_cov, seed=seed)
    tx_cov, _height_to_pval = res[0:2]
    tx_cov = max(tx_cov, min_cov)
    print(datetime.datetime.now(),"\t","full tx coverage: ", tx_cov)
    # print(datetime.datetime.now(),"\t","full tx _height_to_pval: ", _height_to_pval)
    # _, _, rlen, _ = tx_reads_global[0:4]
    rlen = [x[2] for x in tx_reads_global]
    aver_read_len = int(np.mean(np.array(rlen)))
    # print(datetime.datetime.now(),"tx aver read length: ", aver_read_len)
    
    # 0. average signal over bins with 50% overlap, and filter low cov bins (<min_cov)
    half_bin_width = bin_width//2
    n_bins = max(1, gene_length//half_bin_width)
    bin_cov = np.zeros(n_bins)
    for i in range(n_bins):
        bin_cov[i] = np.mean(x[(i*half_bin_width):min((i + 2)*half_bin_width, gene_length)])
    
    cand_bins = np.nonzero(bin_cov > min_cov)[0] # return 0-based half-bin index in bin_cov
    n_cand_bins = cand_bins.shape[0]
    cand_bin_index = 0
    left_bound = 0
    peaks = []
    while (cand_bin_index < n_cand_bins):
        # 1. find local max pos and value (use fixed hard threshold cov: min_cov, eg. 5)
        i = (cand_bins[cand_bin_index]+1)*half_bin_width # seem +1
        if ((cand_bins[cand_bin_index]+1)*half_bin_width >= gene_length):
            i = gene_length-1
        print("\n",datetime.datetime.now(),"\t\t","start bin_center: ",i)
        start = i
        end = i
        start_trend = i
        end_trend = i
        # max_index_tmp = i
        
        # find local max, not limit to this bin (bin center)
        while (start > left_bound) and (x[start - 1] >= x[start]) and (x[start - 1] >= min_cov): 
            start -= 1
            #start already -1
            if (x[start] > x[start+1]): 
                start_trend = start 
        # most_left_start = start
        start = math.ceil((start_trend + start) / 2)
        while (end < (gene_length - 1)) and (x[end + 1] >= x[end]) and (x[end + 1] >= min_cov): 
            end += 1
            #end already +1
            if (x[end] > x[end-1]): 
                end_trend = end 
        # most_right_end = end
        end = math.floor((end_trend + end) / 2)
        
        # max_index = 0
        if x[start] >= x[end]:
            local_max = x[start]
            max_index = start
        else:
            local_max = x[end]
            max_index = end
        # max_index_tmp = max_index

        if x[max_index] <= min_cov:
            print(datetime.datetime.now(),"\t\t","initial cov <= min_cov, skip this bin !")
            cand_bin_index += 1
            continue
        

        # 2. define shuffle local region (use full-len tx/global if tx len shorter than 2*ext)
        #s: local start site; end: local end site; ext: bp extension each side, usually > 2*aver_read_len (eg.50)
        # if local:
        # seem cannot aviod re-cal range (localmax center)
        most_left_start = max_index
        while (most_left_start > left_bound) and (x[most_left_start - 1] >= min_cov): 
            most_left_start -= 1
        most_right_end = max_index
        while (most_right_end < (gene_length - 1)) and (x[most_right_end + 1] >= min_cov): 
            most_right_end += 1

        # if False:
        #     tmp_ext = int(2*max(most_left_start-max_index, most_right_end-max_index))
        #     ext = int(max(tmp_ext, aver_read_len, 15)) 
        #     s = max_index - ext 
        #     e = max_index + ext
        # else:
        #    ext = int(max(0.5*(most_right_end-most_left_start), 10)) 

        span = 2*(most_right_end-most_left_start)
        if span <= 10:
            print(datetime.datetime.now(),"\t\t","initial span <= 10, skip this bin !")
            cand_bin_index += 1
            continue
        
        # s = most_left_start-(max_index-most_left_start)
        # e = most_right_end+(most_right_end-max_index)
        s = int(0.5*(most_left_start+most_right_end)-(most_right_end-most_left_start))
        e = int(s+span)
        print(datetime.datetime.now(),"\t\t","1st perm start pos: ", s)
        print(datetime.datetime.now(),"\t\t","1st perm end pos: ", e)
        
        print(datetime.datetime.now(),"\t\t","tx max_index: ", max_index)
        print(datetime.datetime.now(),"\t\t","tx max_index - most_left_start: ", max_index-most_left_start)
        print(datetime.datetime.now(),"\t\t","tx most_right_end - max_index: ", most_right_end-max_index)
        print(datetime.datetime.now(),"\t\t","tx RL: ", aver_read_len)
        print(datetime.datetime.now(),"\t\t","tx span: ", span)
        print(datetime.datetime.now(),"\t\t",cand_bin_index,"/",n_cand_bins,"localmax:",local_max)
            
        if gene_length <= span:
            s = 0
            e = gene_length
            print(datetime.datetime.now(),"\t\t",cand_bin_index,"/",n_cand_bins,"local mode beyound boundary, use full-len !")
            print(datetime.datetime.now(),"\t\t","tx cov: ", tx_cov)
        else:
            if s < 0: 
                s = 0 
                e = span
                print(datetime.datetime.now(),"\t\t",cand_bin_index,"/",n_cand_bins,"left beyound boundary !")
            elif e > gene_length:
                s = gene_length - span
                e = gene_length
                print(datetime.datetime.now(),"\t\t",cand_bin_index,"/",n_cand_bins,"right beyound boundary !")
            print(datetime.datetime.now(),"\t\t",cand_bin_index,"/",n_cand_bins,"local mode not beyound boundary !")
            
            # 3. shuffle reads at local region
            # tx_id = tx_reads.reference[0]
            tx_local = [tx_id,s,e,"+",tx_id]
            tx_reads_local = read_tid_frag_from_records(tx_local, tx_reads) # get converted tx reads, iteration index at last pos ? 
            print(datetime.datetime.now(),"\t\t",cand_bin_index,"/",n_cand_bins,"local reads number: ", len(tx_reads_local))
            # print(datetime.datetime.now(),cand_bin_index,"/",n_cand_bins,"local reads record: ", tx_reads_local[0])
            #tx_cov = _child_get_permutation_fdr(tx_local, tx_reads_local, pval_cutoff=permutate_pval, min_cov=min_cov, max_iter=max_iter, seed=seed)
            res = do_permutation(transcr=tx_local, read_transcript=tx_reads_local, max_iter=max_iter, pval_cutoff=permutate_pval, min_cov=min_cov, seed=seed)
            tx_cov, _height_to_pval = res[0:2]
            tx_cov = max(tx_cov, min_cov)
            # print(datetime.datetime.now(),"\t\t",cand_bin_index,"/",n_cand_bins,"local res:",res)
            print(datetime.datetime.now(),"\t\t",cand_bin_index,"/",n_cand_bins,"local tx_cov:",tx_cov)
            # print(datetime.datetime.now(),"\t\t",cand_bin_index,"/",n_cand_bins,"local height_to_pval: ", _height_to_pval)
        # else:
        #     s = 0
        #     e = gene_length 
        if local_max-tx_cov <= 0:
            print(datetime.datetime.now(),"\t\t",cand_bin_index,"/",n_cand_bins,"local_max <= tx_cov, skip this bin !")
            cand_bin_index += 1
            continue
        print(datetime.datetime.now(),"\t\t","perm start pos: ", s)
        print(datetime.datetime.now(),"\t\t","perm end pos: ", e)
        
        # 4. extend boundary to local/global cov     
        # if local_max > min_cov:
        # find bounds when input signal drops below tx_cov
        # todo: re-center localmax pos in a plateau
        start = max_index   
        # start_trend = max_index 
        max_index_new = max_index
        while (start > left_bound) and (x[start - 1]-tx_cov >= decay*(local_max-tx_cov)): # while (start > left_bound) and (x[start - 1] >= decay*local_max)
            start -= 1
            if local_max < x[start]: 
                max_index_new = start
                local_max = x[start] # update local_max if x[start] may hihger than localmax, seem needed
                # if (x[start] > x[start+1]): 
                #     start_trend = start 
        # start = math.floor((start_trend + start) / 2)
        end = max_index
        while (end < (gene_length - 1)) and (x[end + 1]-tx_cov >= decay*(local_max-tx_cov)): # decay*local_max
            end += 1
            #local_max = max(local_max, x[end])
            if local_max < x[end]: 
                max_index_new = end
                local_max = x[end] # update local_max if x[end] may hihger than localmax, seem needed
        print(datetime.datetime.now(),"\t\t",cand_bin_index,"/",n_cand_bins,"start pos:",start)
        print(datetime.datetime.now(),"\t\t",cand_bin_index,"/",n_cand_bins,"end   pos:",end)
        
        # add current peak to results
        if ((end - start) < min_peak_length) or ((end - start) > max_peak_length):
            print(datetime.datetime.now(),"\t\t","finished (lenth not pass filter)",cand_bin_index,"/",n_cand_bins,"cand_bins !")
        else:
            #print(datetime.datetime.now(),(start, end, local_max))
            max_index_new,local_max = find_one_maxima_in_center_local(x[start:(end+1)]) # include end site
            max_index_new += start
            if (len(_height_to_pval) >= 1):
                    try:
                        tmp_fdr = _height_to_pval[int(local_max)]
                    except:
                        if int(local_max) >= max(_height_to_pval.keys()):
                            # local_max as key not found in _height_to_pval keys for some reason
                            tmp_fdr = float(0.0)
                        else:
                            # ?
                            tmp_fdr = float('nan')
            else:
                # _height_to_pval is null/blank dict
                tmp_fdr = float('nan')
            print(datetime.datetime.now(),"\t\t",cand_bin_index,"/",n_cand_bins,"tmp_fdr:",tmp_fdr)
            print(datetime.datetime.now(),"\t\t",cand_bin_index,"/",n_cand_bins,"cand_bins:", start, end, local_max, max_index_new, tx_cov, tmp_fdr)
            peaks.append([start, end, local_max, max_index_new, tx_cov, tmp_fdr]) # why not need deduplicate, adjacent peak not overlap: below left_bound/end updated, also next_cand_bin may jump to far away
            print(datetime.datetime.now(),"\t\t","finished",cand_bin_index,"/",n_cand_bins,"cand_bins !")

            

        # find next candidate bin
        left_bound = end
        next_cand_bin = end//half_bin_width
        while (cand_bin_index < n_cand_bins) and (cand_bins[cand_bin_index]+1 < next_cand_bin): # seem +1
            cand_bin_index += 1
        # else:
        #     print(datetime.datetime.now(),"\t\t","local_max <= min_cov, skip this bin !") # before local shuffle, more efficient
        cand_bin_index += 1
    print(datetime.datetime.now(),"\t","this tx all locals finished",cand_bin_index-1,"/",n_cand_bins,"!")
    return peaks



@command_handler
def call_peaks_localmax(args):
    import pyBigWig
    import pysam
    import numpy as np
    
    logger.info('\t' + 'read bam file: ' + args.input_bam)
    # logger.info('\t' + 'read bw file: ' + args.input_bw)
    logger.info("\t" + 'write bed file: ' + args.output_file)

    bam=pysam.AlignmentFile(args.input_bam, 'rb')
    # bigwig = pyBigWig.open(args.input_bw) # read pre-cal bw, or use clipper.readsToWiggle ?
    bed = open(args.output_file, 'w')
    # chroms = bigwig.chroms()
    chroms={x['SN']:x['LN'] for x in bam.header['SQ']}

    valid_chroms=set()
    for read in bam.fetch():
        valid_chroms.add(bam.get_reference_name(read.reference_id))
    print("\n\n",datetime.datetime.now(),"valid_chroms: ",len(valid_chroms),"/",len(chroms))
    chroms = {x:y for x,y in chroms.items() if x in valid_chroms}
    bam.reset()
    
    n_peaks = 0
    for chrom, size in chroms.items(): 
        # if chrom != 'ENST00000385227_____1': # ENST00000362283_____1, NR_146152_____1, 12965
        #     continue
        print("\n\n",datetime.datetime.now(),"start chrom: ",chrom,n_peaks,"/",len(chroms))
        print(datetime.datetime.now(),"size: ",size)
        
        this_tx_reads = bam.fetch(str(chrom), start=0, end=size) # , until_eof=True, BAM files with a large number of reference sequences are slow (https://pysam.readthedocs.io/en/latest/faq.html)
        tx_global = [chrom,0,size,"+",chrom] # 0-based
        this_tx_reads = read_tid_frag_from_bam(tid=tx_global, reads=this_tx_reads, is_unique=True) # convert 0-based bam iterable reads to 0-based reads records
        if len(this_tx_reads) < args.min_cov:
            n_peaks += 1
            print(datetime.datetime.now(),"skip (little tx reads) !")
            continue 
        count_signal = count_pileup_heights(tlen=size, reads=this_tx_reads)
        # count_signal = np.nan_to_num(bigwig.values(chrom, 0, size))
        print(datetime.datetime.now(),"count_signal: ", count_signal)
        print(datetime.datetime.now(),"sum(count_signal): ", sum(count_signal))
        print(datetime.datetime.now(),"len(this_tx_reads): ",len(this_tx_reads))

        peaks_chrom = _call_peaks_localmax(tx_id=chrom, count_signal=count_signal, tx_reads=this_tx_reads, 
            min_peak_length=args.min_peak_length, max_peak_length=200, bin_width=args.bin_width, min_cov=args.min_cov, decay=args.decay, local=True, permutate_pval=0.05, max_iter=50, seed=1234)
        for peak in peaks_chrom:
            n_peaks += 1
            #peaks.append([chrom, peak[0], peak[1], 'peak_%d'%n_peaks, peak[2], '+'])
            bed.write('%s\t%d\t%d\tpeak_%d\t%d\t+\t%d\t%d\t%.10f\n'%(chrom, peak[0], peak[1], n_peaks, peak[2], peak[3], peak[4], peak[5]))
    bam.close()
    # bigwig.close()
    bed.close()



if __name__ == '__main__':
    main_parser = argparse.ArgumentParser(description='Call peaks from exRNA signals')
    subparsers = main_parser.add_subparsers(dest='command')
    
    parser = subparsers.add_parser('call_peaks_localmax')
    parser.add_argument('--input-bam', type=str, required=True,
        help='input Bam file of raw reads')
    # parser.add_argument('--input-bw', type=str, required=True,
    #     help='input BigWig file of raw reads coverage')
    parser.add_argument('--min-peak-length', type=int, default=10,
        help='minimum length required for a peak')
    parser.add_argument('--decay', type=float, default=0.0,
        help='decay factor of peak summit to define peak boundary')
    parser.add_argument('--min-cov', type=float, default=5,
        help='minimum coverage required to define a peak')
    parser.add_argument('--bin-width', type=int, default=10,
        help='bin width to search enriched bins')
    parser.add_argument('--output-file', '-o', type=str, required=True,
        help='output peaks in BED format')


    args = main_parser.parse_args()
    if args.command is None:
        raise ValueError('empty command')
    logger = logging.getLogger('call_peak.' + args.command)

    command_handlers.get(args.command)(args)
    




# def chunkify(a, n):
#     """
# 	Separate a list (a) into consecutive n chunks.
# 	Returns the chunkified index
# 	"""
# 	k, m = len(a) / n, len(a) % n
# 	return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in xrange(n))

# def unpack_tid_to_qval(compact):
# 	"""
# 	Unpacks the returned values from multi-processing.
# 	"""
# 	tid_to_qval=defaultdict(list)
# 	# for item in compact:
# 	# 	unique, combined = item[0], item[1]
#     for tid in compact:
#         if len(compact[tid])>0:
#             tid_to_qval[tid]=compact[tid]
# 	return tid_to_qval

# def read_gtf(fn):
#     """read in the gene annotation from GTF file
# 	"""
# 	logger.info('read GTF from "%s" '% fn)
# 	gene_annot = {}
# 	with open(fn, 'r') as f:
# 		for line in f:
# 			if line.startswith('#'):
# 				continue
# 			ele = line.strip().split('\t')
# 			if ele[2] != 'gene':
# 				continue
# 			chr, start, end, strand = ele[0], int(ele[3]), int(ele[4]), ele[6]
# 			try:
# 				gene_id = re.search(r'gene_id "(.+?)"', ele[-1]).group(1)
# 			except AttributeError:
# 				continue
# 			gene_annot[gene_id] = [chr, start, end, strand, gene_id]
# 	return gene_annot

# def _child_get_permutation_fdr(tx_local, reads, pval_cutoff, min_cov, max_iter, seed):
#     """
#     original: for all tx bam, adapted to tx specific reads
#     General permutation wrapper for a list of genes. Gets called by multi-processing generated by Pool()
#     Returns packed FDRs from each child process.
#     is_stranded=False,
#     """
#     random.seed(seed)
#     # tid_to_qval=defaultdict(list)	
#     # pid = os.getpid()
#     # tot = len(child_gene_list)

#     # for i in range(len(child_gene_list)):
#     # 	if not i % 200:
#     # 		logger.debug('pid %s : %i / %i (%.2f%%)'% (pid, i, tot, float(i)/float(tot)*100))
#         # gene_name = child_gene_list #[i]
#         # gene = gene_annot[gene_name]
#     # chr, start, end, strand, tid = tx_local[0:5]
#     # reads = read_tid_frag_from_bam(tx_local, reads, is_unique=True)
#     h0 = do_permutation(tx_local, reads, max_iter, pval_cutoff=pval_cutoff, min_cov=min_cov, seed = seed)
#     return h0

#clipper
# def get_full_length_cigar(read):
#     positions = enumerate(read.positions)
#     for t in read.cigartuples:
#         value, times = t

#         #value 3 is splice junction value 2 is deletion in read
#         if value == 1 or value == 3 or value == 4:
#             continue

#         if value == 2:  #In the case of a deletion yeild from the previous current position, this could be buggy in the case of deletions occuring after anything but a match
#             for x in xrange(times):
#                 cur_position += 1
#                 yield cur_position

#         for x in xrange(times):
#             cur_position = next(positions)[1]
#             yield cur_position



#WPS
# def isSoftClipped(cigar):
#       #Op BAM Description
#   #M 0 alignment match (can be a sequence match or mismatch)
#   #I 1 insertion to the reference
#   #D 2 deletion from the reference
#   #N 3 skipped region from the reference
#   #S 4 soft clipping (clipped sequences present in SEQ)
#   #H 5 hard clipping (clipped sequences NOT present in SEQ)
#   #P 6 padding (silent deletion from padded reference)  ???
#   #= 7 sequence match
#   #X 8 sequence mismatch
#   for (op,count) in cigar:
#     if op in [4,5,6]: return True
#   return False

# def aln_length(cigarlist):
#   tlength = 0
#   for operation,length in cigarlist:
#     if operation == 0 or operation == 2 or operation == 3 or operation >= 6: tlength += length
#   return tlength

# lseq = aln_length(read.cigar)


# def get_FDR_cutoff_mean(readlengths, genelength,iterations=100, mincut=5,alpha=0.05):
#     """
#     Returns an int, the number of reads needed to meet the FDR cutoff by randomized method
#     TODO: Allow the minimum cutoff to be paramaritizied
#     TODO: double check math on this
#     :param readlengths: list,  list of lengths of aligned portions of reads
#     :param genelength: int, effective gene length (unalignable regions aren't counted)
#     :param iterations: int, default 100
#     :param mincut: int, default 2, min threshold possible to return
#     :param alpha: float, default 0.05, FDR alpha
#     :return: int, min number of reads per position to read FDR
#     """

#     # if you have very few reads on a gene, don't waste time
#     # trying to find a cutoff
#     if len(readlengths) < 20:
#         return mincut

#     results = shuffle(int(genelength), int(iterations), 0, .05, readlengths) # alpha ?
#     total = 0

#     # parses results from peaks script, calculates mean from peaks results
#     # should document peaks function call return value somewhere around here

#     for cut, n_observed in enumerate(results):
#         total += (cut * n_observed)

#     # logic for min cutoffs
#     cutoff = total / iterations
#     if cutoff < mincut:
#         cutoff = mincut
#     return int(round(cutoff, 0))