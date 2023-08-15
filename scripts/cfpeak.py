#!/usr/bin/env python 

import argparse, sys, os, errno
from distutils.spawn import spawn
import logging
from tkinter import N
import numpy as np
import datetime
import pandas as pd
from pandas import DataFrame, Series
from scipy.fftpack import fft
from scipy.signal import convolve
from scipy import stats
import pyBigWig
import pysam
from collections import defaultdict, OrderedDict
import math, random, re, time
from multiprocessing import Pool
from threading import Thread
# import numba
# import matplotlib
# logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')


command_handlers = {}
def command_handler(f):
    command_handlers[f.__name__] = f
    return f



def count_pileup_heights(tlen, reads, downsample=True):
    """
    Partially adapted from CLAM
    Used within do_permutation()
    Counts the distribution of coverage/heights for a given gene/permutation
    Test:
        reads=[['r1', 0, 20, 1], ['r2', 3, 21, 1], ['r2', 4, 21, 1]]
        tlen=80
    """
    loc_heights = np.zeros(tlen)
    if not downsample:
        for id, pos, read_len, score in reads:
            loc_heights[max(0, pos):min(pos+read_len, tlen)] += score
        return loc_heights

    else:
        # count read lengths
        # read_lengths = {}
        # for read in reads:
        #     read_len = read[2]
        #     read_lengths[read_len] = read_lengths.get(read_len, 0) + 1
        # sum_RL = sum(np.array(list(read_lengths.keys())) * np.array(list(read_lengths.values())))
        sum_RL = sum(read[2] for read in reads) #sum_RL = sum(read[2] for read in reads)
        ratio = (1000 * tlen) / sum_RL
        if  ratio<1:  
            reads_dwnsamp = []
            for read in reads:
                if np.random.random() <= ratio:
                    reads_dwnsamp.append(read)
        else:
            reads_dwnsamp = reads

        for id, pos, read_len, score in reads_dwnsamp:
            for i in range(max(0, pos), min(pos+read_len, tlen)):
                loc_heights[i] += score

        if ratio<1: 
            loc_heights = [int(cov/ratio) for cov in loc_heights]

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
def permutate_heights(tlen, reads, downsample=True):
    """
    Partially adapted from CLAM
    Used within do_permutation()
    Randomly allocate the read locations.
    Test:
        reads=[['r1', 0, 20, 1], ['r2', 3, 21, 1], ['r2', 4, 21, 1]]
        tlen=80
    """
    loc_heights = np.zeros(tlen)
    
    if not downsample:
        reads_dwnsamp = reads

        for id, pos, read_len, score in reads_dwnsamp:
            rand_pos=np.random.randint(-read_len, max(1, tlen)) # consider read overlapped with range
            loc_heights[max(0, rand_pos):min(rand_pos + read_len, tlen)] += 1
            # for i in range(max(0,rand_pos), min(rand_pos + read_len, tlen)): # consider read overlapped with range    
            #     loc_heights[i]+=1
                
    else:
        # count read lengths
        # read_lengths = {}
        # for read in reads:
        #     read_len = read[2]
        #     read_lengths[read_len] = read_lengths.get(read_len, 0) + 1
        # sum_RL = sum(np.array(list(read_lengths.keys())) * np.array(list(read_lengths.values())))
        sum_RL = sum(read[2] for read in reads)
        ratio = (1000 * tlen) / sum_RL
        
        if ratio < 1: # dwn samp reads
            reads_dwnsamp = []
            for read in reads:
                if np.random.random() <= ratio:
                    reads_dwnsamp.append(read)
        else: 
            reads_dwnsamp = reads

        for id, pos, read_len, score in reads_dwnsamp:
            rand_pos=np.random.randint(-read_len, max(1, tlen)) # consider read overlapped with range
            loc_heights[max(0, rand_pos):min(rand_pos + read_len, tlen)] += 1

                
        if ratio < 1: # dwn samp reads
            loc_heights = [int(cov/ratio) for cov in loc_heights]
        
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
    Partially adapted from CLAM
    Permutes the reads along a given gene length
    Returns the locally corrected p-values for each observed height on the given gene.
    """
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

    for i in obs_heights_count:
        if i==0:
            continue
        else:
            rand_heights_dist[int(i)]+=1
            rand_sum+=1
    for B in range(max_iter):
        new_heights_count=permutate_heights(tlen=tid_length, reads=read_transcript, downsample=True)
        for i in new_heights_count:
            if i==0:
                continue
            else:
                rand_heights_dist[int(i)]+=1
                rand_sum+=1
    height_to_pval={}
    for h in set(rand_heights_dist.keys()):
        if h < 1:
            continue 
        else:
            lefter=0
            for j in range(int(h), max(rand_heights_dist)+1): 
                lefter+=rand_heights_dist[j]
            height_to_pval[int(h)]=lefter/float(rand_sum) 

    not_sig = {key:val for key,val in height_to_pval.items() if val > pval_cutoff} 
    try:
        h0_tmp = min(not_sig, key=lambda k: not_sig[k]) 
    except:
        h0_tmp = int(np.percentile(list(rand_heights_dist.keys()), (95)))
        print(datetime.datetime.now(),"do_permutation exception !")
    h0_tmp = int(max(h0_tmp,min_cov))

    res = []
    res.append(h0_tmp)
    res.append(height_to_pval)
    
    return res


def find_local_maxima(arr):
    """
    Partially adapted from CLIPper
    Returns a list of values for an array that mark if a value is a local
    maxima (pos with value > 0)
    Test:
        a = np.array([1,2,3,2,1,2,1]) #[-1., -1.,  3., -1., -1.,  2., -1.]
        b = np.array([1,2,3,3,2,2,2]) #[-1., -1.,  3., -1., -1., -1., -1.]
        c = np.array([3,2,1,1,1,2,2]) #[ 3., -1., -1., -1., -1., -1.,  2.]
        c2 = np.array([3,3,1,1,1,2,2]) #[ 3., -1., -1., -1., -1., -1.,  2.]
        c3 = np.array([3,3,3,1,2,2,2]) #[-1.,  3., -1., -1., -1.,  2., -1.]
        c4 = np.array([3,3,3,3,1,2,2,2,2]) #[-1.,  3., -1., -1., -1., -1., -1.,  2., -1.]
        d = np.array([2,1,1,1,1,2,2]) #[ 2., -1., -1., -1., -1., -1.,  2.]
        e = np.array([1,2,3,4,5,6,7]) #[-1., -1., -1., -1., -1., -1.,  7.]
        f = np.array([7,6,5,4,3,2,1]) #[ 7., -1., -1., -1., -1., -1., -1.]
        g = np.array([2,2,2,2,2,2,2]) #[-1., -1., -1.,  2., -1., -1., -1.]
    """
    if len(arr) <= 2:
        print("require arr length >2")
        return -1
    # to initalize a new array to all false
    maxima = np.empty(len(arr)) # , dtype='bool'
    maxima.fill(-1)
    # print("array:",arr)
    max_range_start = 0
    increasing = True
    for i in range(len(arr[:-1])): # -1: bellow i+1

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

    return maxima 

def find_local_minima(arr,include_boundary=False):
    """
    Partially adapted from CLIPper
    Returns a list of values for an array that mark if a value is a local
    minima (pos with value > 0; default not include minima at boundary)
    Test:
        a = np.array([1,2,3,2,1,2,1]) #[-1., -1., -1., -1.,  1., -1., -1.]
        b = np.array([1,2,3,3,2,2,2]) #[-1., -1., -1., -1., -1., -1., -1.]
        c = np.array([3,2,1,1,1,2,2]) #[-1., -1., -1.,  1., -1., -1., -1.]
        d = np.array([2,1,1,1,1,2,2]) #[-1., -1.,  1., -1., -1., -1., -1.]
        d1 = np.array([2,1,1,1,1,1,2]) #[-1., -1., -1.,  1., -1., -1., -1.]
        e = np.array([1,2,3,4,5,6,7]) #[-1., -1., -1., -1., -1., -1., -1.]
        f = np.array([7,6,5,4,3,2,1]) #[-1., -1., -1., -1., -1., -1., -1.]
        g = np.array([2,2,2,2,2,2,2]) #[-1., -1., -1., -1., -1., -1., -1.]
    """
    if len(arr) <= 2:
        print("require arr length >2")
        return -1
    
    minima = np.empty(len(arr)) #(arr == -1) # initialize at all false
    minima.fill(-1) # diff from maxima, cause minima cov may be 0 
    min_range_start = 0
    if include_boundary:
        decreasing = True # include boundary site 
    else:
        decreasing = False
    
    for i in range(len(arr[:-1])): # -1: bellow i+1
        # array needs to be smooth for this to work, otherwise we'll
        # run into odd edge cases
        # update location of minima start until
        if arr[i] > arr[i + 1]:
            min_range_start = i + 1
            decreasing = True
            if include_boundary:
                if min_range_start == len(arr)-1: #include boundary site
                    minima[int(min_range_start)] = int(arr[min_range_start]) #include boundary site
        if (arr[i] < arr[i + 1]) and decreasing is True:
            decreasing = False
            # gets the local minima midpoint
            midpoint = int((min_range_start + i) / 2)
            minima[int(midpoint)] = int(arr[midpoint]) #True
    if np.sum(np.array(minima))==len(arr)*(-1):
        # minima=[]
        print("null minima created for flat cov, this happen not offen, return -1 list !")
        # midpoint = int((0 + len(arr)) / 2)
        # minima[int(midpoint)] = arr[int(midpoint)]
    # print("minima: ", minima)
    return minima


# def find_one_maximum_in_center_local(arr):
#     """
#     Returns only one bigest pos and cov (deprecated)
#     """
#     maxima = find_local_maxima(arr)
    
#     # find one most max cov pos in center
#     indexed_maxima = {i:val for i,val in enumerate(maxima) if val!=-1}
#     maxima_pos = max(indexed_maxima, key=lambda k: indexed_maxima[k]) # return key of max value (first key if multi max exist)
#     maxima_cov = indexed_maxima[maxima_pos]
#     return maxima_pos, int(maxima_cov)

def findMiddle(input_list):
    middle = float(len(input_list))/2
    if middle % 2 != 0:
        return int(middle - .5) # left only
    else:
        return int(middle-1) # left; right: int(middle)
    
def find_one_maximum_in_center_local(arr):
    """
    Returns only one bigest pos and cov
    Test:
        arr2 = [-1, -1, 2, -1, -1] #(2, 2)
        arr2 = [-1, -1, 2, -1, 2] #(2, 2)
        arr2 = [-1, 2, -1, 2, -1] #(1,2)
        arr2 = [2, -1, -1, -1, 2] #(0, 2)
        arr2 = [2, -1, 2, -1, 2] #(2, 2)
        arr2 = [2, -1, 2, -1, 2, -1, -1] #(2, 2) 
        arr2 = [-1, -1, 2, -1, 2, -1, 2] #(2, 2) 
        arr2 = [-1, -1, -1, 2, -1, 2, -1, 2] #(3, 2)
        arr2 = [-1, -1, -1, -1, 2, -1, 2, -1, 2] #(4, 2)
        arr2 = [2, -1, -1, -1, 3] #(4, 3)
        arr2 = [-1, -1, -1, -1, -1] #(2, -1)
    """
    arr2 = find_local_maxima(arr)
    middle_index = len(arr2) / 2 - 0.5 # 0-based idx for center pos, eg: 0.5,1)
    maxDep = max(arr2)
    indexed_maxima = [[i,val] for i, val in enumerate(arr2) if val==maxDep]
    # maxima_rank = len(indexed_maxima) // 2
    
    if len(indexed_maxima) == 1:
        max_pos_val = indexed_maxima[0]
    else:
        maxima_index = [i for i, val in enumerate(arr2) if val==maxDep]
        maxima_index_to_mid = [*map(lambda x: x - middle_index, maxima_index)]
        maxima_index_to_mid = [abs(i) for i in maxima_index_to_mid]
        maxima_index_to_mid_idx = [i for i, val in enumerate(maxima_index_to_mid) if val==min(maxima_index_to_mid)] # len<=2
        maxima_index_to_mid_idx = maxima_index_to_mid_idx[0] # random choose one if two sym maxima exist
        idx = int(maxima_index[maxima_index_to_mid_idx]) # need int to rm float
        max_pos_val = [idx,int(arr2[idx])]
    return max_pos_val[0],max_pos_val[1]

# def find_one_minimum_in_center_local(arr):
#     """
#     Returns only one smallest pos and cov (deprecated)
#     """
#     minima = find_local_minima(arr)
#     try:
#         # find one most min cov pos in center
#         indexed_minima = {i:val for i,val in enumerate(minima) if val!=-1} # mini cov may be zero
#         minima_pos = min(indexed_minima, key=lambda k: indexed_minima[k]) # return key of min value (first key if multi min exist)
#         minima_cov = indexed_minima[minima_pos]
#         return minima_pos, int(minima_cov)
#     except:
#         print("null returned in find_one_minimum_in_center_local !")
#         return -1,-1
    
def find_one_minimum_in_center_local(arr2):
    """
    Returns only one minium pos and cov
    Test:
        arr2 = [-1, -1, 2, -1, -1] #(2, 2)
        arr2 = [-1, -1, 2, -1, 2] #(2, 2)
        arr2 = [-1, 2, -1, 2, -1] #(1,2)
        arr2 = [2, -1, -1, -1, 2] #(0, 2)
        arr2 = [2, -1, 2, -1, 2] #(2, 2)
        arr2 = [2, -1, 2, -1, 2, -1, -1] #(2, 2) 
        arr2 = [-1, -1, 2, -1, 2, -1, 2] #(2, 2) 
        arr2 = [-1, -1, -1, 2, -1, 2, -1, 2] #(3, 2)
        arr2 = [-1, -1, -1, -1, 2, -1, 2, -1, 2] #(4, 2)
        arr2 = [2, -1, -1, -1, 3] #(0, 2) !!!
        arr2 = [-1, -1, -1, -1, -1] #(-1, -1) !!!
    """
    # arr2 = find_local_minima(arr)
    middle_index = len(arr2) / 2 - 0.5 # 0-based idx for center pos, eg: 0.5,1)
    if max(arr2) == -1:
        return -1, -1
    else:
        minDep = min(list(filter(lambda a: a != -1, arr2))) # (-1 exist)
        indexed_minima = [[i,val] for i, val in enumerate(arr2) if val==minDep] 
        if len(indexed_minima) == 1:
            min_pos_val = indexed_minima[0]
        else:
            minima_index = [i for i, val in enumerate(arr2) if val==minDep] 
            minima_index_to_mid = [*map(lambda x: x - middle_index, minima_index)]
            minima_index_to_mid = [abs(i) for i in minima_index_to_mid]
            minima_index_to_mid_idx = [i for i, val in enumerate(minima_index_to_mid) if val==min(minima_index_to_mid)] # len<=2
            minima_index_to_mid_idx = minima_index_to_mid_idx[0] # random choose one if two sym minima exist
            idx = int(minima_index[minima_index_to_mid_idx]) # need int to rm float
            min_pos_val = [idx,int(arr2[idx])]
        return min_pos_val[0],min_pos_val[1]


def merge_maximas_by_valley_cov(maxs_list,mins_list,decay):
    """
    Input:
        decay=0.5
        maxs_list = [[0,6],[2,12],[5,10],[7,9]]
        mins_list = [   [1,2], [3,5], [6,6]   ]
    return:     [[0,6],[2,12],[5,10]      ]
    inspired by NAR,CTK,valley-seeking
    """
    # last_num = len(maxs_list)
    # if last_num <=1:
    #     print("maxs_list size <=1, not iterable, not changed")
    #     return maxs_list
    # else:  
    #     maxs_list2 = maxs_list[:-1] # init, enable 1st iter
    #     while len(maxs_list2) != last_num:
    # sanity check
    maxs_list = sorted(maxs_list)
    mins_list = sorted(mins_list)
    # if len(maxs_list)!=len(mins_list)+1:
    #     print("wrong max and min num !")
    maxs_list2 = [] #maxs_list.copy() # copy as backup
    maxs_set_rm_idx = set() #[]
    keep_second_max = False
    n_iter = 0
    
    for i in range(len(maxs_list)-1):
        # print("i: ",i)
        pos_1,cov_1 = maxs_list[i]
        pos_2,cov_2 = maxs_list[i+1]
        if (mins_list[i][1] <= min(cov_1,cov_2)*decay) or (abs(pos_1-pos_2)>500): #  
            if n_iter==0:
                maxs_list2.append(maxs_list[i])
                # keep_second_max = True
            maxs_list2.append(maxs_list[i+1])
            keep_second_max = True
        elif (cov_1 < cov_2) and ((mins_list[i][1] > cov_1*decay) or (abs(pos_1-pos_2)<=10)):
            if keep_second_max:
                # del maxs_list2[i-1] # ?
                # maxs_list2[i]=[-1,-1] # mark to rm
                maxs_set_rm_idx.add(i) # mark to rm
            maxs_list2.append(maxs_list[i+1])
            keep_second_max = True
        elif (cov_1 >= cov_2) and ((mins_list[i][1] > cov_2*decay) or (abs(pos_1-pos_2)<=10)):
            if n_iter==0:
                maxs_list2.append(maxs_list[i])
            # if not keep_second_max:
            #     maxs_list2.append(maxs_list[i])
            keep_second_max = False
            # continue #maxs_list2.append(maxs_list[i])
        n_iter += 1 
    maxs_list_rm_list = [maxs_list[i] for i in maxs_set_rm_idx]
    maxs_list2 = [i for i in maxs_list2 if i not in maxs_list_rm_list]
    # maxs_list2 = [i for i in maxs_list2 if i[0]!=-1] 
    maxs_list2 = sorted(maxs_list2)
    # uniq_maxs_list2 = list(set(maxs_list2)) # TypeError: unhashable type: 'list'
    # [uniq_maxs_list2.append(i) for i in maxs_list2 if i not in uniq_maxs_list2] # rm potential dup records
    # maxs_list2 = sorted(uniq_maxs_list2)
    # last_num = len(maxs_list2)
    return maxs_list2

def find_one_minima_between_maxima(arr,maxima_list):
    """
    get minima (list) from each pair of maximas in maxima_list
    minima not include maxima site
    """
    maxima_list = sorted(maxima_list)
    # return null list if only one maxima is present
    if len(maxima_list)<=1:
        print("input maxima_list num <= 1 in find_one_minima_between_maxima !")
        return [[]]
    else: 
        minima_list2 = []
        for i in range(len(maxima_list)-1):
            pos_1,cov_1 = maxima_list[i] 
            pos_2,cov_2 = maxima_list[i+1]
            minima_list2_pos,minima_list2_cov = find_one_minimum_in_center_local(arr[pos_1:(pos_2+1)]) # need +1 !
            minima_list2_pos += pos_1# convert coordinate            
            if minima_list2_pos!=-1: # -1: no mini for one maxi
                minima_list2.append([minima_list2_pos,minima_list2_cov])
            else:
                minima_list2.append([])
        minima_list2 = sorted(minima_list2)
        return minima_list2

def find_one_minima_from_multi_minima(maxima_list, mins_list):
    """
    avoid re-cal minimas, filter from original full range minimas
    """
    maxima_list = sorted(maxima_list)
    # return null list if only one maxima is present
    if len(maxima_list)<=1:
        print("input maxima_list num <= 1 in find_one_minima_from_multi_minima !")
        return [[]]
    else: 
        minima_list2 = []
        for i in range(len(maxima_list)-1):
            pos_1,cov_1 = maxima_list[i] 
            pos_2,cov_2 = maxima_list[i+1]
            filter_minimas = {i[0]:i[1] for i in mins_list if i[0]>pos_1 and i[0]<pos_2}
            minima_list2_pos = min(filter_minimas, key=lambda k: filter_minimas[k]) # return key with mini val, only one (first if multi exist)
            minima_list2_cov = filter_minimas[minima_list2_pos]
            # this func. not need convert coordinate            
            minima_list2.append([minima_list2_pos,minima_list2_cov])
        # minima_list2 = sorted(minima_list2)
        return minima_list2




def read_tid_frag_from_bam(tid, reads, is_unique):
    """
    Partially adapted from CLAM
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
    # reads=[x for x in reads if x.pos>=int(start) and x.positions[-1]<=int(end)] # x.pos, x.positions: 0-based; tid: 1-based 
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

def read_tid_frag_from_readList(tid, reads, is_unique):
    """
    convert to relative tx pos
    """
    tid_reads=[]
    #gene, chr, strand, start, end=tid
    chr, start, end, strand, gene = tid[0:5]
    for read in reads:
        r_chr,r_qname,r_start,r_len,r_score=read
        if ( int(r_start)>=int(start) and int(r_start+r_len)<=int(end) and r_chr==chr):    #to aviod long read:  and read_length<200; to avoid junction reads: (not 'N' in read.cigarstring), remove these to keep in line with bw signal 
            tid_reads.append([r_qname, r_start-start, r_len, r_score])
    return tid_reads # [["a",1,2],["b",2,3]...]

def read_tid_frag_from_records(tid, reads, full_reads_within=False):
    """
    Use reads records [read.qname, read.pos-start, read_length, score] info for a given gene and its loci.
    Returns reads, read weights and its mapped loci.
    is_stranded, is_unique
    """
    tid_reads=[]
    #gene, chr, strand, start, end=tid
    chr, start, end, strand, gene = tid[0:5]
    if not full_reads_within:
        reads=[x for x in reads if (x[1]+x[2])>=int(start) and x[1]<=int(end)] # include all overlapped reads, not just full reads inside, to aviod tx_cov/localmax not included in _height_to_pval, should result to fewer nan pvalue
    else: 
        reads=[x for x in reads if (x[1]+x[2])>=int(start) and x[1]<=int(end)] # only include reads that full range within tx, useful for poisson test ? 
    for read in reads:
        tid_reads.append([read[0], read[1]-start, read[2], read[3]])
    return tid_reads # [["a",1,2],["b",2,3]...]

def get_readsList_from_bam(gene_list, bam, is_unique=True):
    # reads=bam.fetch() # all aligned reads
    # bam.reset()
    all_tx_reads_tmp = []
    for gene_list_tmp in gene_list:
        # _size = chroms[gene_list_tmp]
        # _tx_global = [gene_list_tmp,0,_size,"+",gene_list_tmp] # 0-based
        these_tx_reads_tmp = bam.fetch(contig=str(gene_list_tmp)) # , start=0, end=_size
        # _all_tx_reads_tmp = read_tid_frag_from_bam(tid=_tx_global, reads=_this_tx_reads, is_unique=True)
        
        # chr, start, end, strand, gene = _tx_global[0:5]
        # reads=[x for x in reads if x.pos>=int(start) and x.positions[-1]<=int(end)] # x.pos, x.positions: 0-based; tid: 1-based 
        for read in these_tx_reads_tmp:
            # if is_unique:
            score=1
            # try:
            #     read_length = read.opt('RL')
            # except:
            read_length = read.positions[-1] - read.positions[0] + 1 # new pysam: get_reference_positions
            
            # if (read.positions[0]-start>=0 and end-read.positions[-1]>=0):    #to aviod long read:  and read_length<200; to avoid junction reads: (not 'N' in read.cigarstring), remove these to keep in line with bw signal 
            all_tx_reads_tmp.append([bam.get_reference_name(read.reference_id), read.qname, read.pos, read_length, score])                
        # all_tx_reads_tmp.append(_this_tx_reads)
    return all_tx_reads_tmp


def chunkify(a, n):
    """
    Separate a list (a) into consecutive n chunks.
    Returns the chunkified index
    """
    k, m = int(len(a) / n), len(a) % n # py2: no need for int()
    chunkify_tmp = [a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n)] # original: () generator
    return chunkify_tmp # xrange in py2

def get_chunkify_readsList(child_gene_list, all_tx_reads):
    chunkify_gene_reads = []
    n = len(child_gene_list)
    for i in range(n):
        chunkify_gene_reads.append([])
    for read in all_tx_reads:
        ref_id,qname,pos,rlen,score = read
        for i in range(n):
            if ref_id in child_gene_list[i]:
                chunkify_gene_reads[i].append(read) # not add [read]
                break # skip this if
    return chunkify_gene_reads

def poissonP(reads_in_gene, reads_in_peak, gene_length, peak_length):
    """
    Partially adapted from CLIPper
    scipy.stats.poisson.cdf
    compute the p-value for a peak of peak_length length with reads_in_peak reads,
    given that the read is gene_length long and has reads_in_gene reads

    If there are fewer than 3 reads expected to fall in the region, assume there's 3 reads
    expected...

    Paramaters
    ----------
    reads_in_gene: Integer representing number of reads in gene
    reads_in_peak: Integer representing the number of reads in a specific peak
    gene_length: Integer representing length of gene
    peak_length: Integer representing length of peak
    Returns double, the p-value that the peak is significant
    If calculation fails returns 1
    """
    try:
        # lam is estimate of the lambda value
        # poission takes a value and the lambda
        # this is average number of reads per single
        # site in the gene, but
        # a peak is not a single site, so it the average number
        # gets multipled by the peak
        # length as an estimator of the mean
        lam = 1 + ((float(reads_in_gene) / float(gene_length)) * float(peak_length))  # expect at least one read.
        cum_p = stats.poisson.sf(reads_in_peak, lam)
        return cum_p
    except Exception as error:
        # logging.error("Poisson cutoff failled %s " % (error))
        print("Poisson failed")
        return 0


# merge the intervals
def mergeIntervals(intervals):
    intervals.sort(key=lambda x: x[0])
    merged = []
    for interval in intervals:
        if not merged or merged[-1][1] < interval[0]:
            merged.append(interval)
        else:
            merged[-1][1] = max(merged[-1][1], interval[1])
    return merged

def gini(x):
    total = 0
    for i, xi in enumerate(x[:-1], 1):
        total += np.sum(np.abs(xi - x[i:]))
    return total / (len(x)**2 * np.mean(x))

def _call_peaks_localmax(tx_id, count_signal, tx_reads, min_peak_length, max_peak_length, bin_width, min_cov, permutate_pval, poisson_pval, decay, mode, max_iter, seed, boundary,recursive_maxima_elimination):
    '''
    Input: bw signals,
    Return: local or global maxima position
    '''
    np.random.seed(seed)

    x = np.array(count_signal)
    gene_length = x.shape[0] # gene_length = len(x)

    tx_global = [tx_id,0,gene_length,"+",tx_id]
    tx_global_reads = read_tid_frag_from_records(tx_global, tx_reads) # read_tid_frag_from_records(tx_global, tx_reads) # get converted tx reads 
    print(datetime.datetime.now(),"\t","full tx reads num: ", len(tx_global_reads))
    
    res = do_permutation(transcr=tx_global, read_transcript=tx_global_reads, max_iter=max_iter, pval_cutoff=permutate_pval, min_cov=min_cov, seed=seed)
    tx_cov, _height_to_pval = res[0:2]
    tx_cov = max(tx_cov, min_cov)
    print(datetime.datetime.now(),"\t","full tx coverage: ", tx_cov)
    
    # 0. average signal over bins with 50% overlap, and filter low cov bins (<min_cov)
    half_bin_width = bin_width//2
    n_bins = max(1, gene_length//half_bin_width)
    bin_cov = np.zeros(n_bins)
    for i in range(n_bins):
        bin_cov[i] = np.mean(x[(i*half_bin_width):min((i + 2)*half_bin_width, gene_length)]) # not gene_length-1
    
    cand_bins = np.nonzero(bin_cov > min_cov)[0] # return 0-based half-bin index in bin_cov
    # n_cand_bins = cand_bins.shape[0]

    cand_bins_list = [[i*half_bin_width,min((i + 2)*half_bin_width, gene_length)] for i in cand_bins] # not gene_length-1

    cand_bins_merge = mergeIntervals(cand_bins_list) # get merged intervals
    print("cand_bins_merge: ",cand_bins_merge)
    print("decay: ",decay)
    
    # find all local maximas in each merged candidate bins/regions
    # find all local minimas between each adjacent found maximas
    all_maximas = []
    all_minimas = []
    for cand_bin_merge_index in range(len(cand_bins_merge)):
        this_bin_merge_start,this_bin_merge_end = cand_bins_merge[cand_bin_merge_index]
        print("\n","this_bin_merge_start: ",this_bin_merge_start)
        print("this_bin_merge_end: ",this_bin_merge_end)
        this_maxima = find_local_maxima(arr=x[this_bin_merge_start:this_bin_merge_end]) # converted pos
        this_indexed_maxima = sorted([[i,int(val)] for i,val in enumerate(this_maxima) if val>0]) # rm other pos, sort by max pos then cov, small->large
        print("this_indexed_maxima: ",this_indexed_maxima)
        this_indexed_maxima_filter = [i for i in this_indexed_maxima if i[1]>min_cov] # filter low cov
        
        print("this_indexed_maxima_filter: ",this_indexed_maxima_filter)
        this_indexed_minima_filter = find_one_minima_between_maxima(arr=x[this_bin_merge_start:this_bin_merge_end], maxima_list=this_indexed_maxima_filter)
        # this_indexed_minima_filter = sorted(this_indexed_minima_filter)
        print("this_indexed_minima_filter: ",this_indexed_minima_filter)
        
        # iteration
        last_num = len(this_indexed_maxima_filter)
        this_indexed_maxima_filter_merge = this_indexed_maxima_filter # .copy()
        this_indexed_minima_filter_merge = this_indexed_minima_filter # .copy()
        
        if last_num <=1:
            print("maxs_list size <=1, not iterable, not changed")

            # return maxs_list
        elif not recursive_maxima_elimination:
            print("recursive_maxima_elimination == False, not changed")
        else:  
            last_num = 0 # init, enable 1st iter
            iter_num = 0

            # init_num = len(this_indexed_maxima_filter_merge)
            while (len(this_indexed_maxima_filter_merge) != last_num) and (len(this_indexed_maxima_filter_merge) > 1): 
                print("last num: ",last_num)
                print("iter num: ",iter_num)
                last_num = len(this_indexed_maxima_filter_merge)
                this_indexed_maxima_filter_merge = merge_maximas_by_valley_cov(maxs_list=this_indexed_maxima_filter_merge,mins_list=this_indexed_minima_filter_merge,decay=decay) # iterate filter adjacent minor peaks to one major peak
                # this_indexed_minima_filter_merge = find_one_minima_between_maxima(arr=x[this_bin_merge_start:this_bin_merge_end], maxima_list=this_indexed_maxima_filter_merge)
                this_indexed_minima_filter_merge = find_one_minima_from_multi_minima(maxima_list=this_indexed_maxima_filter_merge, mins_list=this_indexed_minima_filter) # avoid re-cal minimas, filter from original full range minimas
                iter_num += 1
                print("this_indexed_maxima_filter_merge: ",this_indexed_maxima_filter_merge)
                print("this_indexed_minima_filter_merge: ",this_indexed_minima_filter_merge)
        # convert to full-len tx coordinate before append to all 
        this_indexed_maxima_filter_merge = [[i[0]+this_bin_merge_start,i[1]]for i in this_indexed_maxima_filter_merge]
        try: #
            this_indexed_minima_filter_merge = [[i[0]+this_bin_merge_start,i[1]]for i in this_indexed_minima_filter_merge] # this_indexed_minima_filter_merge may be null [[]]
        except: 
            this_indexed_minima_filter_merge = [[]]
            
        all_maximas.append(this_indexed_maxima_filter_merge) #find_local_maxima() # find_one_maximum_in_center_local
        all_minimas.append(this_indexed_minima_filter_merge)
    print("all_maximas: ",all_maximas)
    print("all_minimas: ",all_minimas)
        


    all_maximas_flat = [element for sublist in all_maximas for element in sublist]
    n_cand_peaks = len(all_maximas_flat)
    n_cand_peak = 0
    # cand_bin_index = 0
    # left_bound = 0
    locals_list = [[0,0]] # avoid re-calculate local region above min_cov, init with 0,0
    fold_list = [] # expanded local is 2 fold length of original peak length (equally extend from peak center to both side)
    peaks = []    
    for cand_bin_merge_index in range(len(cand_bins_merge)):
        this_bin_merge_start,this_bin_merge_end = cand_bins_merge[cand_bin_merge_index]
        tmp_maximas = all_maximas[cand_bin_merge_index]
        tmp_minimas = all_minimas[cand_bin_merge_index]

        print("tmp_maximas: ",tmp_maximas)
        print("tmp_minimas: ",tmp_minimas)
        

        left_bound = this_bin_merge_start
        left_valley = min_cov 
        for peak_index in range(len(tmp_maximas)): # cand_bin_index: localmax index in a merged bin
            # 1. find local max pos and value (use fixed hard threshold cov: min_cov, eg. 5)
            try:
                right_bound,right_valley = tmp_minimas[peak_index] # may be null list: []
            except:
                right_bound = this_bin_merge_end
                right_valley = min_cov 
            print("left_bound: ",left_bound)
            print("right_bound: ",right_bound)
            try:
                max_index,local_max = tmp_maximas[peak_index]
            except:
                print(datetime.datetime.now(),"\t\t","no tmp_maxima found, skip this localmax !") # why this exist ?
                # peak_index += 1
                left_bound = right_bound
                left_valley = right_valley
                n_cand_peak += 1
                continue
            if max_index >= gene_length:
                max_index = gene_length-1
            print("",datetime.datetime.now(),"\t\t","start bin_center: ",max_index)


            if x[max_index] <= min_cov:
                print(datetime.datetime.now(),"\t\t","initial cov <= min_cov, skip this bin !")
                # peak_index += 1
                left_bound = right_bound
                left_valley = right_valley
                n_cand_peak += 1
                continue
            

            # 2. define shuffle local region (use full-len tx/global if tx len shorter than 2*ext)
            #s: local start site; end: local end site; ext: bp extension each side, usually > 2*aver_read_len (eg.50)
            # 2.1 local mode
            if mode=="local":
                if max_index>locals_list[-1][0] and max_index<locals_list[-1][1]: # re-use last local boundary if possiable
                    s = locals_list[-1][0]
                    e = locals_list[-1][1]
                    print("re-use last local boundary in local mode")
                    fold = fold_list[-1]
                    print("re-use last local fold in local mode")
                else:
                    # seem cannot aviod re-cal range (localmax center) ?
                    most_left_start = max_index
                    while (most_left_start > left_bound) and (x[most_left_start - 1] >= min_cov): # fixed shuf region might be better and faster
                        most_left_start -= 1
                    most_right_end = max_index
                    while (most_right_end < (gene_length - 1)) and (x[most_right_end + 1] >= min_cov): # fixed shuf region might be better and faster
                        most_right_end += 1
                    
                    locals_list.append([most_left_start,most_right_end]) # most_right_end,most_left_start


                    subset = min(500,(most_right_end-most_left_start)) # ctrl complexity of gini calc: no more than 500 pos/values
                    gini_val = np.nan_to_num(gini(np.random.choice(x[most_left_start:most_right_end], subset, replace=False)))
                    print(datetime.datetime.now(),"\t\t","gini: ",gini_val) # 85% <0.1, 0% >0.9
                    # fold = int(3+2/(0.05+gini_val)) # to avoid zero err: 2/(1+gini) replace 1/gini, hard limit:[3,40]
                    # fold = int(3+40*(1-math.log(gini_val+1))) # before 1229; hard limit:[15,43]
                    fold = int(3+47*(1-math.log2(gini_val+1))) # after 0104; hard limit:[3,50] (no much change for RNAs boundary)
                    fold_list.append(fold)
                    span = fold*(most_right_end-most_left_start)
                    
                    if most_right_end-most_left_start <= 10:
                        print(datetime.datetime.now(),"\t\t","initial local <= 10, skip this bin !")
                        # peak_index += 1
                        left_bound = right_bound
                        left_valley = right_valley
                        n_cand_peak += 1
                        continue

                    s = int(0.5*(most_left_start+most_right_end)-0.5*span)
                    e = int(s+span)
                    print(datetime.datetime.now(),"\t\t","1st perm start pos: ", s)
                    print(datetime.datetime.now(),"\t\t","1st perm end pos: ", e)
                    
                    print(datetime.datetime.now(),"\t\t","tx max_index: ", max_index)
                    print(datetime.datetime.now(),"\t\t","tx max_index - most_left_start: ", max_index-most_left_start)
                    print(datetime.datetime.now(),"\t\t","tx most_right_end - max_index: ", most_right_end-max_index)
                    # print(datetime.datetime.now(),"\t\t","tx RL: ", aver_read_len)
                    print(datetime.datetime.now(),"\t\t","tx span: ", span)
                    print(datetime.datetime.now(),"\t\t",n_cand_peak,"/",n_cand_peaks,"localmax:",local_max)
                        
                    if gene_length <= span:
                        s = 0
                        e = gene_length
                        print(datetime.datetime.now(),"\t\t",n_cand_peak,"/",n_cand_peaks,"local mode beyound boundary, use full-len !")
                        # print(datetime.datetime.now(),"\t\t","local tx_cov: ", tx_cov)
                    else:
                        if s < 0: 
                            s = 0 
                            e = span
                            print(datetime.datetime.now(),"\t\t",n_cand_peak,"/",n_cand_peaks,"left beyound boundary !")
                        elif e > gene_length:
                            s = gene_length - span
                            e = gene_length
                            print(datetime.datetime.now(),"\t\t",n_cand_peak,"/",n_cand_peaks,"right beyound boundary !")
                        print(datetime.datetime.now(),"\t\t",n_cand_peak,"/",n_cand_peaks,"local mode not beyound boundary !")
                        print(datetime.datetime.now(),"\t\t","2nd perm start pos: ", s)
                        print(datetime.datetime.now(),"\t\t","2nd perm end pos: ", e)
                        
                    # 2.1.1 shuffle reads at local region
                    # tx_id = tx_reads.reference[0]
                    tx_local = [tx_id,s,e,"+",tx_id]
                    tx_local2 = [tx_id,most_left_start,most_right_end,"+",tx_id]
                    tx_local_reads = read_tid_frag_from_records(tx_local2, tx_reads) # tx_local: original clipper mode, tx_local2: new mode !!!
                    print(datetime.datetime.now(),"\t\t",n_cand_peak,"/",n_cand_peaks,"local reads number: ", len(tx_local_reads))
                    # print(datetime.datetime.now(),cand_bin_index,"/",n_cand_bins,"local reads record: ", tx_local_reads[0])
                    #tx_cov = _child_get_permutation_fdr(tx_local, tx_local_reads, pval_cutoff=permutate_pval, min_cov=min_cov, max_iter=max_iter, seed=seed)
                    res = do_permutation(transcr=tx_local, read_transcript=tx_local_reads, max_iter=max_iter, pval_cutoff=permutate_pval, min_cov=min_cov, seed=seed) # pseudo-shuffling
                    tx_cov, _height_to_pval = res[0:2]
                    tx_cov = max(tx_cov, min_cov)
                    # print(datetime.datetime.now(),"\t\t",cand_bin_index,"/",n_cand_bins,"local res:",res)
                    print(datetime.datetime.now(),"\t\t",n_cand_peak,"/",n_cand_peaks,"local tx_cov:",tx_cov)
                    



            # 2.2 global mode
            elif mode=="global":
                s = 0
                e = gene_length 
            else:
                raise Exception("only global and local mode for background cov estimation supported !")
            

            # 3. filter bin/localmax
            if local_max-tx_cov <= 0:
                print(datetime.datetime.now(),"\t\t",n_cand_peak,"/",n_cand_peaks,"local_max <= tx_cov, skip this bin !")
                # peak_index += 1
                left_bound = right_bound
                left_valley = right_valley
                n_cand_peak += 1
                continue
            print(datetime.datetime.now(),"\t\t","perm start pos: ", s)
            print(datetime.datetime.now(),"\t\t","perm end pos: ", e)
            
            
            # 4. extend boundary to local/global bg cov     
            # if local_max > min_cov:
            # find bounds when input signal drops below tx_cov
            start = end = max_index   
            # start_trend = max_index 
            max_index_new = max_index
            if boundary == "localmaxdecay":
                while (start > left_bound) and (x[start - 1] >= decay*local_max) and (x[start - 1] >= left_valley): #  decay*(local_max-tx_cov)
                    start -= 1
                # end = max_index
                while (end < right_bound-1) and (x[end + 1] >= decay*local_max) and (x[end + 1] >= right_valley): # (gene_length - 1), decay*(local_max-tx_cov),
                    end += 1
            elif boundary == "background": # default
                while (start > left_bound) and (x[start - 1] >= tx_cov) and (x[start - 1] >= left_valley): #  decay*(local_max-tx_cov)
                    start -= 1

                while (end < right_bound-1) and (x[end + 1] >= tx_cov) and (x[end + 1] >= right_valley): # (gene_length - 1), decay*(local_max-tx_cov),
                    end += 1
                    # #local_max = max(local_max, x[end])
                    # if local_max < x[end]: 
                    #     max_index_new = end
                    #     local_max = x[end] # update local_max if x[end] may hihger than localmax, seem needed
            else:
                raise Exception("only localmaxdecay and background mode for boundary are supported !")
            print(datetime.datetime.now(),"\t\t",n_cand_peak,"/",n_cand_peaks,"start pos:",start)
            print(datetime.datetime.now(),"\t\t",n_cand_peak,"/",n_cand_peaks,"end   pos:",end)



            ## 5. filter peak by Poisson dist. significance test
            tx_peak = [tx_id,start,end,"+",tx_id]
            tx_peak_reads = read_tid_frag_from_records(tx_peak, tx_reads)
            if mode=="local":
                fold = max(10,fold) # require broader local range when poisson test
                half_peak_len = int((fold-1)*0.5*(end-start)) # 0.5 --> 10 (if 0.5, many peaks lost for high poisson_val)
                tmp_local_start = max(0,start-half_peak_len)
                tmp_local_end = min(end+half_peak_len,gene_length)
                tx_peak_local = [tx_id,tmp_local_start,tmp_local_end,"+",tx_id]
                tx_peak_local_reads = read_tid_frag_from_records(tx_peak_local, tx_reads)
                print("tx_peak_local_reads num:",len(tx_peak_local_reads))
                print("reads_in_peak num:",len(tx_peak_reads))
                print("gene_length:",tmp_local_end-tmp_local_start)
                print("peak_length:",end-start)
                poisson_p = poissonP(reads_in_gene=len(tx_peak_local_reads), reads_in_peak=len(tx_peak_reads), gene_length=tmp_local_end-tmp_local_start, peak_length=end-start)
                # lam = min(1, int(len(tx_peak_local_reads)*(end-start)/(tmp_local_end-tmp_local_start)))
            elif mode=="global":
                poisson_p = poissonP(reads_in_gene=len(tx_global_reads), reads_in_peak=len(tx_peak_reads), gene_length=gene_length, peak_length=end-start)
                # lam = min(1, int(len(tx_global_reads)*(end-start)/(gene_length)))
            else:
                raise Exception("only global and local mode supported !")
            print(datetime.datetime.now(),"\t\t",n_cand_peak,"/",n_cand_peaks,"poisson_p:",poisson_p)
            
            # add current peak to results
            if ((end - start) < min_peak_length) or ((end - start) > max_peak_length): #  max_peak_length
                print(datetime.datetime.now(),"\t\t","finished (lenth",(end - start),"not pass filter)",n_cand_peak,"/",n_cand_peaks,"cand_bins !")
            elif (poisson_p > poisson_pval):
                print(datetime.datetime.now(),"\t\t","finished (poisson_pval",poisson_p,"not pass filter",poisson_pval,")",n_cand_peak,"/",n_cand_peaks,"cand_bins !")
            else:
                #print(datetime.datetime.now(),(start, end, local_max))
                # max_index_new,local_max = find_one_maximum_in_center_local(x[start:end]) # 22/10/19: change end+1 to end, 22/10/20 comment
                # max_index_new += start # 22/10/20 comment
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
                print(datetime.datetime.now(),"\t\t",n_cand_peak,"/",n_cand_peaks,"tmp_fdr:",tmp_fdr)
                print(datetime.datetime.now(),"\t\t",n_cand_peak,"/",n_cand_peaks,"cand_bins:", start, end, local_max, max_index_new, tx_cov, tmp_fdr, poisson_p)
                peaks.append([tx_id, start, end, local_max, max_index_new, tx_cov, tmp_fdr, poisson_p]) # why not need deduplicate, adjacent peak not overlap: below left_bound/end updated, also next_cand_bin may jump to far away
                print(datetime.datetime.now(),"\t\t","finished",n_cand_peak,"/",n_cand_peaks,"cand_bins !")

            # find next candidate bin
            # left_bound = end
            left_bound = right_bound
            left_valley = right_valley
            n_cand_peak += 1
        # cand_bin_merge_index += 1    
    print(datetime.datetime.now(),"\t","this tx all locals finished !")
    # print(datetime.datetime.now(),"\t","this tx all peaks:",peaks)
    return peaks


def single_process_get_chrom_peaks(list_):
    single_process_gene_list,single_process_gene_readList,min_cov,bin_width,min_peak_length,max_peak_length,decay,permutate_pval,poisson_pval,max_iter,seed,mode,boundary,nthread,recursive_maxima_elimination = list_
    threads = []
    global peaks_chrom # inside one process
    peaks_chrom = [[] for i in range(nthread)] # [[], [], [], []]     # [None] * nthread # [None, None, None, None]
    start_time = time.time()

    multi_thread_gene_list = [x for x in chunkify(single_process_gene_list, nthread)]
    multi_thread_gene_list = [ele for ele in multi_thread_gene_list if ele != []] # remove potential blank elements: chunkify([1,3,5],2)
    nthread = max(1,len(multi_thread_gene_list)) # update nthreads if blank elements exist
    multi_thread_gene_reads = get_chunkify_readsList(child_gene_list=multi_thread_gene_list, all_tx_reads=single_process_gene_readList)
    
    for thread in range(nthread):
        l1 = [multi_thread_gene_list[thread],multi_thread_gene_reads[thread],min_cov,bin_width,min_peak_length,max_peak_length,decay,permutate_pval,poisson_pval,max_iter,seed,mode,boundary,thread,recursive_maxima_elimination] # seem no map for multi-threads, unlike multi-proc
        t = Thread(target=single_thread_get_chrom_peaks, args=(l1,))
        threads.append(t)
        t.start()
    # wait for the threads to complete
    for t in threads:
        t.join()
    peaks_chrom = [element for sublist in peaks_chrom for element in sublist] # flattern list of diff thread 
    peaks_chrom = [element for sublist in peaks_chrom for element in sublist] # flattern list of diff chr 
    end_time = time.time()
    cur_pid = os.getpid()
    print(f'sub process {cur_pid: 0.0f} run-time: {end_time- start_time: 0.2f} seconds')
    return peaks_chrom

def single_thread_get_chrom_peaks(l):
    # if chrom != 'ENST00000385227_____1': # ENST00000362283_____1, NR_146152_____1, 12965
    #     continue
    single_thread_gene_list,single_thread_gene_readList,min_cov,bin_width,min_peak_length,max_peak_length,decay,permutate_pval,poisson_pval,max_iter,seed,mode,boundary,thread,recursive_maxima_elimination = l
    # min_cov = 5
    # bin_width = 10
    # min_peak_length = 15
    # max_peak_length=200
    # decay = 0.0
    # permutate_pval = 0.05
    # max_iter=50
    # seed=1234
    # local= True

    pid = os.getpid()
    tot = len(single_thread_gene_list)
    # global n_peaks
    n_peaks = 0
    # global bam
    # global chroms # no need to define global, cause not changed, just read
    # peaks_chrom = []
    for i in range(len(single_thread_gene_list)):
        if not i % 200:
            # logger.debug('\n\n','pid %s, threadid %s : %i / %i (%.2f%%)'% (pid, thread, i, tot, float(i)/float(tot)*100))
            print('\n\n','pid %s, threadid %s : %i / %i (%.2f%%)'% (pid, thread, i, tot, float(i)/float(tot)*100))
        chrom = single_thread_gene_list[i]
        size = chroms[chrom]
        print("\n\n",datetime.datetime.now(),"start chrom: ",chrom,n_peaks,"/",len(single_thread_gene_list))
        print(datetime.datetime.now(),"size: ",size)
        
        # this_tx_reads = bam.fetch(str(chrom), start=0, end=size) # , until_eof=True, BAM files with a large number of reference sequences are slow (https://pysam.readthedocs.io/en/latest/faq.html)
        tx_global = [chrom,0,size,"+",chrom] # 0-based
        # this_tx_reads = read_tid_frag_from_bam(tid=tx_global, reads=this_tx_reads, is_unique=True) # convert 0-based bam iterable reads to 0-based reads records
        this_tx_reads = read_tid_frag_from_readList(tid=tx_global, reads=single_thread_gene_readList, is_unique=True) # convert 0-based bam iterable reads to 0-based reads records
        count_signal = count_pileup_heights(tlen=size, reads=this_tx_reads, downsample=False) # should not downsample ?
        if max(count_signal) <= min_cov:
            n_peaks += 1
            print(datetime.datetime.now(),"skip, low tx reads depth: ",max(count_signal))
            print("",datetime.datetime.now(),"end chrom: ",chrom,n_peaks,"/",len(single_thread_gene_list))
            continue 
        # count_signal = np.nan_to_num(bigwig.values(chrom, 0, size))
        print(datetime.datetime.now(),"count_signal: ", count_signal)
        print(datetime.datetime.now(),"sum(count_signal): ", sum(count_signal))
        print(datetime.datetime.now(),"len(this_tx_reads): ",len(this_tx_reads))

        _peaks_chrom = _call_peaks_localmax(tx_id=chrom, count_signal=count_signal, tx_reads=this_tx_reads, 
            min_peak_length=min_peak_length, max_peak_length=max_peak_length, bin_width=bin_width, min_cov=min_cov, decay=decay, mode=mode, boundary=boundary, permutate_pval=permutate_pval, poisson_pval=poisson_pval, max_iter=max_iter, seed=seed, recursive_maxima_elimination=recursive_maxima_elimination)
        if len(_peaks_chrom)>0:
            peaks_chrom[thread].append(_peaks_chrom) # _peaks_chrom is a list
        print("",datetime.datetime.now(),"end chrom: ",chrom,n_peaks,"/",len(single_thread_gene_list))
    # return peaks_chrom


@command_handler
def call_peaks_localmax(args):
    # import pyBigWig
    import pysam
    import numpy as np

    input_bam = args.input_bam
    min_cov = args.min_cov # 5
    bin_width = args.bin_width #10
    min_peak_length = args.min_peak_length #15
    max_peak_length = args.max_peak_length #200
    decay = args.decay #0.0
    permutate_pval = args.permutate_pval #0.05
    poisson_pval = args.poisson_pval #0.05
    max_iter = args.max_iter #100
    seed = args.seed #1234
    mode = args.mode # local
    boundary = args.boundary # background
    nprocess = max(int(args.process),1)
    nthread = max(int(args.thread),1)
    recursive_maxima_elimination = args.recursive_maxima_elimination
    output_file = args.output_file
    random.seed(seed)
    
    print('\t','read bam file: ',input_bam)
    print('\t','permutate_pval: ',permutate_pval)
    print('\t','poisson_pval: ',poisson_pval)
    print('\t','mode: ',mode)
    print('\t','boundary: ',boundary)
    print('\t','read seed: ',seed)
    print('\t','recursive_maxima_elimination: ',recursive_maxima_elimination)
    print('\t','min_cov: ',min_cov)
    print('\t','nprocess: ',nprocess)
    print('\t','nthread: ',nthread)
    # logger.info('\t' + 'read bw file: ' + args.input_bw)
    # print("",datetime.datetime.now(),"nthread: ",nthread)
    # print("",datetime.datetime.now(),"nprocess: ",nprocess)
    bam=pysam.AlignmentFile(input_bam, 'rb')
    # bigwig = pyBigWig.open(args.input_bw) # read pre-cal bw, or use clipper.readsToWiggle ?
    # chroms = bigwig.chroms()
    global chroms
    chroms={x['SN']:x['LN'] for x in bam.header['SQ']}
    
    valid_chroms=set()

    # valid_chroms.add("ENST00000385227_____1") # miR eg
    # valid_chroms.add("NR_023363_____1") # missed 1st peak by broad local definition
    # valid_chroms.add("NR_146144_____1") 
    # valid_chroms.add("24596")
    # valid_chroms.add("25812")
    # valid_chroms.add("ENST00000385243_____1")
    # valid_chroms.add("promoter__chr5___140657847____140660847_pos")


    
    for read in bam.fetch():
        valid_chroms.add(bam.get_reference_name(read.reference_id))
    print("",datetime.datetime.now(),"valid_chroms: ",len(valid_chroms),"/",len(chroms))
    chroms = {x:y for x,y in chroms.items() if x in valid_chroms}
    bam.reset()  
    gene_list = list(chroms.keys()) # py2: list chroms.keys()
    multi_process_gene_list = [x for x in chunkify(gene_list, nprocess)]
    multi_process_gene_list = [ele for ele in multi_process_gene_list if ele != []] # remove potential blank elements: chunkify([1,3,5],2)
    nprocess = max(1, len(multi_process_gene_list)) # update nprocess 
    all_tx_reads = get_readsList_from_bam(gene_list=gene_list, bam=bam) 
    multi_process_gene_reads = get_chunkify_readsList(child_gene_list=multi_process_gene_list, all_tx_reads=all_tx_reads)
    bam.close()
    # print(datetime.datetime.now(),"\t","all_tx_reads: ", all_tx_reads)
    # print(datetime.datetime.now(),"\t","child_gene_reads: ", child_gene_reads)
    print(datetime.datetime.now(),"\t","len(child_gene_list): ", len(multi_process_gene_list))
    print(datetime.datetime.now(),"\t","len(multi_process_gene_reads): ", len(multi_process_gene_reads))

    if nprocess>1:
        pool = Pool(processes=nprocess)
        pool_res = pool.map(
            single_process_get_chrom_peaks, 
            [[multi_process_gene_list[process], multi_process_gene_reads[process],  min_cov,bin_width,min_peak_length,max_peak_length,decay,permutate_pval,poisson_pval,max_iter,seed,mode,boundary,nthread,recursive_maxima_elimination] for process in range(nprocess)]
        ) # map returns an iterator
        pool.terminate() # close ?
        pool.join()
        all_peaks_chrom=[]
        for process in range(nprocess):
            for j in range(len(pool_res[process])):
                all_peaks_chrom.append(pool_res[process][j]) # pool_res[i][j]: list
    else:
        all_peaks_chrom = single_process_get_chrom_peaks([gene_list, all_tx_reads,  min_cov,bin_width,min_peak_length,max_peak_length,decay,permutate_pval,poisson_pval,max_iter,seed,mode,boundary,nthread,recursive_maxima_elimination])
    # peaks_chrom = [element for sublist in peaks_chrom for element in sublist] # flattern list of diff process 
    all_peaks_chrom = [ele for ele in all_peaks_chrom if ele != []] # drop blank ele
    # print(datetime.datetime.now(),"\t","peaks_chrom: ", peaks_chrom)
    
    bed = open(output_file, 'w')
    # logger.info("\t" + 'write bed file: ' + output_file)
    print("\t",'write bed file: ',output_file)
    n_peaks = 0 # sample-wise
    for peak in all_peaks_chrom:
        n_peaks += 1
        #peaks.append([chrom, peak[0], peak[1], 'peak_%d'%n_peaks, peak[2], '+'])
        # bed.write('%s\t%d\t%d\tpeak_%d\t%d\t+\t%d\t%d\t%.10f\n'%(chrom, peak[0], peak[1], n_peaks, peak[2], peak[3], peak[4], peak[5]))
        # tx_id,start,end,peak_n,localmax,(+),maxpos,bg_cov,permutate_p,poisson_p
        bed.write('%s\t%d\t%d\tpeak_%d\t%d\t+\t%d\t%d\t%.10f\t%.10f\n'%(peak[0], peak[1], peak[2], n_peaks, peak[3], peak[4], peak[5], peak[6], peak[7]))
    # bigwig.close()
    bed.close()


if __name__ == '__main__':
    main_parser = argparse.ArgumentParser(description='Call peaks from exRNA tbam')
    subparsers = main_parser.add_subparsers(dest='command')
    
    parser = subparsers.add_parser('call_peaks_localmax')
    parser.add_argument('--input-bam', type=str, required=True,
        help='input Bowtie2 Bam file of raw reads')
    # parser.add_argument('--input-bw', type=str, required=True,
    #     help='input BigWig file of raw reads coverage')
    parser.add_argument('--min-peak-length', type=int, default=10,
        help='minimum length required for a peak')
    parser.add_argument('--max-peak-length', type=int, default=500,
        help='maximum length required for a peak')
    parser.add_argument('--decay', type=float, default=0.5,
        help='decay factor of peak summit to define whether two adjacent localmax belong to one')
    parser.add_argument('--permutate-pval', type=float, default=0.05,
        help='pval for bg cov permutation estimation')
    parser.add_argument('--poisson-pval', type=float, default=0.05,
        help='pval for poisson dist significance test, set it to 1 if not filter peak by count poisson test (and this may find more peaks with low coverage)')
    parser.add_argument('--min-cov', type=float, default=5,
        help='minimum coverage required to define a peak')
    parser.add_argument('--bin-width', type=int, default=10,
        help='bin width to search enriched bins')
    parser.add_argument('--max-iter', type=int, default=100,
        help='max iteration num to estimate bg cov')
    parser.add_argument('--seed', type=int, default=1234,
        help='seed for random permutation')
    parser.add_argument('--output-file', '-o', type=str, required=True,
        help='output peaks in BED format')
    parser.add_argument('--mode', type=str, default="local",
        help='call peak permutation mode when estimate bg cov: local (default) or global')
    parser.add_argument('--recursive-maxima-elimination', type=bool, default=True, 
        help='whether or not to recursive-maxima-elimination each local maxima. default: True; optional: False')
    parser.add_argument('--boundary', type=str, default="background",
        help='call peak boundary position: background (default) or localmaxdecay')
    parser.add_argument('--process', type=int, default=1,
        help='num of processes to call peak in parallel, default: 1')
    parser.add_argument('--thread', type=int, default=max(int(0.5*os.cpu_count()), 1),
        help='num of threads to call peak in parallel (in each process), default: half num of cpu dectected')

    # input_bam = args.input_bam
    # min_cov = args.min_cov # 5
    # bin_width = args.bin_width #10
    # min_peak_length = args.min_peak_length #15
    # max_peak_length = args.max_peak_length #200
    # decay = args.decay #0.0
    # permutate_pval = args.permutate_pval #0.05
    # max_iter = args.max_iter #50
    # seed = args.seed #1234
    # mode = args.local # True
    # nthread = max(int(args.threads),1)
    # output_file = args.output_file

    args = main_parser.parse_args()
    if args.command is None:
        raise ValueError('empty command')
    # logger = logging.getLogger('call_peak.' + args.command)

    command_handlers.get(args.command)(args)