#!/usr/bin/env python3

# input fa id rows should not contain blank space
# todo: add two layer parallel in priority like IntaRNA, current at records. at RNAfold if limited records

import os
import shutil
import subprocess
import re
import sys
from glob import glob
import numpy as np
from scipy.stats import norm
from scipy.stats.mstats import zscore
import multiprocess
import time
import datetime
import csv
start_time=time.time()
sys.stderr.write('Start time: %s\n' %(datetime.datetime.now()))
#sys.stderr.write('Pysam version used: %s\n' %(pysamVersion))

FAFILE = str(sys.argv[1])
DIR = os.path.dirname(FAFILE)
SHUF_TIMES = sys.argv[2] # 100
SEED = sys.argv[3] # 1234
OUTFILE = sys.argv[4] # 
CORES = sys.argv[5] # cores
TRASH = open(os.devnull)


##################################
def file2string(filename):
    seqtitle=[]
    primary_sequence=[]
    lines=open(filename,'r+').readlines()
    for line in lines:
        if '>' in line:
            seqtitle.append(line.strip())
        else:
            primary_sequence.append(line.strip())
    return seqtitle,primary_sequence

###################################
#cmd = f"multiperm -n 1000 -w {FAFILE}"
#/BioII/lulab_b/jinyunfan/anaconda3/envs/bioinfo_py36/bin/fasta-dinucleotide-shuffle-py3 -c 100 -t shuffle -f /BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/lulab/cf_FTA_only.txt.preMIR.fa
#cmd = f"fasta-dinucleotide-shuffle-py3 -s {SEED} -c {SHUF_TIMES} -f {FAFILE} > {FAFILE}_perm"
cmd = f"/BioII/lulab_b/jinyunfan/anaconda3/envs/bioinfo_py36/bin/fasta-dinucleotide-shuffle-py3 -s {SEED} -c {SHUF_TIMES} -f {FAFILE} > {FAFILE}_perm"
p = subprocess.Popen(cmd, shell=True)
p.wait()
#shuf each fa record until SHUF_TIMES, then next fa record !


###################################
def vienna_rnafold(seq):
    command=['/BioII/lulab_b/jinyunfan/anaconda3/envs/bioinfo_py36/bin/RNAfold','-d2','--noLP'] # '--jobs 4 '
    rna_fold=subprocess.Popen(command,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
    rna_fold_output=rna_fold.communicate(seq.encode('utf-8'))[0] # seem no need to encode
    #print(rna_fold_output)
    if len(rna_fold_output) > len(seq):
        energy=re.findall("(-?\d+\.\d+)", rna_fold_output.decode('UTF-8') ) # re.findall return list of str: ['-18.0']
        dotbracket=re.findall("(?<=\n)[(.][(.)]*(?=\s\(.*\d)", rna_fold_output.decode('UTF-8') ) # ['..(..(...)).']
        #(?<=\n)[(.][(.)]*(?=\s\(.*\d)
    else:
        print('Error encountered while doing rnafold. Skipping...')
        energy=None
        dotbracket=None
    return energy,dotbracket

###################################
seqtitle,primary_sequence=file2string(FAFILE)
shufffledseqtitle,shuffled_sequence=file2string(FAFILE+"_perm")
pool=multiprocess.Pool(int(CORES)) # int(multiprocess.cpu_count()*0.5)) # CORES
#energy_list=pool.map(vienna_rnafold,primary_sequence) # single return: [['-10.2'],['-14.00'],['0.1']]
energy_list,dotbracket_list=zip(*pool.map(vienna_rnafold,primary_sequence)) # multi return: [(['-10.2'],['..(..)..']),(['-31.0'],['.(...)..'])]
#energy_list: (['-10.2'], ['-31.0'])
#dotbracket_list: (['..(..)..'], ['.(...)..'])
energy_list=list(energy_list)
dotbracket_list=list(dotbracket_list)
shuffled_energy_list,_=zip(*pool.map(vienna_rnafold,shuffled_sequence)) # a list of result-tuples from your pool.map() call
shuffled_energy_list=list(shuffled_energy_list)
# shuffled_energy_list,_=[x[0] for x in tmp_list],[x[1] for x in tmp_list]
#print(shuffled_energy_list)




####################################
# calculate z-score and p-values
flat_energy=[x for xs in energy_list for x in xs] # energy_list=[['-10.2'],['-14.00'],['0.1']] --> ['-10.2', '-14.00', '0.1']
flat_dotbracket=[x for xs in dotbracket_list for x in xs] # dotbracket_list=(['..(..)..'], ['.(...)..']) --> 
flat_shuffle_energy=[x for xs in shuffled_energy_list for x in xs]
index=0
print(len(flat_energy))
print(len(flat_dotbracket))
print(len(flat_shuffle_energy))
print(flat_energy)
print(flat_dotbracket)


pvalue_list=[]
rnafold_list=[]
for energy in flat_energy:
    energies=flat_shuffle_energy[index:index+int(SHUF_TIMES)]
    energies.append(energy)
    a = np.array(energies).astype(float) 
    z = np.nan_to_num(zscore(a))[-1]  # avoid zero std, stil warning: RuntimeWarning: invalid value encountered in true_divide
    #pvalue = norm.sf(abs(z)) * 2
    # pvalue = norm.cdf(-1*abs(z)) * 2 
    pvalue = norm.cdf(z)
    pvalue_list.append(pvalue)
    rnafold_list.append(energy)
    index=index+int(SHUF_TIMES)
#print(len(pvalue_list))

# dotbracket_list2=[]
# for dotbracket in flat_dotbracket:
# 	dotbracket_list2.append(dotbracket)

# print(len(dotbracket_list2))

# res = {}
# for key in seqtitle:
#     for value in pvalue_list:
#         res[key] = value
#         pvalue_list.remove(value)
#         break  

# open file for writing, "w" is writing
w = csv.writer(open(OUTFILE, "w"))
# loop over dictionary keys and values
# for key, val in res.items():
#     # write every key and value to file
#     w.writerow([key, val])
for i in range(0,len(seqtitle)):
    # write every key and value to file
    w.writerow([seqtitle[i], pvalue_list[i], rnafold_list[i], flat_dotbracket[i]])
    
#print(pvalue_list)
#print(f"original_file: {ORIGINAL}")
#print(f"original_MFE: {ENERGY} kcal/mol")
#print(f"shuffle_times: {SHUF_TIMES}")
#print(f"p-value: {pvalue}")
#print("--- %s seconds ---" % (time.time() -start_time))
sys.stderr.write('Time (s): %s\n' %(time.time() -start_time))
sys.stderr.write('End time: %s\n' %(datetime.datetime.now()))

for out_junk in glob('*.ps'):
    os.remove(out_junk)
TRASH.close()
#os.remove(FAFILE+"_perm")
