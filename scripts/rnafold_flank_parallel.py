#!/usr/bin/env python3
# not finished yet !!!!!!
# todo: use flank/shuf regions as bg dist, cal peak mfe pval

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
import csv
start_time=time.time()


# DIR = os.path.dirname(FAFILE)
FAFILE = str(sys.argv[1])
FAFILE_SHUF = sys.argv[2] # shuf bed

SHUF_TIMES = sys.argv[3] # 100
SEED = sys.argv[4] # 1234
OUTFILE = sys.argv[5] # 
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
# cmd = f"/BioII/lulab_b/jinyunfan/anaconda3/envs/bioinfo_py36/bin/fasta-dinucleotide-shuffle-py3 -s {SEED} -c {SHUF_TIMES} -f {FAFILE} > {FAFILE}_perm" # shuffle peak self
# cmd = f"/BioII/lulab_b/jinyunfan/anaconda3/envs/bioinfo_py36/bin/fasta-dinucleotide-shuffle-py3 -s {SEED} -c {SHUF_TIMES} -f {FAFILE_SHUF} > {FAFILE_SHUF}_perm" # shuffle peak flank/bg
# p = subprocess.Popen(cmd, shell=True)
# p.wait()

###################################
def vienna_rnafold(seq):
    command=['/BioII/lulab_b/jinyunfan/anaconda3/envs/bioinfo_py36/bin/RNAfold','-d2','--noLP']
    rna_fold=subprocess.Popen(command,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
    rna_fold_output=rna_fold.communicate(seq.encode('utf-8'))[0]
    #print(rna_fold_output)
    if len(rna_fold_output) > len(seq):
        energy=re.findall("(-?\d+\.\d+)",str(rna_fold_output))
    else:
        print('Error encountered while doing rnafold. Skipping...')
        energy=None
    return energy

###################################
seqtitle,primary_sequence=file2string(FAFILE)
shufffledseqtitle,shuffled_sequence=file2string(FAFILE_SHUF)
pool=multiprocess.Pool(multiprocess.cpu_count())
energy_list=pool.map(vienna_rnafold,primary_sequence)
shuffled_energy_list=pool.map(vienna_rnafold,shuffled_sequence)
#print(shuffled_energy_list)

####################################
# calculate z-score and p-values
flat_energy=[x for xs in energy_list for x in xs]
flat_shuffle_energy=[x for xs in shuffled_energy_list for x in xs]
index=0

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
    w.writerow([seqtitle[i], pvalue_list[i], rnafold_list[i]])
    
#print(pvalue_list)
#print(f"original_file: {ORIGINAL}")
#print(f"original_MFE: {ENERGY} kcal/mol")
#print(f"shuffle_times: {SHUF_TIMES}")
#print(f"p-value: {pvalue}")
#print("--- %s seconds ---" % (time.time() -start_time))


for out_junk in glob('*.ps'):
    os.remove(out_junk)
TRASH.close()
#os.remove(FAFILE+"_perm")
