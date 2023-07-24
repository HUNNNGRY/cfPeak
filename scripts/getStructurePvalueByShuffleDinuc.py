#!/usr/bin/env python3

"""
Wrapper for ViennaRNA:RNAfold and meme:fasta-dinucleotide-shuffle-py3 to check for significance of a fasta file (one record).
original adapted from https://github.com/klamkiew/cov_trs_structure

Usage:
    python3 bin/getStructurePvalueByShuffleDinuc.py <input.fa> <shuffle_times>

Output:
    original_MFE
    p-value
    ...
    
Dependencies:
    ViennaRNA: RNAfold installed in your $PATH variable.
    (conda) meme:fasta-dinucleotide-shuffle-py3 installed in your $PATH variable.

Note:
- currently only support input fa file with one record

Test run:
% cat test.fa
>hsa-let-7d
CCUAGGAAGAGGUAGUAGGUUGCAUAGUUUUAGGGCAGGGAUUUUGCCCACAAGGAGGUAACUAUACGACCUGCUGCCUUUCUUAGG

% RNAfold -d2 --noLP test.fa # -p -d2 --noLP
>hsa-let-7d
CCUAGGAAGAGGUAGUAGGUUGCAUAGUUUUAGGGCAGGGAUUUUGCCCACAAGGAGGUAACUAUACGACCUGCUGCCUUUCUUAGG
(((((((.((((((((((((((.((((((...((((((.....))))))..........)))))).))))))))))))))))))))) (-42.70)

% python3 bin/getStructurePvalueByShuffleDinuc.py test.fa 500
original_file: test.fa
original_MFE: -42.70 kcal/mol
shuffle_times: 500
p-value: 6.588972465456299e-05

"""

import os
import shutil
import subprocess
import re
import sys

from glob import glob

import numpy as np
from scipy.stats import norm
from scipy.stats.mstats import zscore


ORIGINAL = sys.argv[1]
DIR = os.path.dirname(ORIGINAL)
SHUF_TIMES =  sys.argv[2] # 100
RNAFOLD_PARAMS = "-d2 --noLP" # "-p -d2 --noLP"

TRASH = open(os.devnull)


###################################
# original MFE calculation
#cmd = f"bash -c 'RNAalifold -r --cfactor 0.6 --nfactor 0.5 --noLP --noPS {ORIGINAL}'"
cmd = f"bash -c 'RNAfold {RNAFOLD_PARAMS} {ORIGINAL}'"
p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=TRASH)
output = p.stdout.read().decode('ascii').split('\n')[2]

ENERGY = re.findall("(-?\d+\.\d+)", output)[0]
#https://stackoverflow.com/questions/28467446/problems-with-using-re-findall-in-python
#\d+\.\d+   matches one or more digit characters pus a literal dot plus one or more digits

energies = [ENERGY]

###################################
#cmd = f"multiperm -n 1000 -w {ORIGINAL}"
#/BioII/lulab_b/jinyunfan/anaconda3/envs/bioinfo_py36/bin/fasta-dinucleotide-shuffle-py3 -c 100 -t shuffle -f /BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/lulab/cf_FTA_only.txt.preMIR.fa
cmd = f"/BioII/lulab_b/jinyunfan/anaconda3/envs/bioinfo_py36/bin/fasta-dinucleotide-shuffle-py3 -c {SHUF_TIMES} -f {ORIGINAL} > {ORIGINAL}_perm"
p = subprocess.Popen(cmd, shell=True)
p.wait()

###################################
# move tmp alignments to the input directory
# for data in glob("perm*.aln"):
#   if os.path.exists(f"{DIR}/{data}"):
#     os.remove(f"{DIR}/{data}")
#   shutil.move(data, DIR)

###################################
# for each shuffled alignment
# determine the MFE value

#for aln in glob(f"{ORIGINAL}.perm"):
aln=ORIGINAL+"_perm"
cmd = f"bash -c 'RNAfold {RNAFOLD_PARAMS} {aln}'"
p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=TRASH)
output = p.stdout.read().decode('ascii').split('\n')
#for i in range(2,(SHUF_TIMES-1)*3,3):
energies.extend(re.findall("(-?\d+\.\d+)", str(output)))

os.remove(aln)  # removing all tmp alignments.
# for out_junk in glob("*.ps"):
#     os.remove(out_junk)
TRASH.close()

###################################
# calculate z-score and p-values
a = np.array(energies).astype(float)
z = zscore(a)[0]
#pvalue = norm.sf(abs(z)) * 2
#pvalue = norm.cdf(-1*abs(z)) * 2
pvalue = norm.cdf(z)

# print(f"output: {energies}")
print(f"original_file: {ORIGINAL}")
print(f"original_MFE: {ENERGY} kcal/mol")
print(f"shuffle_times: {SHUF_TIMES}")
print(f"p-value: {pvalue}")
