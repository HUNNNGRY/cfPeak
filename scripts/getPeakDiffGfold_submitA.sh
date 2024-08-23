#!/bin/bash
#SBATCH -J Gfold
#SBATCH -p Acluster
#SBATCH -n 6
#SBATCH --cpus-per-task 1
#SBATCH --output=%j.out
#SBATCH --error=%j.err

source activate cfDNA_base
dst="DATASET" #"TCGA_small7"
outSuf="SURFIX" #"_tumorigenesis"
cell="CELLTYPE" # HNSC

cd /BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/$dst/gfold/count
mkdir -p ../diff/log

for smp in `cat ../${cell}${outSuf}_s2.csv | tr "," " "`
do
  if [ -s ../diff/${cell}_${smp}${outSuf}.diff ];
  then
    echo "skip $smp"
  else
    echo "start $smp at `date`"
    gfold diff -s1 `cat ../${cell}${outSuf}_s1.csv` -s2 $smp -sc 0.01 -suf .read_cnt -o ../diff/${cell}_${smp}${outSuf}.diff > ../diff/log/${cell}${outSuf}_${smp}.log 2>&1
  fi
done

