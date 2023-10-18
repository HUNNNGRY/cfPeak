import gzip
import sys
import re
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input_file',help='Path and prefix of  input *_[12].fastq.gz file',required=True)
parser.add_argument('-o','--output_file',help='Path and prefix of output  *_[12].fastq.gz file',required=True)
parser.add_argument('-m','--min_length',help='Reads shorter than this value will be discarded',type=int,default=30)
parser.add_argument('-s','--strandness',help='Strandness of input reads',default="forward")
args = parser.parse_args()

matemap = {'1':'1','2':'2'} if args.strandness=="forward" else {'1':'2','2':'1'}

output_file_1 = args.output_file+'_{}.fastq.gz'.format(matemap["1"])
output_file_2 = args.output_file+'_{}.fastq.gz'.format(matemap["2"])
input_file_1 = args.input_file+'_{}.fastq.gz'.format(matemap["1"])
input_file_2 = args.input_file+'_{}.fastq.gz'.format(matemap["2"])

min_length=args.min_length
patterns = ['^[GC]+','[GC]+$'] 

input_1 =  gzip.open(input_file_1)
input_2 =  gzip.open(input_file_2)
output_1 = gzip.open(output_file_1,'w')
output_2 = gzip.open(output_file_2,'w')

count = 0
too_short = 0
lengths_1 = []
lengths_2 = []

for label_1,label_2 in zip(input_1,input_2):
    count += 1
    if count%40000 == 0:
        sys.stderr.write('{} reads processed ... '.format(count/4))
    sequence_1 = next(input_1).strip()
    sequence_2 = next(input_2).strip()
    match_1 = re.search(patterns[0],sequence_1.decode())
    match_2 = re.search(patterns[1],sequence_2.decode())
    length_1 = 0 if match_1 is None else len(match_1.group(0))
    length_2 = 0 if match_2 is None else len(match_2.group(0))
    if len(sequence_1)-length_1 < 30 or len(sequence_2)-length_2 < 30:
        too_short += 1
        _ = next(input_1)
        _ = next(input_2)
        _ = next(input_1)
        _ = next(input_2)
        continue
    sequence_1 = sequence_1[length_1:]
    sequence_2 = sequence_2 if length_2 == 0 else sequence_2[:-length_2]
    output_1.write(label_1)
    output_2.write(label_2)
    output_1.write(sequence_1+b'\n')
    output_2.write(sequence_2+b'\n')
    output_1.write(next(input_1))
    output_2.write(next(input_2))
    quality_1 = next(input_1).strip()[length_1:]
    quality_2 = next(input_2).strip() if length_2 == 0 else next(input_2).strip()[:-length_2]
    output_1.write(quality_1+b'\n')
    output_2.write(quality_2+b'\n')
    lengths_1.append(length_1)
    lengths_2.append(length_2)
lengths_1 = np.array(lengths_1)
lengths_2 = np.array(lengths_2)
stats_1 = np.unique(lengths_1,return_counts=True)
stats_2 = np.unique(lengths_2,return_counts=True)
print('Input reads\tOutput reads\tDisard reads')
print('{}\t{}\t{}\n'.format(count,count-too_short,too_short))
print('length\tnumber_1')
for n,c in zip(stats_1[0],stats_1[1]):
    print('{}\t{}'.format(n,c))
print('length\tnumber_2')
for n,c in zip(stats_2[0],stats_2[1]):
    print('{}\t{}'.format(n,c))
output_1.close()
output_2.close()


        
       
 
       

