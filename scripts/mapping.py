#!/usr/bin/env python
import subprocess
import argparse
import os
import io
import time
# os.environ['JOBLIB_TEMP_FOLDER'] = '/BioII/lulab_b/baopengfei/tmp' # seems not working

def align(fastq,index,bam,unmapped,log,threads,mode,multimap,clipmode,maxIns,minIns):
    # original version: --sensitive not --very-fast
    mapping_cmd = ["bowtie2","-p",str(threads),"-I",str(minIns),"-X",str(maxIns),"-k",str(multimap),"--norc",str("--"+mode),str("--"+clipmode),"--no-unal","--un-gz",unmapped,"-x",index,fastq,"-S","-"]
    #/BioII/lulab_b/jinyunfan/anaconda3/envs/bioinfo_py36/bin/bowtie2 (version 2.3.5)
    print(" ".join(mapping_cmd))
    ps = subprocess.Popen(mapping_cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    sort_cmd = ["samtools","sort","-@",str(threads),"--output-fmt","BAM","-o",bam] # "-@",str(threads),
    print(" ".join(sort_cmd))
    _ = subprocess.check_output(sort_cmd, stdin=ps.stdout)
    ps.wait()
    
    os.system(" ".join(list(map(str, ["samtools","index","-@",str(threads),bam]))))
    # _ = subprocess.run(["samtools","index","-@",str(threads),bam])
    with open(log,"w") as f:
        for line in io.TextIOWrapper(ps.stderr, encoding="utf-8"):
            f.write(line)
    return ps.poll() 

#spikein, univec, rRNA, lncRNA, miRNA, mRNA, piRNA, snoRNA, snRNA, srpRNA, tRNA, tucpRNA, Y RNA, circRNA, genome
def main():
    parser = argparse.ArgumentParser(description='Align/map SE reads to gn/tx reference using bowtie2')
    parser.add_argument('--fastq','-f',required=True,help="Input single end cleaned fastq")
    parser.add_argument('--bam-dir','-bd',required=True,help="Dir for output bam files")
    parser.add_argument('--unmapped-dir','-fd',required=True,help="Dir for output unmapped fastq files")
    parser.add_argument('--index-dir','-id',default="index",help="Dir contains bowtie index")
    parser.add_argument('--log-dir','-ld',required=True,help="Dir for output log files")
    parser.add_argument('--priority','-p',
                        default="spikein_small,univec,rRNA,miRNA,lncRNA,mRNA,piRNA,snoRNA,snRNA,srpRNA,tRNA,tucpRNA,Y_RNA,intron.for,intron.rev,promoter.for,promoter.rev,enhancer.for,enhancer.rev",help="Priority of mapping in left-->right order")
    parser.add_argument('--threads','-t',type=int,default=1,help="Number of threads for mapping")
    parser.add_argument('--mode','-m',default="very-fast",help="map stringent mode, very-fast (default) , fast, sensitive, OR very-sensitive")
    parser.add_argument('--clip-mode',default="end-to-end",
                        help="end clip mod. end-to-end OR local (this default: end-to-end, bowtie2 default: --end-to-end")
    parser.add_argument('--multimap-max',type=int,default=100,
                        help="max number of multimap align record each read threads for mapping by bowtie(2) -k. (this default: 100, bowtie2 default: search for multi-align but only report best single align)")
    parser.add_argument('--maxIns',type=int,default=500,
                        help="max insert size. (valid for PE, this default: 500, bowtie2 default: 500")
    parser.add_argument('--minIns',type=int,default=10,
                        help="min insert size. (valid for PE, this default: 10, bowtie2 default: 0")
    # better not filter insert length in bam, you can filter after this step
    
    args = parser.parse_args()
    fastq = args.fastq
    if not os.path.exists(args.bam_dir):
        os.makedirs(args.bam_dir)
    if not os.path.exists(args.unmapped_dir):
        os.makedirs(args.unmapped_dir)
    if not os.path.exists(args.log_dir):
        os.makedirs(args.log_dir)
    for sequence in args.priority.split(","):
        print(time.ctime()+": Start align {}".format(sequence))
        bam = os.path.join(args.bam_dir,sequence+".bam")
        unmapped = os.path.join(args.unmapped_dir,sequence+".fastq.gz")
        index = os.path.join(args.index_dir,sequence)
        log = os.path.join(args.log_dir,sequence+".log")
        assert align(fastq,index,bam,unmapped,log,threads=args.threads,mode=args.mode,multimap=args.multimap_max,clipmode=args.clip_mode,maxIns=args.maxIns,minIns=args.minIns) == 0
        print(time.ctime()+": Done .")
        fastq = unmapped
        
if __name__ == "__main__":
    main() 
