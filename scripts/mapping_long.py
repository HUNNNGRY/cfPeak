#!/usr/bin/env python
import subprocess
import argparse
import os
import io
import time
# os.environ['JOBLIB_TEMP_FOLDER'] = '/BioII/lulab_b/baopengfei/tmp' # seems not working

def align(fastq1,fastq2,index,bam,unmapped,log,threads,lib,mode,multimap,clipmode,maxIns,minIns):
    # https://likit.github.io/running-bowtiebowtie2-rsem-and-tophat-on-dutp-strand-specific-reads.html
    # for reverse stranded library like dUTP, use --fr --nofw 
    # in paired-end mode, --nofw and --norc mean the fragments (R1 strand?)
    if lib == "reverse":
        s = "--nofw"
    elif lib == "forward":
        s = "--norc"
    mapping_cmd = ["bowtie2","-p",str(threads),"-I",str(minIns),"-X",str(maxIns),"-k",str(multimap),"--fr",s,"--no-discordant","--no-mixed",str("--"+mode),str("--"+clipmode),"--no-unal","--un-conc-gz",unmapped,"-x",index,"-1",fastq1,"-2",fastq2,"-S","-"]
    #/BioII/lulab_b/jinyunfan/anaconda3/envs/bioinfo_py36/bin/bowtie2 (version 2.3.5)
    print(" ".join(mapping_cmd))
    ps = subprocess.Popen(mapping_cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    sort_cmd = ["samtools","sort","-@",str(threads),"--output-fmt","BAM","-o",bam]
    print(" ".join(sort_cmd))
    _ = subprocess.check_output(sort_cmd, stdin=ps.stdout)
    ps.wait()

    # os.system(" ".join(list(map(str, ["samtools","sort","-@",str(threads),"--output-fmt","BAM","-o",bam]))))
    index_cmd = ["samtools","index","-@",str(threads),bam]
    print(" ".join(index_cmd))
    os.system(" ".join(list(map(str, index_cmd))))
    with open(log,"w") as f:
        for line in io.TextIOWrapper(ps.stderr, encoding="utf-8"):
            f.write(line)
    return ps.poll() 




#spikein, univec, rRNA, lncRNA, miRNA, mRNA, piRNA, snoRNA, snRNA, srpRNA, tRNA, tucpRNA, Y RNA, circRNA, genome
def main():
    parser = argparse.ArgumentParser(description='Align/map PE reads to gn/tx reference using bowtie2')
    parser.add_argument('--fastq1','-f1',required=True,help="Input paired end cleaned fastq1")
    parser.add_argument('--fastq2','-f2',required=True,help="Input paired end cleaned fastq2")
    parser.add_argument('--bam-dir','-bd',required=True,help="Dir for output bam files")
    parser.add_argument('--unmapped-dir','-fd',required=True,help="Dir for output unmapped fastq files")
    parser.add_argument('--index-dir','-id',default="index",help="Dir contains bowtie index")
    parser.add_argument('--log-dir','-ld',required=True,help="Dir for output log files")
    parser.add_argument('--priority','-p',default="spikein,univec,rRNA,miRNA,lncRNA,mRNA,piRNA,snoRNA,snRNA,srpRNA,tRNA,tucpRNA,Y_RNA,intron.for,intron.rev,promoter.for,promoter.rev,enhancer.for,enhancer.rev",help="Priority of mapping in left-->right order")
    parser.add_argument('--threads','-t',type=int,default=1,help="Number of threads for mapping")
    parser.add_argument('--strand','-s',default="reverse",help="library strand type, reverse (default) OR forward") # required=False,
    parser.add_argument('--mode','-m',default="very-fast",help="map stringent mode, very-fast (default) , fast, sensitive, OR very-sensitive")
    parser.add_argument('--clip-mode',default="end-to-end",help="end clip mode, end-to-end (default) , OR local")
    parser.add_argument('--multimap-max',type=int,default=100,help="max number of multimap align record each read threads for mapping, for bowtie(2) -k. (default:100)")
    parser.add_argument('--maxIns',type=int,default=500,help="valid for PE, max insert size. (this default: 500, bowtie2 default: 500")
    parser.add_argument('--minIns',type=int,default=10,help="valid for PE, min insert size. (this default: 10, bowtie2 default: 0")
    args = parser.parse_args()
    fastq1 = args.fastq1
    fastq2 = args.fastq2
    if not os.path.exists(args.bam_dir):
        os.makedirs(args.bam_dir)
    if not os.path.exists(args.unmapped_dir):
        os.makedirs(args.unmapped_dir)
    if not os.path.exists(args.log_dir):
        os.makedirs(args.log_dir)
    for sequence in args.priority.split(","):
        print(time.ctime()+": Start align {}".format(sequence))
        bam = os.path.join(args.bam_dir,sequence+".bam")
        unmapped1 = os.path.join(args.unmapped_dir,sequence+"_1.fastq.gz")
        unmapped2 = os.path.join(args.unmapped_dir,sequence+"_2.fastq.gz")
        unmapped = os.path.join(args.unmapped_dir,sequence+"_%.fastq.gz")
        index = os.path.join(args.index_dir,sequence)
        log = os.path.join(args.log_dir,sequence+".log")
        assert align(fastq1,fastq2,index,bam,unmapped,log,threads=args.threads,lib=args.strand,mode=args.mode,multimap=args.multimap_max,clipmode=args.clip_mode,maxIns=args.maxIns,minIns=args.minIns) == 0
        print(time.ctime()+": Done .")
        fastq1 = unmapped1
        fastq2 = unmapped2
        

if __name__ == "__main__":
    main() 
