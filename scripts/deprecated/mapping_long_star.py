#!/usr/bin/env python
import subprocess
import argparse
import os
import io
import time
os.environ['JOBLIB_TEMP_FOLDER'] = '/BioII/lulab_b/baopengfei/tmp' # seems not working

# todo: currently STAR subprocess error
def align(fastq1,fastq2,index,bam,unmapped1,unmapped2,log,threads,output_prefix,multimap,clipmode):
    # original version: --sensitive not --very-fast
    mapping_cmd = ["/BioII/lulab_b/baopengfei/gitsoft/STAR-2.5.4a/bin/Linux_x86_64/STAR","--outReadsUnmapped","Fastx","--readFilesCommand","gzip","-d","-c","--outSAMtype", "BAM","Unsorted",
        "--genomeDir",index,"--readFilesIn",fastq1,fastq2,
        "--runThreadN",str(threads),"--outFileNamePrefix",output_prefix,
        "--alignEndsType",clipmode,"--outFilterMultimapNmax",str(multimap),"--seedPerWindowNmax",str(20),">",log,"2>&1"]
    print(" ".join(mapping_cmd))
    os.system(" ".join(list(map(str, mapping_cmd))))
    # ps = subprocess.Popen(mapping_cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    # ps.wait()
    
    # STAR tend to meet subprocess error, use os.system instead
    
    tmpBam = output_prefix+"Aligned.out.bam"
    tmpMate1 = output_prefix+"Unmapped.out.mate1"
    tmpMate2 = output_prefix+"Unmapped.out.mate2"
    tmpMate3 = output_prefix+"Unmapped.out.mate3"
    tmpMate4 = output_prefix+"Unmapped.out.mate4"
    
    sort_cmd = ["/BioII/lulab_b/jinyunfan/anaconda3/envs/bioinfo_py36/bin/samtools","sort","-@",str(threads),"-o",bam,tmpBam]
    print(" ".join(sort_cmd))
    # subprocess.run(sort_cmd)
    os.system(" ".join(list(map(str, sort_cmd))))

    index_cmd = ["/BioII/lulab_b/jinyunfan/anaconda3/envs/bioinfo_py36/bin/samtools","index","-@",str(threads),bam]
    print(" ".join(index_cmd))
    os.system(" ".join(list(map(str, index_cmd))))
    # subprocess.run(index_cmd)

    #https://linuxhint.com/python-subprocess-pipes/
    logrepair = log+"repair"
    repair_cmd = ["/BioII/lulab_b/baopengfei/biosoft/bbmap/repair.sh","overwrite=t","in=",tmpMate1,"in2=",tmpMate2,"out=",tmpMate3,"out2=",tmpMate4,">>",logrepair,"2>&1"] 
    print(" ".join(repair_cmd))
    os.system(" ".join(repair_cmd))
    # subprocess.run(repair_cmd)

    pgz1_cmd = ["pigz","-c","-p",str(threads),tmpMate3,">",unmapped1]
    print(" ".join(pgz1_cmd))
    os.system(" ".join(pgz1_cmd))
    # subprocess.run(pgz1_cmd)
    pgz2_cmd = ["pigz","-c","-p",str(threads),tmpMate4,">",unmapped2]
    print(" ".join(pgz2_cmd))
    os.system(" ".join(pgz2_cmd))
    # subprocess.run(pgz2_cmd)
    
    os.system("rm -f "+tmpBam+" "+tmpMate1+" "+tmpMate2+" "+tmpMate3+" "+tmpMate4)

    # with open(log,"w") as f:
    #     for line in io.TextIOWrapper(ps.stderr, encoding="utf-8"):
    #         f.write(line)
    # return ps.poll() 
    return (not (os.path.exists(unmapped2) & os.path.exists(bam+".bai")))  # return 0 mean success


#spikein, univec, rRNA, lncRNA, miRNA, mRNA, piRNA, snoRNA, snRNA, srpRNA, tRNA, tucpRNA, Y RNA, circRNA, genome
def main():
    parser = argparse.ArgumentParser(description='Sequential Tx Map/Align by STAR')
    parser.add_argument('--fastq1','-f1',required=True,help="Input paired end cleaned fastq1")
    parser.add_argument('--fastq2','-f2',required=True,help="Input paired end cleaned fastq2")
    parser.add_argument('--bam-dir','-bd',required=True,help="Dir for output bam files")
    parser.add_argument('--unmapped-dir','-fd',required=True,help="Dir for output unmapped fastq files")
    parser.add_argument('--index-dir','-id',default="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/index/star",
                        help="Dir contains STAR index")
    parser.add_argument('--log-dir','-ld',required=True,help="Dir for output log files")
    parser.add_argument('--priority','-p',default="spikein_small,univec,rRNA,miRNA,lncRNA,mRNA,piRNA,snoRNA,snRNA,srpRNA,tRNA,tucpRNA,Y_RNA,intron.for,intron.rev,promoter.for,promoter.rev,enhancer.for,enhancer.rev",help="Priority of mapping")
    parser.add_argument('--threads','-t',type=int,default=1,help="Number of threads for mapping")
    # parser.add_argument('--mode','-m',default="very-fast",help="map stringent mode, very-fast (default) , fast, sensitive, OR very-sensitive")
    parser.add_argument('--clip-mode',default="EndToEnd",help="end clip mode for STAR --alignEndsType, (this default: EndToEnd, STAR default: Local)")
    parser.add_argument('--multimap-max',type=int,default=100,help="max number of multimap align record each read threads for mapping for STAR --outFilterMultimapNmax. (default:100)")

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
        index = os.path.join(args.index_dir,sequence)
        log = os.path.join(args.log_dir,sequence+".log")
        output_prefix_dir = os.path.join(args.log_dir,sequence)
        if not os.path.exists(output_prefix_dir):
            os.makedirs(output_prefix_dir)
        output_prefix = output_prefix_dir +"/"
        assert align(fastq1,fastq2,index,bam,unmapped1,unmapped2,log,threads=args.threads,output_prefix=output_prefix,multimap=args.multimap_max,clipmode=args.clip_mode) == 0
        print(time.ctime()+": Done .")
        fastq1 = unmapped1
        fastq2 = unmapped2
        

if __name__ == "__main__":
    main() 




