import re
import sys
#from functools import partial
import pysam



def add_umi_tag(in_bam, out_bam, tag="RX", delim="_", frag=0):
    # cdef:
    #     AlignmentFile inbam, outbam
    #     AlignedSegment aln
    #     int aln_count
    #     str id, umi

    # print('Parsing from %s to %s' %(in_bam, out_bam), file = sys.stderr)
    with pysam.Samfile(in_bam, 'rb') as inbam:
        with pysam.Samfile(out_bam,'wb', template=inbam) as outbam:
            for aln_count, aln in enumerate(inbam):
                id = aln.query_name
                splitted_id = id.split(delim)
                umi = splitted_id[frag]
                aln.query_name = id.replace(umi + delim, '').replace(delim+umi, '')
                aln.tags += [(tag,umi)]
                outbam.write(aln)
    # print('Parsed %i alignments from %s to %s' %(aln_count, in_bam, out_bam), file = sys.stderr)
    return 0

def main():
    add_umi_tag(in_bam=sys.argv[1], out_bam=sys.argv[2])

if __name__ == '__main__':
    main()