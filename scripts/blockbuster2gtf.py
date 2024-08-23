#!/usr/bin/env python
import sys
from optparse import OptionParser
# from pybedtools import BedTool

def blockbuster_to_gtf(input_file, output_file, minTxLength, maxTxLength, gnBlacklistBed=None):
    gtf_records = []
    if input_file:
        f_in = open(input_file, 'r')
    else:
        f_in = sys.stdin
        
    with f_in:
        gene_id = None
        for line in f_in:
            if line.startswith('>'):
                gene_id = line.split()[0][1:]
                continue
            fields = line.split()
            transcript_id = fields[0]
            seqname = fields[1]
            start = int(fields[2]) + 1  # Convert from 0-based to 1-based
            end = int(fields[3])
            length = end - start + 1  # Adjust length calculation
            if length < minTxLength or length > maxTxLength:
                continue
            strand = fields[4]
            score = fields[5]
            attribute = f'gene_id "{gene_id}"; transcript_id "{gene_id}_{transcript_id}"'
            gtf_line = f'{seqname}\thg38\texon\t{start}\t{end}\t{score}\t{strand}\t.\t{attribute}\n'
            # gtf_records.append(gtf_line)

    # gtf_records_str = "\n".join(gtf_records)
    # gtf = BedTool(gtf_records_str, from_string=True)

    # if gnBlacklistBed:
    #     blacklist = BedTool(gnBlacklistBed)
    #     gtf = gtf.subtract(blacklist)

            if output_file:
                with open(output_file, 'w') as f_out:
                    f_out.write(str(gtf_line))
            else:
                sys.stdout.write(str(gtf_line))

def main():
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="input_file", default=None, help="input blockbuster file")
    parser.add_option("-o", "--output", dest="output_file", default=None, help="output GTF file")
    parser.add_option("--minTxLength", dest="minTxLength", default=10, type="int", help="minimum transcript length")
    parser.add_option("--maxTxLength", dest="maxTxLength", default=200, type="int", help="maximum transcript length")
    # parser.add_option("--gnBlacklistBed", dest="gnBlacklistBed", default=None, help="genome blacklist BED file")
    (options, args) = parser.parse_args()

    # if not options.input_file or not options.output_file:
    #     parser.error("Both input and output files must be provided")

    blockbuster_to_gtf(options.input_file, options.output_file, options.minTxLength, options.maxTxLength) # , options.gnBlacklistBed

if __name__ == "__main__":
    main()