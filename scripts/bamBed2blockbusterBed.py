#!/usr/bin/env python
import pandas as pd
import sys
from optparse import OptionParser
#adapted from github.com/GaspardR/snakemake_blockbuster

def Dupplicate_count_column(dataf, ID_column, dup_column):
    '''
    change to more compacted format of read bed, not de-dup
    '''
    dataf = dataf[dataf.name.str.contains("/2") == False] # rm R2 for PE

    dupplicate_serie = dataf.groupby(dataf[ID_column], as_index=False).size()
    dupplicate_quants = list(dupplicate_serie.values)

    unique_keys = list(dupplicate_serie.keys()) # original
    
    dataf_dup = pd.DataFrame(
        data={
            ID_column: unique_keys,
            dup_column: dupplicate_quants
        }
    )
    dataf = pd.merge(dataf, dataf_dup, on=ID_column, how='left')
    dataf = dataf.drop_duplicates(subset=['chr', 'strand', 'start', 'end'])
    return dataf


def main(bedpath, output):
    bed_df = pd.read_csv(
        filepath_or_buffer=bedpath,
        index_col=False,
        sep='\t',
        header=None,
        names=['chr', 'start', 'end', 'name', 'score', 'strand']
    )
    ID_column = 'ID'
    dup_column = 'quantity'
    for col in ['chr', 'start', 'end', 'name', 'strand']:
        bed_df[col] = bed_df[col].map(str)
    bed_df[ID_column] = bed_df[['chr', 'start', 'end', 'strand']].apply(
        lambda x: '_'.join(x),
        axis=1
    )
    bed_df = Dupplicate_count_column(bed_df, ID_column, dup_column)
    bed_df = bed_df.sort_values(['chr', 'strand', 'start', 'end'], axis=0)
    cols = ['chr', 'start', 'end', 'name', dup_column, 'strand']
    bed_df = bed_df[cols]
    bed_df.to_csv(path_or_buf=output, index=False, sep='\t', header=None)


#main(snakemake.input.bed, snakemake.output.bed)
if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-b", "--bedpath", dest="bedpath", help="path to bed file")
    parser.add_option("-o", "--output", dest="output", help="output file path")

    (options, args) = parser.parse_args()

    if not options.bedpath or not options.output:
        parser.error("Both bedpath and output must be provided")
    print(options)
    main(options.bedpath, options.output)
