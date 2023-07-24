#! /usr/bin/env python

# todo: not finished !!! 
# to replace count_matrix_* rules in peak_common.smk

import re
bin_size = 5
decay = 50 
pvalue_expeak = 1
sample_id = "test"
pat_cov = re.compile(r'/expeakCNN_counts/b(?P<bin_size>[^/]+)_d(?P<decay>[^/]+)_p(?P<pvalue_expeak>[^/]+)/(?P<sample_id>[^\.]+).bed')
filename = "/expeakCNN_counts/b5_d50_p1/test.bed"

command_handlers = {}
def command_handler(f):
    command_handlers[f.__name__] = f
    return f

@command_handler
def count_peaks_localmax(args):
    import pandas as pd
    import re
    import numpy as np
    input_transcript_table = ""
    input_chrom_sizes = ""
    input_domains = ""
    domain_dir = ""
    domain_dir_sub = "domain_localmax_counts_EM2"
    input_peaks = ""
    
    transcript_table = pd.read_table(input_transcript_table, sep='\t', dtype='str')
    transcript_table.drop_duplicates(['transcript_id'], inplace=True)
    transcript_table.set_index('transcript_id', drop=False, inplace=True)
    transcript_table = transcript_table.loc[:, ['gene_id', 'gene_name', 'gene_type', 'transcript_id', 'transcript_type', 'start', 'end']].copy()
    # extend transcript_table with genome regions
    chrom_sizes = pd.read_table(input_chrom_sizes, sep='\t', names=['chrom', 'end'])
    chrom_sizes.set_index('chrom', drop=False, inplace=True)
    domains = pd.read_table(input_domains, sep='\t', header=None,
        names=['chrom', 'start', 'end', 'domain_id', 'score', 'strand'], dtype='str')

    pat_cov = re.compile(r'{domain_dir}/{domain_dir_sub}/b(?P<bin_size>[^/]+)_d(?P<decay>[^/]+)_p(?P<pvalue>[^/]+)/(?P<sample_id>[^\.]+).bed'.format(domain_dir=domain_dir,domain_dir_sub=domain_dir_sub))
    # pat_cov = re.compile(r'{domain_dir}/domain_localmax_counts_EM/(?P<sample_id>[^\.]+).bed'.format(domain_dir=domain_dir))
    mat = []
    peak_labels = None
    for filename in input_peaks:
        sample_id = pat_cov.match(filename).groupdict()['sample_id']
        df = pd.read_table(filename, header=None)
        if peak_labels is None:
            peak_labels = df.iloc[:, 3].values
        df.index = df.iloc[:, 3]
        cov = df.iloc[:, 4].copy()
        cov.name = sample_id
        mat.append(cov)
    mat = pd.concat(mat, axis=1)
    # get seq
    seq_ids = domains['chrom'].values
    # get transcript peaks
    is_genome_peaks = np.isin(seq_ids, chrom_sizes['chrom'].values)
    seq_ids_genome = seq_ids[is_genome_peaks]
    seq_ids_transcript = seq_ids[~is_genome_peaks]
    # annotate transcript peaks with gene information
    # feature name format: gene_id|gene_type|gene_name|domain_id|transcript_id|start|end
    feature_names = np.empty(mat.shape[0], dtype='object')
    print(np.sum(~is_genome_peaks), seq_ids_transcript.shape, transcript_table.loc[seq_ids_transcript, 'gene_name'].values.shape)
    feature_names[~is_genome_peaks] = transcript_table.loc[seq_ids_transcript, 'gene_id'].values \
        + '|' + transcript_table.loc[seq_ids_transcript, 'gene_type'].values \
        + '|' + transcript_table.loc[seq_ids_transcript, 'gene_name'].values \
        + '|' + domains['domain_id'].values[~is_genome_peaks] \
        + '|' + transcript_table.loc[seq_ids_transcript, 'transcript_id'].values \
        + '|' + domains['start'].values[~is_genome_peaks] \
        + '|' + domains['end'].values[~is_genome_peaks]
    # annotate genome peaks
    print(seq_ids_genome.shape, np.sum(is_genome_peaks))
    gene_ids_genome = seq_ids_genome + '_' + domains['start'].values[is_genome_peaks] \
        + '_' + domains['end'].values[is_genome_peaks] + '_' + domains['strand'].values[is_genome_peaks]
    feature_names[is_genome_peaks] = gene_ids_genome \
        + '|' + transcript_table.loc[seq_ids_genome, 'transcript_type'].values \
        + '|' + transcript_table.loc[seq_ids_genome, 'transcript_id'].values \
        + '|' + domains['domain_id'].values[is_genome_peaks] \
        + '|' + domains['chrom'].values[is_genome_peaks] \
        + '|' + domains['start'].values[is_genome_peaks] \
        + '|' + domains['end'].values[is_genome_peaks]
    mat.index = feature_names
    mat.index.name = 'feature'
    mat.to_csv(output[0], sep='\t', header=True, index=True)



if __name__ == '__main__':
    main_parser = argparse.ArgumentParser(description='Count peaks from exRNA tbed.gz')
    subparsers = main_parser.add_subparsers(dest='command')
    
    parser = subparsers.add_parser('count_peaks_localmax')
    # parser.add_argument('--input-bam', type=str, required=True,
    #     help='input Bowtie2 Bam file of raw reads')
    # # parser.add_argument('--input-bw', type=str, required=True,
    # #     help='input BigWig file of raw reads coverage')
    # parser.add_argument('--min-peak-length', type=int, default=10,
    #     help='minimum length required for a peak')
    # parser.add_argument('--max-peak-length', type=int, default=200,
    #     help='maximum length required for a peak')
    # parser.add_argument('--decay', type=float, default=0.0,
    #     help='decay factor of peak summit to define peak boundary')
    # parser.add_argument('--permutate-pval', type=float, default=0.05,
    #     help='pval for bg cov permutation estimation')
    # parser.add_argument('--poisson-pval', type=float, default=0.05,
    #     help='pval for poisson dist significance test')
    # parser.add_argument('--min-cov', type=float, default=5,
    #     help='minimum coverage required to define a peak')
    # parser.add_argument('--bin-width', type=int, default=10,
    #     help='bin width to search enriched bins')
    # parser.add_argument('--max-iter', type=int, default=50,
    #     help='max iteration num to estimate bg cov')
    # parser.add_argument('--seed', type=int, default=1234,
    #     help='seed for random permutation')
    parser.add_argument('--output-file', '-o', type=str, required=True,
        help='output peaks in BED format')
    # parser.add_argument('--mode', type=str, default="local",
    #     help='call peak permutation mode when estimate bg cov: local (default) or global')
    # parser.add_argument('--process', type=int, default=1,
    #     help='num of processes to call peak in parallel, default: 1')
    # parser.add_argument('--thread', type=int, default=max(int(0.5*os.cpu_count()), 1),
    #     help='num of threads to call peak in parallel (in each process), default: half num of cpu dectected')

    args = main_parser.parse_args()
    if args.command is None:
        raise ValueError('empty command')
    # logger = logging.getLogger('call_peak.' + args.command)

    command_handlers.get(args.command)(args)
    

