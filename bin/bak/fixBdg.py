import sys

#if len(sys.argv) == 1:
#    print "Usage: python.py <argument>"
#else:
#    use sys.argv[1]
file_path = sys.argv[1]

with open(file_path, 'r') as f:
    data = f.readline()
if len(data) == 0:
    with open(input.chrom_sizes, 'r') as f:
            chrom, size = f.readline().strip().split('\t')
            size = int(size)
    with open(file_path, 'w') as f:
        f.write('{0}\t{1}\t{2}\t0'.format(chrom, 0, size))
