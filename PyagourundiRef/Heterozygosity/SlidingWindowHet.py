#From https://github.com/jarobin/vaquitagenomics2022/blob/main/genome-wide_diversity/SlidingWindowHet.py
#minimal modifications

import sys
import pysam
import os
import gzip

# Open input file and make sure the VCF file is indexed (if not, create index)
filename = sys.argv[1]
VCF = gzip.open(filename, 'rt')  # Use 'rt' for text mode in Python 3

if not os.path.exists("%s.tbi" % filename):
    pysam.tabix_index(filename, preset="vcf")
parsevcf = pysam.Tabixfile(filename)

# Set variables
chrom_lengths = sys.argv[2]
window_size = int(sys.argv[3])
step_size = int(sys.argv[4])
chrom = sys.argv[5]

# Generate a dictionary with chromosomes and chromosome lengths
cc = open(chrom_lengths, 'r')
chrom_size = {line.strip().split('\t')[0]: line.strip().split('\t')[1] for line in cc}
cc.close()

# Get list of samples from VCF file header
samples = []
for line in VCF:
    if line.startswith('##'):
        pass
    else:
        for i in line.split()[9:]:
            samples.append(i)
        break

# Get start and end positions of chromosome
VCF.seek(0)  # Reset file pointer to the beginning
for line in VCF:
    if line[0] != '#':
        start_pos = int(line.strip().split()[1])
        end_pos = int(chrom_size[chrom])
        break

# Create output file
output = open(filename + '_het_%skbWin_%sKbstep.txt' % (int(window_size/1000), int(step_size/1000)), 'w')
output.write('chrom\twindow_start\tsites_total\tcalls_%s\thets_%s\n' % ('\tcalls_'.join(samples), '\thets_'.join(samples)))


# Fetch a region, ignore sites that fail filters, tally genotype calls and heterozygotes
def snp_cal(chrom, window_start, window_end):
   #print("%s:%s" % (chrom, window_start))
    rows = tuple(parsevcf.fetch(region="%s:%s-%s" % (chrom, window_start, window_end), parser=pysam.asTuple()))
    sites_total = 0
    calls = [0] * len(samples)
    hets = [0] * len(samples)
    for line in rows:
        if line[6] != "PASS":
            continue
        sites_total += 1
        for i in range(0, len(samples)):
            if line[i + 9][:1] == '.':
                continue
            calls[i] += 1
            GT = line[i + 9].split(':')[0]
            if '/' in GT:
                sp = '/'
            if '|' in GT:
                sp = '|'
            if GT.split(sp)[0] != GT.split(sp)[1]:
                hets[i] += 1
    output.write('%s\t%s\t%s\t%s\t%s\n' % (chrom, window_start, sites_total, '\t'.join(map(str, calls)),
                                             '\t'.join(map(str, hets))))


# Initialize window start and end coordinates
window_start = start_pos
window_end = start_pos + window_size - 1

# Calculate stats for window, update window start and end positions,
# repeat to the end of the chromosome
while window_end <= end_pos:
    if window_end < end_pos:
        snp_cal(chrom, window_start, window_end)
        window_start = window_start + step_size
        window_end = window_start + window_size - 1
    else:
        snp_cal(chrom, window_start, window_end)
        break
else:
    window_end = end_pos
    snp_cal(chrom, window_start, window_end)

# Close files and exit
VCF.close()
output.close()

exit()

