#! /usr/bin/env python2.5
"""
A reimplementation of load-probes-by-gene using pygr.
"""
import sys
from pygr import seqdb, cnestedlist
import math

REGION_SIZE=10000
WINDOW_SIZE=250
PROBE_SIZE=50
SIGNAL_CUTOFF=1.25

scaffold_file, gene_starts, probe_data, output_file = sys.argv[1:5]

### load probe data

class SignalPair:
    def __init__(self, gene, seq_id, probe, position, cy3, cy5):
        self.gene = gene
        self.seq_id = seq_id
        self.probe = probe
        self.start, self.stop = position, position + PROBE_SIZE
        self.cy3 = cy3
        self.cy5 = cy5
        self.id = self.scaffold_name()
        self.name = probe
        self.logratio = math.log(cy5 / cy3, 2)

    def scaffold_name(self):
        return '_'.join(self.seq_id.split('_')[:2])

def read_signal_pairs(filename):
    d = {}
    
    first = True
    for n, line in enumerate(open(filename)):
        if first:
            first = False
            continue

        if line.startswith('RANDOM'):
            continue

        (_, gene, seq_id, probe, position, cy3, cy5) = line.split()

        signal = SignalPair(gene, seq_id, probe, int(position),
                            float(cy3), float(cy5))

        assert not probe in d
        d[probe] = signal

    return d

###

print 'reading signal pairs'
probe_dict = read_signal_pairs(probe_data)

print 'reading genome scaffolds'
scaffolds = seqdb.BlastDB(scaffold_file)

print 'constructing probe annotations'
annotations = seqdb.AnnotationDB(dict(probe_dict), scaffolds,
                                 annotationType='probe:')

probe_map = cnestedlist.NLMSA('probe_map', mode='memory',
                              use_virtual_lpo=True, bidirectional=False)

for annotation in annotations.values():
    probe_map.addAnnotation(annotation)

probe_map.build()

#probe = probe_dict.values()[0]
#test = scaffolds[probe.id]
#k = probe_map[test].keys()[0]


###

signal_bins = [ 0 ] * (2*REGION_SIZE / WINDOW_SIZE + 1)
count_bins = [ 0 ] * (2*REGION_SIZE / WINDOW_SIZE + 1)
bin_mid = REGION_SIZE/WINDOW_SIZE

print 'reading gene starts'
lines = open(gene_starts).readlines()[1:]
for line in lines:
    gene_name, gene_start, scaffold_name, orient = line.split('\t')
    gene_start = int(gene_start)
    orient = int(orient)

    region = scaffolds[scaffold_name]

    start = max(0, gene_start - REGION_SIZE)
    stop = gene_start + REGION_SIZE

    before_region = region[start:gene_start]
    after_region = region[gene_start:stop]

    if orient < 0:
        before_region, after_region = -after_region, -before_region

    for i in range(1, len(before_region), WINDOW_SIZE):
        ival = before_region[-i - WINDOW_SIZE:-i]
        try:
            signals = [ p.logratio for p in probe_map[ival].keys() ]
        except KeyError:
            continue

        distance = -(i + WINDOW_SIZE / 2)
        signal = sum(signals)
        count = len(signals)

        bin = distance / WINDOW_SIZE + bin_mid
        signal_bins[bin] += signal
        count_bins[bin] += count

    for i in range(0, len(after_region), WINDOW_SIZE):
        ival = after_region[i:i+WINDOW_SIZE]
        try:
            signals = [ p.logratio for p in probe_map[ival].keys() ]
        except KeyError:
            continue

        distance = i + WINDOW_SIZE / 2
        signal = sum(signals)
        count = len(signals)

        bin = distance / WINDOW_SIZE + bin_mid
        signal_bins[bin] += signal
        count_bins[bin] += count

fp = open(output_file, 'w')
for i in range(0, len(signal_bins)):
    if count_bins[i]:
        print >>fp, i, signal_bins[i] / float(count_bins[i])
