#! /usr/bin/env python2.5
import sys
from pygr import cnestedlist
import motility
import bndarray

#
# first, load in the ecoli/salm alignments.
#

# note: use_virtual_lpo was set to True on save => use to load as well.
alignment = cnestedlist.NLMSA('pairbac', 'r', use_virtual_lpo=True)

# retrieve the E. coli genome from the alignment.
ecoli_genome = alignment.seqDict['ecoliK12']

# load energy operator
op_en = bndarray.parse_as_motility_operator(open('crp_init.open'))
op_en = motility.EnergyOperator(op_en)

#
# now, search for motif matches and iterate over the results.
#

results = op_en.find(str(ecoli_genome), 7.0)

count = 0
diffs = [0] * len(op_en)
for (start, stop, orient, _) in results:

    #
    # convert each motif match into a sliced sequence suitable for
    # querying the NLMSA.
    #

    ecoli_site = ecoli_genome[start:stop]
    if orient == -1:
        ecoli_site = -ecoli_site        # reverse complement

    #
    # get an edge to the aligned salmonella sequences (if any)
    #

    edge = alignment[ecoli_site]

    #
    # now *retrieve* aligned sequences that have precisely the right size
    # (no gaps/insertions, by default, and with length equal to the query)
    #

    salm_sites = edge.keys(minAlignSize=len(ecoli_site))

    for n, salm_site in enumerate(salm_sites):
        count += 1
        for j in range(0, len(ecoli_site)):
            if str(ecoli_site)[j] != str(salm_site)[j]:
                diffs[j] += 1

# print mutation profile
for n, val in enumerate(diffs):
    print n, val / float(count)
