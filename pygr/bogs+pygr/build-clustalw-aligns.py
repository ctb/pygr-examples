#! /usr/bin/env python2.5
import os

# pygr imports
from pygr import cnestedlist
import pygr.seqdb

# PTT file utilities for NCBI annotation files
import cogs2

# import clustalw utilities, too
from clustalw_utils import *

### a utility function to get only "interesting" (named, with length)
### intergenic regions from the PTT file.

def get_intergenic_intervals(ptt_info):
    intergenic = cogs2.IntergenicRegionsByFootprint(ptt_info)

    interval_dict = {}
    for (start, end, name) in intergenic:
        if '-' in name:
            continue

        key = tuple(name)
        if end - start > 0:
            interval_dict[key] = (start, end) # @CTB check end!

    return interval_dict

###

## get the two genomes; build abspath to DNA db.

thisdir = os.path.abspath(os.path.dirname(__file__))

bothdb = pygr.seqdb.BlastDB(os.path.join(thisdir, 'data/both.fna'))
ecoli_genome = bothdb['ecoliK12']
salm_genome = bothdb['salmLT2']

# load in the PTT files

ecoli_info = cogs2.CogsFileContent('data/NC_000913.ptt', len(ecoli_genome))
salm_info = cogs2.CogsFileContent('data/NC_003197.ptt', len(salm_genome))

# build dictionaries of intergenic stuff

ecoli_dict = get_intergenic_intervals(ecoli_info)
salm_dict = get_intergenic_intervals(salm_info)

# find intersection

common_keys = set(ecoli_dict.keys())
common_keys.intersection_update(salm_dict.keys())

#
# create the NLMSA object to hold the alignments.  Note that use_virtual_lpo
# will be automatically set to True (because this is a true pairwise
# alignment; see docs) but I'm specifying it *explicitly* because I want
# to remind myself to set it 
#

alignment = cnestedlist.NLMSA('pairbac', mode='w', seqDict=bothdb,
                              use_virtual_lpo=True)

# explicitly add ecoli, for some reason.  weird syntax?!
alignment += ecoli_genome

#
# iterate over all intergenic regions in common, build clustalw alignments,
# and save into NLMSA.
#

for n, key in enumerate(common_keys):
    if n % 100 == 0:
        print '...', n

    # find coords
    ecoli_start, ecoli_stop = ecoli_dict[key]
    salm_start, salm_stop = salm_dict[key]

    # extract intervals
    ecoli_ival = ecoli_genome[ecoli_start:ecoli_stop]
    salm_ival = salm_genome[salm_start:salm_stop]

    # run clustalw
    a, b = run_pair_clustalw(ecoli_ival, salm_ival)

    # build list of aligned sub-intervals
    interval_list = build_interval_list(a, b)

    # save!
    for (a, b, x, y) in interval_list:
        ec = ecoli_ival[a:b]
        sa = salm_ival[x:y]

        alignment[ec] += sa

    if n > 500:
       break

# "build" NLMSA object (this saves it to disk, too)
alignment.build(saveSeqDict=True)

# done!
