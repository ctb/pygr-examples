#! /usr/bin/env python
import sys
from pygr.seqdb import BlastDB

# create indices
db = BlastDB(sys.argv[1])

# iterate over keys & print
for name in db:
    seq = db[name]
    print '%s\t%d' % (name, len(seq))
