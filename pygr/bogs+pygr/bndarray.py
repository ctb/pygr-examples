#! /usr/bin/env python
"""
This is a library for parsing binding arrays from files.  See the 'motility'
project for more information; this is a simple & slow standalone utility
for dealing with similar stuff in a particular format.
"""
import sys
from biolib import fasta

#
# BindingMatrix
#

class BindingMatrix:
    """
    A class to hold a representation of a binding matrix.
    """
    def __init__(self, length, arr):
        "Initialize."
        self.length = int(length)
        assert len(arr) == self.length
        self.arr = arr

    def match_energy(self, site):
        "Calculate the match of the site under this matrix."
        
        assert len(site) == self.length

        conv = { 'A' : 0, 'C' : 1, 'G' : 2, 'T' : 3 }

        total = 0.
        for i in range(0, self.length):
            base = site[i]
            base_i = conv[base]
            strength = (self.arr[i])[base_i]
            total += strength

        return total

    def get_as_motility_operator(self):
        """
        Tack on a 1000 for an 'N' match.  The resulting operator
        can be used directly by the _motility functions.
        """
        
        operator = []
        for (a,c,g,t) in self.arr:
            operator.append((a, c, g, t, 1000.,))

        return operator

#
# parse_bndarray
#

def parse_bndarray(fp):
    """
    A function to load a bndarray from a file.  Returns an object of
    type BindingMatrix.
    """
    length = int(fp.readline())
    
    arr = []
    for l in fp:
        (A, C, G, T) = l.split()
        (A, C, G, T) = map(float, (A, C, G, T))

        arr.append((A, C, G, T,))

    return BindingMatrix(length, arr)

#
# parse_as_motility_operator
#

def parse_as_motility_operator(fp):
    """
    A function to load a bndarray from a file.  Returns a (5xN) array
    that can be fed directly into motility.EnergyOperator.
    """
    length = int(fp.readline())
    
    arr = []
    for l in fp:
        (A, C, G, T) = l.split()
        (A, C, G, T) = map(float, (A, C, G, T))

        arr.append((A, C, G, T, 1000.))

    return arr
