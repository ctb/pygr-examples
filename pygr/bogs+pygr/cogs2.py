"""
This module contains a bunch of utilities for iterating through NCBI lists
of coordinates of genes (colloquially called "COGS files" by us).

See GenomeWrangler for usage details.
"""

import spanning

#
# NamedSpan
#

class NamedSpan(spanning.Span):
    """
    A subclass of spanning.Span that carries with it a list of names.
    """
    def __init__(self, start, end, name):
        spanning.Span.__init__(self, start, end)
        self.name = list(name)

    def join(self, other, within=0):
        new_span = spanning.Span.join(self, other, within)
        new_name = []
        new_name.extend(self.name)
        new_name.extend(other.name)
        return NamedSpan(new_span.start, new_span.end, new_name)

    def __hash__(self):
        """
        dangerous: this could conceivably change.
        """
        return hash((self.start, self.end,))

#
# CogsLine
#

class CogsLine:
    def __init__(self, start, end, strand, gene, syn, code, cog, prod):
        self.start = start
        self.end = end
        self.strand = strand.strip()
        self.gene = gene.strip()
        self.syn = syn.strip()
        self.code = code.strip()
        self.cog = cog.strip()
        self.prod = prod.strip()

#
# CogsFileContent
#

class CogsFileContent:
    def __init__(self, filename, chr_len):
        self.filename = filename
        self.chr_len = chr_len

        self._load_cogs(filename)

    def _load_cogs(self, filename):
        #
        # read in the file & find the header
        #
        lines = open(filename).readlines()
        i = 0
        while i < 15:                       # heuristic: usually < 15
            if lines[i].find('Location') >= 0:
                break

            i += 1

        assert i < 15
        line_n = i + 1

        #
        # ok, now run through the file & get all of the gene entries.
        #

        l = []
        cogs_by_syn = {}
        while line_n < len(lines):
            # parse each line

            try:
                (loc, strand, length, pid, gene, syn, code, cog, prod) = \
                      lines[line_n].split('\t', 9)
            except ValueError:
                print lines[line_n]
                print ":::".join(lines[line_n].split('\t', 9))
                raise

            # convert coords into ints
            (start, junk, end) = loc.split('.')
            (start, end) = (int(start), int(end))

            # store
            cogs_line = CogsLine(start, end, strand, gene,syn, code, cog, prod)
            l.append(cogs_line)

            assert cogs_by_syn.get(syn) is None
            cogs_by_syn[syn] = cogs_line

            line_n += 1

        self._cogs_lines = l
        self._cogs_by_syn = cogs_by_syn

    def get_cogs_lines(self):
        return self._cogs_lines

    def get_syn(self, syn):
        return self._cogs_by_syn.get(syn)

##############

class _RegionsList:
    """
    Parent class of the various region list containers.
    """
    def __iter__(self):
        return _Regions_Iterator(self)

class _Regions_Iterator:
    """
    An iterator over a list of regions; returns tuples (start, end, name).
    """
    def __init__(self, parent):
        self.parent = parent
        self.index = 0
        self.length = len(self.parent.spans)

    def next(self):
        i = self.index
        if i >= self.length:
            raise StopIteration("end of spanningregions")

        self.index += 1
        span = self.parent.spans[i]
        return (span.start, span.end, span.name)

#
# CodingRegionsByGene
#

class CodingRegionsByGene(_RegionsList):
    """
    Iterate over a list of regions, one for each gene.
    """
    def __init__(self, cogs_file_content):
        assert isinstance(cogs_file_content, CogsFileContent)
        self.cogs_content = cogs_file_content

        spans = [ NamedSpan(a.start - 1, a.end - 1, (a.gene,)) for \
                  a in self.cogs_content.get_cogs_lines() ]
        self.spans = spans

##############

#
# CodingRegionsByFootprint
#

class CodingRegionsByFootprint(_RegionsList):
    """
    Iterate over a list of regions, one for each continuous set of
    coding sequences.
    """
    def __init__(self, cogs_file_content):
        assert isinstance(cogs_file_content, CogsFileContent)
        self.cogs_content = cogs_file_content

        spans = [ NamedSpan(a.start - 1, a.end - 1, (a.gene,)) for \
                  a in self.cogs_content.get_cogs_lines() ]
        spans = spanning.join(spans)
        self.spans = spans

##############

class IntergenicRegionsByFootprint(_RegionsList):
    """
    Iterate over a list of regions, one for each intergenic region.
    """
    def __init__(self, cogs_file_content, overlap=0):
        assert isinstance(cogs_file_content, CogsFileContent)
        self.cogs_content = cogs_file_content
        self.overlap = overlap

        chr_len = cogs_file_content.chr_len

        # first get the coding spans & join 'em
        coding_spans = [ NamedSpan(a.start - 1, a.end, (a.gene,)) for \
                  a in self.cogs_content.get_cogs_lines() ]
        coding_spans = spanning.join(coding_spans)

        # then run through and construct the list of complemented spans.
        spans = []
        prev_span = coding_spans[0]

        # for the first span, start at 0 & go to the beginning of the
        # first gene.
        
        new_name = ["",]
        new_name.extend(prev_span.name)
        spans.append(NamedSpan(0, prev_span.start + overlap - 1, new_name))

        # for each successive region, take the end of the previous gene
        # and the start of the next & make a span out of it.  account
        # for overlaps...

        for i in range(1, len(coding_spans)):
            span = coding_spans[i]
            new_name = []
            new_name.extend(prev_span.name)
            new_name.extend(span.name)
            
            spans.append(NamedSpan(prev_span.end - overlap + 1,
                                   span.start + overlap - 1,
                                   new_name))
            prev_span = span

        new_name = []
        new_name.extend(prev_span.name)
        new_name.append("")
        spans.append(NamedSpan(prev_span.end - overlap + 1, chr_len, new_name))
        self.spans = spans

##############

#
# PutativeOperons
#

class PutativeOperons(_RegionsList):
    """
    Iterate over a list of putative operons.
    """

    def __init__(self, cogs_file_content, overlap=1):
        assert isinstance(cogs_file_content, CogsFileContent)
        self.cogs_content = cogs_file_content

        cogs_lines = cogs_file_content.get_cogs_lines()

        def make_name(cogs_line):
            if a.gene != '-':
                return "%s (%s)" % (a.gene, a.syn,)
            else:
                return a.syn
            
        spans = []

        a = cogs_lines[0]
        last_span = NamedSpan(a.start - 1, a.end, (make_name(a),))
        last_orient = a.strand

        orients = {}
        
        for a in cogs_lines[1:]:
            span = NamedSpan(a.start - 1, a.end, (make_name(a.gene),))
            if a.strand == last_orient and last_span.overlaps(span, overlap):
                last_span = last_span.join(span, overlap)
            else:
                spans.append(last_span)
                orients[last_span] = last_orient
                
                last_span = span
                last_orient = a.strand

            ## either way, we've consumed last_span.

        spans.append(last_span)
        self.spans = spans
        self.orients = orients

    def get_orientation(self, span):
        return self.orients.get(span)

##############

#
# PutativePromoters
#

class PutativePromoters(_RegionsList):
    """
    Iterate over a list of putative promoters.
    """

    def __init__(self, cogs_file_content, chr_len,
                 coding_overlap=1, nc_overlap=0):
        assert isinstance(cogs_file_content, CogsFileContent)
        self.cogs_content = cogs_file_content

        cogs_lines = cogs_file_content.get_cogs_lines()

        #
        # first, construct a list of operons + directions.
        #

        def make_name(cogs_line):
            if a.gene != '-':
                return "%s (%s)" % (a.gene, a.syn,)
            else:
                return a.syn

        operons = []

        a = cogs_lines[0]
        last_span = NamedSpan(a.start - 1, a.end - 1, (make_name(a),))
        last_orient = a.strand
        
        for a in cogs_lines[1:]:
            name = make_name(a)
            span = NamedSpan(a.start - 1, a.end - 1, (name,))
            if a.strand == last_orient and last_span.overlaps(span,
                                                              coding_overlap):
                last_span = last_span.join(span, coding_overlap)
            else:
                operons.append((last_span, last_orient))

                last_span = span
                last_orient = a.strand

        operons.append((last_span, last_orient))
        

        #
        # now, for each operon, assign an intergenic region (if possible)
        #

        last_end = 0
        promoters = []
        next_name = None

        for (operon, direction) in operons:
            next_start = operon.start

            # the next intergenic region belongs to the last operon...
            if next_name:
                promoter = NamedSpan(last_end - nc_overlap + 1,
                                     next_start + nc_overlap - 1, next_name)
                promoters.append(promoter)
                next_name = None

            # does this intergenic region belong to this operon?
            if direction == '+':
                promoter = NamedSpan(last_end - nc_overlap + 1,
                                     next_start + nc_overlap - 1, operon.name)
                promoters.append(promoter)
                next_name = None
            elif direction == '-':      # nope, the next.
                next_name = operon.name

            last_end = operon.end

        # redux
        if next_name:
            promoter = NamedSpan(last_end - nc_overlap + 1,
                                 next_start + nc_overlap - 1, next_name)
            promoters.append(promoter)
            
        self.spans = promoters
