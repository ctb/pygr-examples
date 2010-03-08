"""
Span objects and code to manipulate them.

Span objects represent a region in a sequence and are quite useful
for iterating over genomes, picking out regions, and extending and
comparing them.

**Classes:**

* `Span(start, end)` -- a class representing intervals from [start:end+1].

**Functions:**

* `join(list_of_spans, within_distance=0)` -- combine spans that overlap
  or are within 'within_distance' of each other.  Returns list
  of combined spans.
           
* `cover_spans(list_of_spans)` -- create a span that contains all of the
  spans in the list.

* `complement(start, end, list_of_spans)` -- create the complement of the
  spans in the list, such that the union of the complement and
  the original list will cover the interval (start, end).
"""

#
# Span
#

class Span:
    """
    A sequence span, representing an interval from [start:end+1].

    Create by calling `Span(start, end)`.
    """
    def __init__(self, start, end):
        self.start = min(start, end)
        self.end = max(start, end)

    def extend(self, point):
        """
        Extend to cover point.
        """
        if point < self.start:
            self.start = point
        elif point > self.end:
            self.end = point

    def within(self, point, dist=0):
        """
        Returns 1 if point is contained within, or is within dist (default 0),
        of this span.
        """
        if point > (self.start - dist) and point < (self.end + dist):
            return 1

        return 0

    def overlaps(self, other, within=0):
        """
        Returns 1 if the given span overlaps with this span.
        """
        assert isinstance(other, Span)
        if other.contains(self):
            return 1
        if other.start >= (self.start - within) and \
               other.start <= (self.end + within):
            return 1
        if other.end <= (self.end + within) and \
               other.end >= (self.start - within):
            return 1

        return 0

    def intersect(self, other):
        """
        Return the intersection of the two spans.
        """
        assert self.overlaps(other)
        
        start = max(self.start, other.start)
        end = min(self.end, other.end)

        return Span(start, end)

    def contains(self, other):
        """
        Returns 1 if the given span is entirely within this span.
        """
        assert isinstance(other, Span)
        if other.start >= self.start and other.end <= self.end:
            return 1
        return 0

    def join(self, other, within=0):
        """
        Returns the join of the two spans.
        """
        if not self.overlaps(other, within):
            raise Exception("spans to join must overlap within %d" % (within,))
        
        new_start = min(self.start, other.start)
        new_end = max(self.end, other.end)

        return Span(new_start, new_end)

    def __cmp__(self, other):
        if self.start == other.start:
            if self.end == other.end:
               return 0

            return cmp(self.end, other.end)
        return cmp(self.start, other.start)

    def __repr__(self):
        return "[span: %d, %d]" % (self.start, self.end,)

    def __len__(self):
        return self.end - self.start + 1

#
# join
#

def join(l, within=0):
    """
    Given a list of spans, join all of those that overlap.  Return a
    list sorted by start of span.
    """
    if len(l) == 0:
        return []

    # sort by 'begin' position.
    l.sort()

    # now, iteratively join overlapping spans.
    new_l = []
    cur = l[0]
    for i in range(1, len(l)):
        next_span = l[i]
        if cur.overlaps(next_span, within):
            cur = cur.join(next_span, within)
        else:
            new_l.append(cur)
            cur = next_span

    new_l.append(cur)

    return new_l

#
# cover_spans
#

def cover_spans(l):
    """
    Construct a single span covering all of the spans in the given list.
    """
    if not len(l):
        return None
    
    span = l[0]
    
    new_span = Span(span.start, span.end)
    for span in l[1:]:
        new_span.extend(span.start)
        new_span.extend(span.end)

    return new_span

#
# complement
#

def complement(start, end, l):
    """
    Return the complement of the spans in l from start to end.

    If you want the complement to be a strict complement, use `join`
    first.
    """
    if not len(l):
        return [Span(start, end)]

    l.sort()

    if start != l[0].start:
        complement = [Span(start, l[0].start)]
    else:
        complement = []
        
    last = l[0].end
    
    for sp in l[1:]:
        complement.append(Span(last, sp.start))
        last = sp.end

    if last != end:
        complement.append(Span(last, end))

    return complement
