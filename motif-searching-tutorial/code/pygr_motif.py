from pygr import seqdb, cnestedlist

class SiteMatch:
    def __init__(self, name, id, start, stop, orientation):
        self.name = name
        self.id = id
        if orientation == -1:
            self.start, self.stop = -stop, -start
        else:
            self.start, self.stop = start, stop

def map_matches(genome, region, sites, prefix=''):
    d = {}
    for n, (start, stop, orientation, site) in enumerate(sites):
        n = 'prefix' + str(n)
        o = SiteMatch(n, region.id, start, stop, orientation)
        d[n] = o

    annodb = seqdb.AnnotationDB(d, genome)
    map = cnestedlist.NLMSA('', 'memory', pairwiseMode=True)

    for k in annodb:
        map.addAnnotation(annodb[k])
        
    map.build()

    return map
