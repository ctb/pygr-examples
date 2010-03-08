from pygr import sequence, cnestedlist, seqdb

def extract_closest(features, position):
    sort_by_dist = lambda x: min(abs(x.sequence.start - position),
                                 abs(x.sequence.stop - 1 - position))

    features = list(features)
    features.sort(key=sort_by_dist)
    return features[0]

def find_nearest_feature(map, sequence, position, span=1024, factor=10):
    seqlen = len(sequence)

    features = map[sequence[position:position+1]]
    if features:
        features = [ f.pathForward for f in features ]
        features = [ (len(f), f) for f in features ]
        features.sort()
        return features[0][1]
        
    while 1:
        start = max(position - span, 0)
        stop = min(position + span, seqlen)
        ival = sequence[start:stop]
        features = map[ival]

        if len(features):
            break
        
        if start == 0 and stop == seqlen:
            break
        
        span = span*factor

    if not len(features):
        return None

    return extract_closest(features, position)

def test():

    NAME='test'
    SPACING=5000

    class Spot:
        def __init__(self, name, start, stop):
            self.id = NAME
            self.name = name
            self.start = int(start)
            self.stop = int(stop)


    seq_dict = {}
    seq_dict[NAME] = s = sequence.Sequence('A'*5000000, NAME)

    annot_d = {}
    for i in range(0, len(s), SPACING):
        spot = Spot(str(i), i, i+100)
        annot_d[str(i)] = spot

    annot_db = seqdb.AnnotationDB(annot_d, seq_dict)
    annot_map = cnestedlist.NLMSA('spots', mode='memory', pairwiseMode=True)

    for v in annot_db.values():
        annot_map.addAnnotation(v)

    annot_map.build()

    print 'Built map, hoo-ha'

    assert find_nearest_feature(annot_map, s, 0).sequence.start == 0
    assert find_nearest_feature(annot_map, s, 5).sequence.start == 0
    assert find_nearest_feature(annot_map, s, 100).sequence.start == 0

    assert find_nearest_feature(annot_map, s, 2000).sequence.start == 0
    assert find_nearest_feature(annot_map, s, 2549).sequence.start == 0
    assert find_nearest_feature(annot_map, s, 2550).sequence.start == 5000
    assert find_nearest_feature(annot_map, s, 2651).sequence.start == 5000
    assert find_nearest_feature(annot_map, s, 2600).sequence.start == 5000

if __name__ == '__main__':
    test()
