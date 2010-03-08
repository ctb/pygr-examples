import subprocess

def read_pair_clustalw(fp):
    """
    Read aligned sequences from a pairwise CLUSTALW alignment.
    """
    lines = fp.readlines()
    assert lines[0].startswith('CLUSTAL '), lines[0]

    lines = lines[3:]
    a = ''
    b = ''
    for i in range(0, len(lines), 4):
        x, y, = lines[i:i+2]
        a += x[16:].strip()
        b += y[16:].strip()

    return a, b

def run_pair_clustalw(top, bot):
    """
    Run a pairwise CLUSTALW alignment on the two sequences & returned
    aligned sequences.
    """
    fp = open('out', 'w')
    print >>fp, '>ecoli'
    print >>fp, top
    print >>fp, '>salm'
    print >>fp, bot
    fp.close()

    subprocess.check_call('clustalw out', shell=True, stdout=subprocess.PIPE)

    a, b = read_pair_clustalw(open('out.aln'))
    return a, b

def build_interval_list(a, b):
    """
    Hacky code to extract all ungapped aligned subintervals from a
    pair of aligned sequences.
    """
    interval_list = []

    a_start = None
    b_start = None

    a_count = b_count = 0
    for i in range(0, len(a)):
        if a[i] == '-' or b[i] == '-':
            if a_start is not None:           # want to end at i-1
                interval_list.append((a_start, a_count, b_start, b_count))

                a_start = b_start = None
        else:
            if a_start is None:
                a_start = a_count
                b_start = b_count

        if a[i] != '-':
            a_count += 1
        if b[i] != '-':
            b_count += 1

    if a_start is not None:
        interval_list.append((a_start, a_count, b_start, b_count))

    assert a_count == len(a.replace('-', ''))
    assert b_count == len(b.replace('-', ''))
        
    return interval_list

