from csv import DictReader
from collections import defaultdict
import random

def get_tpos(hgvs):
    p = hgvs.split(':')[1].split('.')[1]
    try:
        n = 1
        while True:
            tpos = int(p[:n])
            n+=1
    except:
        pass
    return tpos

def get_hgmd_positions(gname):
    posset = set()
    inf = file('results_%s.txt' % gname,'r')
    rows = DictReader(inf,dialect='excel-tab')
    for row in rows:
        hgvs = row['hgvstranscript']
        ve = row['variant_effect']
        tag = row['tag']
        if ve == 'missense' and tag == 'DM':
            tpos = get_tpos(hgvs)
            posset.add( tpos )
    inf.close()
    positions = []
    positions.extend(posset)
    positions.sort()
    return positions

def get_cds_length(gname):
    inf = file('geography%s.txt' % gname, 'r')
    for i in range(6):
        line = inf.readline()
    inf.close()
    x = line.strip().split('\t')
    if len(x) != 2:
        print 'whoops'
        print line
        exit()
    l = int(x[1]) - int(x[0])
    return l

def generate_distances(np,positions):
    dists = []
    ranges = defaultdict(list)
    for i in range(len(positions)-(np-1)):
        p0 = positions[i]
        p1 = positions[i+np-1]
        dists.append( p1 - p0 )
        ranges[p1-p0].append( (p0,p1) )
    return dists,ranges
         
def make_permutations(npoints,interval_size,nperm,rlow,rhigh):
    distribution = defaultdict(list)
    allpos = range(interval_size)
    for perm_i in range(nperm):
        random.shuffle(allpos)
        points = allpos[:npoints]
        points.sort()
        for d in range(rlow,rhigh):
            ds,rs = generate_distances( d, points)
            distribution[d].append( ds[0] )
            distribution[d].sort()
    return distribution

def get_p_value(value,distlist):
    nless = 0
    for perm_d in distlist:
        if perm_d <= value:
            nless += 1
        else:
            break
    return 1. * nless / len(distlist)
                

def go(gname,nperm):
    print gname
    outf = file('hotspots_%s.txt' % gname,'w')
    outf.write('Npoints\tDistance\tp_value\tstart\tend\tN_perm\n')
    pos = get_hgmd_positions(gname)
    print len(pos)
    clength = get_cds_length(gname)
    rlow = 2
    rhigh = min( [100, len(pos) - 1] )
    pdists = make_permutations(len(pos),clength,nperm,rlow,rhigh)
    for i in range(rlow,rhigh):
        distances_i,ranges = generate_distances(i,pos)
        distances_i.sort()
        for dist in distances_i:
            p = get_p_value(dist,pdists[i])
            if p < 0.001:
                for r in ranges[dist]:
                    outf.write('%d\t%d\t%f\t%d\t%d\t%d\n' % ( i,dist,p,r[0],r[1],nperm ) )
    outf.close()

if __name__ == '__main__':
    #go('GRN',100000)
    go('MAPT',100000)
    #go('SOD1',10000)


