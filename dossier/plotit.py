import numpy as np
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from csv import DictReader
import math

def read_results(gname):
    return rows

def read_geography(gene):
    inf = file( 'geography%s.txt' % gene)
    transcript = inf.readline()
    strand = inf.readline()
    for i in range(3):
        l = inf.readline() #crap
    cdsline = inf.readline()
    x = cdsline.strip().split('\t')
    cstart = int(x[0])
    cend = int(x[1])
    tstart = 1
    for line in inf:
        pass #just go to the end one.
    x = line.strip().split('\t')
    tend = int(x[4])
    inf.close()    
    return tstart,tend,cstart,cend

def read_hotspots(gname):
    fname = 'hotspots_%s.txt' % gname
    inf = file(fname,'r')
    rows = DictReader(inf,dialect='excel-tab')
    rlist = []
    for row in rows:
        rlist.append(row)
    inf.close()
    return rlist

def read_regions(gname):
    fname = 'regions%s.txt' % gname
    inf = file(fname,'r')
    rows = DictReader(inf,dialect='excel-tab')
    regions = []
    for row in rows:
        if row['start_type'].startswith('#'):
            continue
        re = int(row['region_end'])
        rs = int(row['region_start'])
        points = (re,rs)
        note = row['note']
        name = row['feature_type_type_name']
        regions.append( ( min(points), max(points), note, name ) )
    inf.close()
    return regions
    
ve2col = {'synonymous': 'yellow', 'missense': 'k', 'nonsense': 'r', 'frameshifting indel': 'r', \
          'non-frameshifting indel': 'orange', 'UTR-5': 'yellow', 'UTR-3': 'yellow'}

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

def drawgene(ts,te,cs,ce,gene,regions,hspots):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    gheight = 0
    plt.plot( (ts,cs), (gheight,gheight), color = 'blue' )
    plt.plot( (ce,te), (gheight,gheight), color = 'blue' )
    plt.plot( (cs,ce), (gheight,gheight), color = 'blue' , lw=4)
    #i = 0
    for r in regions:
        plt.plot( (r[0], r[1]), (gheight-0.1, gheight-0.1) , lw= 2, color = 'g')
        #plt.plot( (r[0], r[1]), (gheight-0.005*i, gheight-0.005*i) , lw= 2, color = 'g')
        #i += 1
    h_height = 0.2
    h_factor = 0.12
    for h in hspots:
        pv = float(h['p_value'])
        n = float(h['N_perm'])
        hstart = int(h['start'])
        hend = int(h['end'])
        if pv != 0:
            y = h_height + h_factor * math.log10(pv)
            plt.plot( (hstart, hend), (y,y) , lw= 2, color = 'b')
        else:
            pv = 1./(n+1)
            y = h_height + h_factor * math.log10(pv)
            plt.plot( (hstart, hend), (y,y) , lw= 2, color = 'r')
    fname = 'results_%s.txt' % gene
    inf = file(fname,'r')
    rows = DictReader(inf,dialect='excel-tab')
    badheight = -0.05
    hgmdheight = -0.1
    for row in rows:
        hgvs = row['hgvstranscript']
#        print 'hgvs is %s' %hgvs
        if hgvs != '' and hgvs != '?' and hgvs != 'None':
            tpos = get_tpos(hgvs)
            ve = row['variant_effect']
            tag = row['tag']
            if ve == 'nonsense':
                plt.plot( (tpos,), (badheight,), 'o', color = 'orange' )
            elif ve == 'frameshifting indel':
                plt.plot( (tpos,), (badheight,), 'x', color = 'red' )
            elif ve == 'missense' and tag == 'DM':
                plt.plot( (tpos,), (hgmdheight,), '+', color = 'purple' )
    ax.set_title(gene)
    ax.set_ylim(-0.5,1)

    #ax.annotate("Outside %s,\nInside %s" % (wise(outside), wise(inside)),
    #            (i * 2.5, -1.5), va="top", ha="center")

    #ax.set_xlim(-2,10)
    #ax.set_title('Mmm, donuts!')
    #ax.set_aspect(1.0)
    plt.show()



def go(gname):
    tstart,tend,cstart,cend = read_geography(gname)
    features = read_regions(gname)
    hotties = read_hotspots(gname)
    drawgene(tstart,tend,cstart,cend,gname,features,hotties)

if __name__ == '__main__':
    go('GRN')
