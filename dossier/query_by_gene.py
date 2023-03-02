from collections import defaultdict
import datetime
import logging
import psycopg2
import signal
from psycopg2.extras import RealDictCursor

#this will include both exons and introns
gd1='''CREATE TABLE bizon.gene_dossier_one AS
(
SELECT * 
 FROM refseq.variants_48_2 rsv
 WHERE  rsv.hgnc_gene = %s 
   AND  rsv.transcr = %s
)
'''

gd2 = '''CREATE TABLE bizon.gene_dossier_freq AS
( 
 select * from clinbin.max_freq where loc_var_id in ( 
  SELECT loc_var_id from bizon.gene_dossier_one
 )
)
'''

gd3 = '''CREATE TABLE bizon.gene_dossier_hgmd AS 
(
  SELECT one.loc_var_id, h.acc_num, h.tag 
  FROM bizon.gene_dossier_one one
  LEFT JOIN hgmd.hgmd_loc_var h 
  ON h.loc_var_id = one.loc_var_id)'''


gd4 = '''SELECT * 
FROM bizon.gene_dossier_one one,
     bizon.gene_dossier_freq freq,
     bizon.gene_dossier_hgmd hgmd
WHERE one.loc_var_id = freq.loc_var_id
  AND one.loc_var_id = hgmd.loc_var_id
ORDER BY pos;
'''

regionquery='''SELECT * 
FROM refseq.region_group rg, 
     refseq.region_group_regions rgr, 
     refseq.feature f 
WHERE rg.transcr_ver_id = %s 
  AND rg.region_group_id = f.loc_region_group_id 
  AND rg.region_group_id = rgr.region_group_id'''

def write_results(fname,results):
    of = file(fname,'w')
    res = results.fetchone()
    print res
    keys = [ x for x in res ]
    of.write('\t'.join(keys))
    of.write('\n')
    while res != None:
        ol = [ str(res[x]) for x in keys ]
        of.write('\t'.join(ol))
        of.write('\n')
        res = results.fetchone()
    of.close()

def get_geography(genename,transcript):
    of = file('geography%s.txt' % genename,'w')
    inf = file('../bin_runner/details.txt')
    line = inf.readline()
    print transcript
    while line != '' and not line.startswith(transcript):
        line = inf.readline()
    if line == '':
        print 'crap'
        exit()
    while len(line.strip()) > 0:
        of.write(line)
        line = inf.readline()
    inf.close()
    of.close()

def run_by_gene(genename,transcript,db,logger):
    logger.debug( 'drop tables' )
    otime = datetime.datetime.now()
    db.runNoResultSql('DROP TABLE IF EXISTS bizon.gene_dossier_one '  ,())
    db.runNoResultSql('DROP TABLE IF EXISTS bizon.gene_dossier_freq '  ,())
    db.runNoResultSql('DROP TABLE IF EXISTS bizon.gene_dossier_hgmd '  ,())
    dtime = datetime.datetime.now()
    logger.debug( 'to delete tables: %d'% (dtime - otime).seconds)

    db.runNoResultSql(gd1,(genename,transcript))
    stime = datetime.datetime.now()
    logger.debug( 'to build first: %d'% (stime - dtime).seconds)

    db.runNoResultSql(gd2,())
    ttime = datetime.datetime.now()
    logger.debug( 'to build second: %d'  % (ttime - stime).seconds)

    db.runNoResultSql(gd3,())
    utime = datetime.datetime.now()
    logger.debug( 'to build third: %d'  % (utime - ttime).seconds)

    results = db.runselect(gd4,())
    write_results('results_%s.txt' % genename, results)

    regres = db.runselect(regionquery,(transcript,))
    write_results('regions%s.txt' % genename, regres)

    get_geography(genename,transcript)

class mydb:
    def __init__(self,logger):
        self.logger = logger
        dbname = 'vardb_berg'
        ref_user = 'bizon'
        db_ip = 'genomicsdb.renci.unc.edu'
        pw = 'MOugRibLfLs4'
        cstring = 'dbname=%s user=%s host=%s password=%s' % \
            (dbname, ref_user, db_ip, pw)
        self.conn = psycopg2.connect(cstring)
        signal.signal(signal.SIGTERM, self.sig_handler)
        signal.signal(signal.SIGINT, self.sig_handler)
    def sig_handler(self,sig,frm):
        self.logger.info('SIGNAL %s in JobDB. Canceling Connection' % sig)
        self.conn.cancel()
        self.logger.info('Canceled Connection.')
    def runselect(self,statement, vals):
        cur = self.conn.cursor(cursor_factory = RealDictCursor)
        cur.execute(statement, vals)
        return cur
    def runNoResultSql(self,statement, vals):
        cur = self.conn.cursor()
        cur.execute(statement, vals)
        self.conn.commit()

def test():
    lo = logging.getLogger()
    lo.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    lo.addHandler(ch)
    db =  mydb(lo)
    run_by_gene('GRN','NM_002087.2',db,lo)
    run_by_gene('MAPT','NM_005910.5',db,lo)
    run_by_gene('SOD1','NM_000454.4',db,lo)

if __name__ == '__main__':
    test()
