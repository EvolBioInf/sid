#!/usr/bin/env python3
import mysql.connector as mc
import sys

DBUSER = 'anonymous'
DBHOST = 'ensembldb.ensembl.org'
DBDB   = 'mus_musculus_core_90_38'

columns = 'site.chrom,site.pos,gene.stable_id,gene.biotype,exon.exon_id,exon.seq_region_start,exon.seq_region_end,exon.seq_region_strand,exon.phase,exon.end_phase'

query_template = ('SELECT DISTINCT {cols}'
        ' FROM gene as gene'
        '   JOIN exon_transcript AS et ON (gene.canonical_transcript_id = et.transcript_id)'
        '   JOIN exon AS exon USING (exon_id)'
        '   JOIN seq_region AS r ON (exon.seq_region_id = r.seq_region_id)'
        '   JOIN translation AS tr ON (gene.canonical_transcript_id = tr.transcript_id)'
        '   JOIN exon AS first ON (first.exon_id = tr.start_exon_id)'
        '   JOIN exon AS last ON (last.exon_id = tr.end_exon_id)'
        '   JOIN ({positions}) as site'
        ' WHERE exon.seq_region_start <= site.pos'
        '   AND exon.seq_region_end >= site.pos'
        '   AND (exon.seq_region_strand != 1 OR (site.pos >= first.seq_region_start + tr.seq_start - 1 AND site.pos <= last.seq_region_start + tr.seq_end -1))'
        '   AND (exon.seq_region_strand != -1 OR (site.pos >= last.seq_region_end - tr.seq_end + 1 AND site.pos <= first.seq_region_end - tr.seq_start + 1))'
        '   AND r.name = CONVERT(site.chrom USING latin1) AND r.coord_system_id = 3'
        ' ORDER BY site.chrom, site.pos ASC')

def check_position(chrom, pos, cursor):
    position_query = 'select \'{}\' as chrom, {} as pos'.format(chrom, pos)
    query = query_template.format(cols=columns, positions=position_query)
    cursor.execute(query)
    for row in cursor.fetchall():
        print(','.join(map(str, row)), flush=True)

def check_positions(sites):
    cn = mc.connect(user=DBUSER, host=DBHOST, db=DBDB)
    cursor = cn.cursor()

    print('Executing queries...', file=sys.stderr)
    print(columns)
    for chrom, site in sites:
        check_position(chrom, site, cursor)

    cursor.close()
    cn.close()

if __name__ == '__main__':
    if len(sys.argv) > 1:
        print('Usage: supply newline-separated pairs of chrom,position on stdin')
        sys.exit()

    print('Parsing input... ', file=sys.stderr)
    sites = (line.strip().split(',')[:2] for line in sys.stdin if line[0] != '#')
    header = next(sites)
    check_positions(sites)
