#!/usr/bin/env python3
import csv
import gzip
import sys
import collections as col
import itertools as it

SitePhase = col.namedtuple('Phase', ['chrom', 'position', 'gene_id', 'phase', 'strand'])

sid_chrom_col = 'chrom'
sid_position_col = 'pos'
sid_genotype_col = 'gt'

def line_context(iterable, context=2):
    buffer = ['' for _ in range(1 + 2*context)]
    center_index = context
    for i in range(center_index + 1, len(buffer)):
        buffer[i] = next(iterable)
    for line in iterable:
        buffer[:-1] = buffer[1:]
        buffer[-1] = line
        yield tuple(buffer)
    for _ in range(context):
        buffer[:-1] = buffer[1:]
        buffer[-1] = ""
        yield tuple(buffer)

Translation = col.namedtuple('Translation', ['chrom', 'pos', 'gene_id', 'strand', 'phase', 'codons', 'translations'])

def generate_site_codons(sid_lines, phase_records):
    phase_records = sorted(phase_records, reverse=True)
    if len(phase_records) < 1:
        return
    # with open(sid_output_path) as sid_file:
    header = next(sid_lines)
    genotype_col_index = header.split(',').index(sid_genotype_col)
    # position, phase, strand = position_records.pop(0)
    current = phase_records.pop()
    # for context in line_context(sid_file, 2):
    for context in line_context(sid_lines, 2):
        # (hopefully) fast heuristic test before parsing line
        if str(current.position) in context[2]:
            ref_chrom, ref_pos, *_ = context[2].split(',')
            if ref_chrom != current.chrom or int(ref_pos) != current.position:
                # false positive by "in" test
                continue
            if current.strand == 1:
                offset = 2 - current.phase
            elif current.strand == -1:
                offset = current.phase
            multicodon = [line.split(',')[genotype_col_index] for line in context[offset:offset+3]]
            if current.strand == -1:
                multicodon = reverse_complement(multicodon)

            translations = sorted(set(map(translate, all_combinations(multicodon))))
            yield Translation(current.chrom, current.position, current.gene_id, current.strand, current.phase, multicodon, translations)
            if len(phase_records) == 0:
                return
            old = current
            current = phase_records.pop()
            while current.chrom == old.chrom and current.position == old.position:
                yield Translation(current.chrom, current.position, current.gene_id, current.strand, current.phase, multicodon, translations)
                current = phase_records.pop()

def reverse_complement(codon):
    def complement(s):
        return s.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()
    return list(map(complement, codon))[::-1]

def all_combinations(pairs):
    return sorted(set(it.product(*pairs)))

def translate(codon):
    codon = ''.join(codon)
    if len(codon) != 3 or not set(codon).issubset(set('ACGT')):
        raise Exception('Invalid base in codon {}'.format(codon))

    if codon == 'TTT' or codon == 'TTC':
        return 'F'
    if codon.startswith('CT') or codon == 'TTA' or codon == 'TTG':
        return 'L'
    if codon.startswith('AT'):
        if codon == 'ATG':
            return 'M'
        else:
            return 'I'
    if codon.startswith('GT'):
        return 'V'
    if codon.startswith('TC'):
        return 'S'
    if codon.startswith('CC'):
        return 'P'
    if codon.startswith('AC'):
        return 'T'
    if codon.startswith('GC'):
        return 'A'
    if codon.startswith('TA'):
        if codon.endswith('T') or codon.endswith('C'):
            return 'Y'
        else:
            return 'stop'
    if codon.startswith('CA'):
        if codon.endswith('T') or codon.endswith('C'):
            return 'H'
        else:
            return 'Q'
    if codon.startswith('AA'):
        if codon.endswith('T') or codon.endswith('C'):
            return 'N'
        else:
            return 'K'
    if codon.startswith('GA'):
        if codon.endswith('T') or codon.endswith('C'):
            return 'D'
        else:
            return 'E'
    if codon.startswith('TG'):
        if codon.endswith('T') or codon.endswith('C'):
            return 'C'
        elif codon.endswith('A'):
            return 'stop'
        else:
            return 'W'
    if codon.startswith('CG'):
        return 'R'
    if codon.startswith('AG'):
        if codon.endswith('T') or codon.endswith('C'):
            return 'S'
        else:
            return 'R'
    if codon.startswith('GG'):
        return 'G'
    raise Exception('This code should be unreachable')

def parse_ensembl_data(ensembl_output_path):
    # records = []
    with open(ensembl_output_path) as ensembl_file:
        ensembl_data = csv.DictReader(ensembl_file)
        for row in ensembl_data:
            chrom = row['site.chrom']
            pos = int(row['site.pos'])
            exon_phase = int(row['exon.phase'])
            exon_end_phase = int(row['exon.end_phase'])
            exon_start = int(row['exon.seq_region_start'])
            exon_end = int(row['exon.seq_region_end'])
            strand = int(row['exon.seq_region_strand'])
            gene_id = row['gene.stable_id']

            if strand == 1:
                if exon_phase != -1:
                    site_phase = (pos - exon_start + exon_phase) % 3
                elif exon_end_phase != -1:
                    # tested
                    site_phase = (exon_end - pos + exon_end_phase + 1) % 3
                else:
                    site_phase = (pos - exon_start) % 3
            elif strand == -1:
                if exon_phase != -1:
                    site_phase = (exon_end - pos + exon_phase) % 3
                elif exon_end_phase != -1:
                    # tested
                    site_phase = (pos - exon_start + exon_end_phase + 1) % 3
                else:
                    site_phase = (exon_end - pos) % 3
            else:
                continue
            yield SitePhase(chrom, pos, gene_id, site_phase, strand)
            # records.append((pos, site_phase, strand))
    # return records

def output_record(record):
    codons = ':'.join(record.codons)
    translations = ':'.join(record.translations)
    if len(record.translations) > 1:
        print(record.chrom, record.pos, record.gene_id, "nonsyn", record.strand, record.phase, codons, translations, sep=',')
    else:
        print(record.chrom, record.pos, record.gene_id, "syn", record.strand, record.phase, codons, translations,  sep=',')


if __name__ == '__main__' and len(sys.argv) > 1:
    print('parsing positions and phase data...', file=sys.stderr)
    ensembl_data = parse_ensembl_data(sys.argv[1])

    print('generating site triplets...', file=sys.stderr)
    sid_path = sys.argv[2]
    codon_data = []
    if sid_path.endswith('.gz'):
        with gzip.open(sid_path) as f:
            for r in generate_site_codons(map(bytes.decode, f), ensembl_data):
                output_record(r)
    else:
        with open(sid_path) as f:
            for r in generate_site_codons(f, ensembl_data):
                output_record(r)

