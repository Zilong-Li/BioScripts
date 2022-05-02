#!/usr/bin/env python3

import sys
import click
import warnings

from cyvcf2 import VCF
from scipy.stats.stats import pearsonr


def _get_samples(f_sam):
    '''
    verify the input sample file
    return list of id_imp and id_true, respectively
    '''
    id_imps, id_trues = [], []
    id1, id2 = 'id_imp', 'id_true'
    with open(f_sam) as I:
        head = I.readline().rstrip().split()
        d = dict((k, i) for i, k in enumerate(head))
        if id1 in d and id2 in d:
            for row in I:
                t = row.rstrip().split()
                id_imps.append(t[d[id1]])
                id_trues.append(t[d[id2]])
        else:
            sys.exit(
                'Error: first row of the sample file must include id_imp and id_true')
    return id_imps, id_trues


def _dosage(gps):
    '''
    gps: v.genotypes
    return a list of dosage
    '''
    dos = [1 if gp[0] == -1 or gp[1] == -1 else gp[0] + gp[1] for gp in gps]
    return dos


def _dosage2(gps, alt_i):
    '''
    gps: Nx3 mat
    alt_i: index of alt
    return a list of dosage
    '''
    dos = [gp[1] + gp[alt_i] * 2 for gp in gps]
    return dos


def _get_index(ids1, ids2):

    d = dict((k, i) for i, k in enumerate(ids2))
    idx = []
    for k in ids1:
        if k in d:
            idx.append(d[k])
        else:
            sys.exit(k + ' does not exist in ' + str(ids2))

    return idx


def _by_index(lst, idx):

    return [lst[i] for i in idx]


def _trueset(f_vcf, region, id_trues, tds):
    '''
    return a dict
    key: chr_pos
    value: a list of dosage
    '''
    vcf = VCF(f_vcf, samples=id_trues)
    samples = vcf.samples
    idx = _get_index(id_trues, samples)
    dd = {}
    for v in vcf(region):
        if not v.is_snp or v.FILTER:
            continue
        k = v.CHROM + '_' + str(v.end)
        alt = v.ALT[0]
        if tds:
            dos = v.format('DS').T[0]
        else:
            dos = _dosage(v.genotypes)
        dd[k] = {'dos': _by_index(dos, idx), 'alt': alt}
    vcf.close()

    return dd


def corr_vcf(f_vcf, region, id_imps, d_true, OUT, info_t):
    vcf = VCF(f_vcf, samples=id_imps)
    samples = vcf.samples
    idx = _get_index(id_imps, samples)
    if info_t:
        OUT.write("ID\tR2\tP\tEAF\tINFO_SCORE\tHWE\n")
    else:
        OUT.write("ID\tR2\tP\n")

    for v in vcf(region):
        if not v.is_snp:
            continue
        k = v.CHROM + '_' + str(v.end)
        if k in d_true:
            dos1 = d_true[k]['dos']
            dos2 = _by_index(v.format('DS').T[0], idx)
            assert len(dos1) == len(dos2)
            r, p = pearsonr(dos1, dos2)
            r2 = r * r
            if info_t:
                info = str(v.INFO['EAF']) + '\t' + \
                    str(v.INFO['INFO_SCORE']) + '\t' + str(v.INFO['HWE'])
                OUT.write(k + '\t' + str(r2) + '\t' +
                          str(p) + '\t' + info + '\n')
            else:
                OUT.write(k + '\t' + str(r2) + '\t' + str(p) + '\n')
        else:
            sys.stderr.write(k + ' dose not exist in the trueset\n')

    vcf.close()


def corr_bgen(f_bgen, id_imps, d_true, OUT):
    from bgen_reader import read_bgen

    OUT.write("ID\tR2\tP\n")
    bgen = read_bgen(f_bgen, verbose=False)
    samples = bgen["samples"].tolist()
    idx = _get_index(id_imps, samples)
    variants = bgen["variants"]
    i = -1
    for row in variants.itertuples():
        k = getattr(row, 'chrom') + '_' + str(getattr(row, 'pos'))
        allele_ids = getattr(row, 'allele_ids')
        alleles = allele_ids.split(',')
        i += 1
        if k in d_true:
            alt_i = 2 if alleles[1] == d_true[k]['alt'] else 0
            g = bgen["genotype"][i].compute()
            dos2 = _dosage2(g["probs"][idx], alt_i)
            dos1 = d_true[k]['dos']
            assert len(dos1) == len(dos2)
            r, p = pearsonr(dos1, dos2)
            r2 = r * r
            OUT.write(k + '\t' + str(r2) + '\t' + str(p) + '\n')
        else:
            sys.stderr.write(k + ' dose not exist in the trueset\n')


@click.command()
@click.argument('input_type', type=click.Choice(['vcf', 'bgen']))
@click.option('-i', '--input', help='specify the imputed dataset, either vcf or bgen')
@click.option('-o', '--out', help='output file')
@click.option('-t', '--true', help='specify a tabixed vcf as the trueset')
@click.option('-s', '--sam', help='specify a sample file whose header include id_imp and id_true')
@click.option('-r', '--region', help='specify a samtools-like region<chr:start-end>')
@click.option('--tds', is_flag=True, help='use DS field rather than GT from the trueset')
@click.option('--info', is_flag=True, help='output INFO field in the imputed vcf')
def main(input_type, input, out, sam, true, region, tds, info):
    id_imps, id_trues = _get_samples(sam)
    d_true = _trueset(true, region, id_trues, tds)
    OUT = open(out, 'w')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if input_type == 'bgen':
            corr_bgen(input, id_imps, d_true, OUT)
        elif input_type == 'vcf':
            corr_vcf(input, region, id_imps, d_true, OUT, info)
        else:
            sys.stderr.write('Error: input file must be vcf or bgen\n')


if __name__ == '__main__':
    main()
