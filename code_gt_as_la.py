#!/usr/bin/env python3

import sys
import argparse
from cyvcf2 import VCF, Writer


def _get_rfmix_samples_idx(header):

    d = {it[:-2]: i-1 for i, it in enumerate(header.rstrip().split("\t")[6:])}

    return d


def _parse_msp_line(line, msp_idx):
    """
    return : [spos, epos, [gts]]
    """
    line = line.rstrip().split("\t")
    spos = int(line[1])
    epos = int(line[2])
    line = line[6:]
    gts = [[int(line[i]), int(line[i+1]), True] for i in msp_idx]

    return list([spos, epos, gts])


def _code_ref_as_inuit(gts, gts_la, idx_la):

    out = []
    for i in range(len(gts)):
        if i in idx_la:
            out.append(gts_la[idx_la[i]])
        else:
            out.append([1, 1, True])

    return out


def recode_gts(theMsp, theInVcf, theOutVcf, chrom=None):

    MSP = open(theMsp, 'r')
    MSP.readline()
    d_idx = _get_rfmix_samples_idx(MSP.readline())

    vcfin = VCF(theInVcf)
    samples = vcfin.samples
    # remove samples not in the header of msp.tsv, ie. samples as ref
    keeps = []
    msp_idx = []
    idx_la = {}
    for i in range(len(samples)):
        if samples[i] in d_idx:
            keeps.append(samples[i])
            msp_idx.append(d_idx[samples[i]])
            idx_la[i] = len(keeps)-1
        else:
            pass

    # vcfin.set_samples(keeps)
    if chrom is not None:
        vi = next(vcfin(chrom))
    else:
        vi = next(vcfin)

    vcfout = Writer(theOutVcf, vcfin)
    line = MSP.readline()
    msp = _parse_msp_line(line, msp_idx)
    eof = False
    try:
        while True:
            if vi.start + 1 < msp[0]:
                vi = next(vcfin)
            elif vi.start + 1 >= msp[0] and vi.start + 1 <= msp[1]:
                # vi.genotypes = msp[2]
                vi.genotypes = _code_ref_as_inuit(vi.genotypes, msp[2], idx_la)
                # todo: remove info
                vcfout.write_record(vi)
                vi = next(vcfin)
            elif vi.start + 1 > msp[1]:
                if eof is False:
                    line = MSP.readline()
                    if line != '':
                        msp = _parse_msp_line(line, msp_idx)
                    else:
                        sys.stdout.write("reach the end of msp.csv file!\n")
                        eof = True
                else:
                    break
    except StopIteration:
        sys.stdout.write("reach the end of vcf!\n")


def main():

    parser = argparse.ArgumentParser(
        description='coding rfmix local ancestry as fake genotype.')
    parser.add_argument("-msp", metavar="FILE",
                        help="msp.csv file by rfmix.")
    parser.add_argument("-vcf", metavar="FILE",
                        help="input vcf to be recoded")
    parser.add_argument("-out", metavar="FILE",
                        help="output vcf")
    parser.add_argument("-chrom", metavar="STRING",
                        help="chromosome to use")
    args = parser.parse_args()
    assert (args.msp is not None) or (args.vcf is not None) or (args.out is not None), \
        "\n please feed the correct params! try -h,--help option!\n"

    recode_gts(args.msp, args.vcf, args.out, args.chrom)


if __name__ == '__main__':
    main()
