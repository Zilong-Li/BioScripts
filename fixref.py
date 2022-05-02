#!/usr/bin/env python3

import argparse

from cyvcf2 import VCF, Writer


def _read_ids(theInVcf):
    """
    k = 1:13417:C:CGAGA
    v = C:CGAGA
    """
    d = {}
    vcfin = VCF(theInVcf)
    for v in vcfin:
        k = (
            v.CHROM.replace("chr", "")
            + ":"
            + str(v.start + 1)
            + ":"
            + v.REF
            + ":"
            + v.ALT[0]
        )
        d[k] = v.REF + ":" + v.ALT[0]
    vcfin.close()
    return d


def _flip_genotypes(gts):
    """
    v.genotypes returns a list Indicating the allele and phasing.
    e.g. [0, 1, True] corresponds to 0|1 while [1, 2, False] corresponds to 1/2
    also, -1 means missing
    """
    d = {0: 1, 1: 0, -1: -1}
    for i in range(len(gts)):
        gts[i][0] = d[gts[i][0]]
        gts[i][1] = d[gts[i][1]]

    return gts


def fix(theVcfRef, theInVcf, theOutVcf):

    ids = _read_ids(theVcfRef)
    vcfin = VCF(theInVcf)
    # create a new vcf Writer using the input vcf as a template.
    vcfout = Writer(theOutVcf, vcfin)
    for v in vcfin:
        k1 = (
            v.CHROM.replace("chr", "")
            + ":"
            + str(v.start + 1)
            + ":"
            + v.REF
            + ":"
            + v.ALT[0]
        )
        k2 = (
            v.CHROM.replace("chr", "")
            + ":"
            + str(v.start + 1)
            + ":"
            + v.ALT[0]
            + ":"
            + v.REF
        )
        if k1 in ids or k2 in ids:
            k = v.ALT[0] + ":" + v.REF
            if any(k == ids.get(i) for i in (k1, k2)):
                # ref and alt are switched
                ta = v.ALT[0]
                v.ALT = [v.REF]
                v.REF = ta
                v.genotypes = _flip_genotypes(v.genotypes)
        vcfout.write_record(v)

    vcfout.close()
    vcfin.close()


def main():
    parser = argparse.ArgumentParser(
        description="flip ref and alt and corresponding genotypes using the ref vcf."
    )
    parser.add_argument("ref", metavar="REF_VCF", help="vcf file used as ref.")
    parser.add_argument("vcf", metavar="IN_VCF", help="input vcf to be fixed")
    parser.add_argument("out", metavar="OUT_VCF", help="output vcf")
    args = parser.parse_args()

    fix(args.ref, args.vcf, args.out)


if __name__ == "__main__":
    main()
