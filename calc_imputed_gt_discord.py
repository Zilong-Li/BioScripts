#!/usr/bin/env python3

import sys
import argparse
import gzip
import numpy as np
import matplotlib.pyplot as plt

from cyvcf2 import VCF


def _get_genotype_list(gts):

    # half call is set to missing
    res = [gt[0] + gt[1] if gt[0] != -1 and gt[1] != -1 else -1 for gt in gts]
    return res


def _count_discord_each_type(gts1, gts2, dct):

    assert(len(gts1) == len(gts2))
    out = []
    for i in range(len(gts1)):
        if gts1[i] == 0 and gts2[i] == 1:
            dct["0->1"] += 1
            out.append("01")
        elif gts1[i] == 0 and gts2[i] == 2:
            dct["0->2"] += 1
            out.append("02")
        elif gts1[i] == 1 and gts2[i] == 0:
            dct["1->0"] += 1
            out.append("10")
        elif gts1[i] == 1 and gts2[i] == 2:
            dct["1->2"] += 1
            out.append("12")
        elif gts1[i] == 2 and gts2[i] == 0:
            dct["2->0"] += 1
            out.append("20")
        elif gts1[i] == 2 and gts2[i] == 1:
            dct["2->1"] += 1
            out.append("21")
        elif gts1[i] == gts2[i]:
            dct["equal"] += 1
            out.append("11")
        elif gts1[i] == -1 or gts2[i] == -1:
            dct["missing"] += 1
            out.append("00")
        else:
            sys.stderr.write(
                f"something wrong when parsing the genotypes!\n{gts1[i]}\t{gts2[i]}")
            sys.exit(1)

    return out


def _count_discord_each_type2(dps, gts1, gts2, dct, dct2):

    assert(len(gts1) == len(gts2))
    out = []
    for i in range(len(gts1)):
        if gts1[i] == 0 and gts2[i] == 1:
            dct["0->1"] += 1
            dct2["0->1"].append(dps[i][0])
            out.append("01")
        elif gts1[i] == 0 and gts2[i] == 2:
            dct["0->2"] += 1
            dct2["0->2"].append(dps[i][0])
            out.append("02")
        elif gts1[i] == 1 and gts2[i] == 0:
            dct["1->0"] += 1
            dct2["1->0"].append(dps[i][0])
            out.append("10")
        elif gts1[i] == 1 and gts2[i] == 2:
            dct["1->2"] += 1
            dct2["1->2"].append(dps[i][0])
            out.append("12")
        elif gts1[i] == 2 and gts2[i] == 0:
            dct["2->0"] += 1
            dct2["2->0"].append(dps[i][0])
            out.append("20")
        elif gts1[i] == 2 and gts2[i] == 1:
            dct["2->1"] += 1
            dct2["2->1"].append(dps[i][0])
            out.append("21")
        elif gts1[i] == gts2[i]:
            dct["equal"] += 1
            out.append("11")
        elif gts1[i] == -1 or gts2[i] == -1:
            dct["missing"] += 1
            dct2["missing"].append(dps[i][0])
            out.append("00")
        else:
            sys.stderr.write(
                f"something wrong when parsing the genotypes!\n{gts1[i]}\t{gts2[i]}")
            sys.exit(1)

    return out


def calc_beagle_glgt_discord_rate(theInVcf, theOutVcf, theOutPref, chrom=None):

    vcfin = VCF(theInVcf)
    vcfout = VCF(theOutVcf)
    out = gzip.open(f"{theOutPref}.log.gz", 'wt')
    # can not just run vcfin(chrom)
    if chrom is not None:
        vi = next(vcfin(chrom))
        vo = next(vcfout(chrom))
    else:
        vi = next(vcfin)
        vo = next(vcfout)

    dct = {"0->1": 0, "0->2": 0, "1->0": 0, "1->2": 0,
           "2->0": 0, "2->1": 0, "equal": 0, "missing": 0}
    dct2 = {"0->1": [], "0->2": [], "1->0": [],
            "1->2": [], "2->0": [], "2->1": [], "missing": []}
    try:
        while True:
            if vi.start == vo.start:
                stats = _count_discord_each_type2(vi.format("DP"), _get_genotype_list(
                    vi.genotypes), _get_genotype_list(vo.genotypes), dct, dct2)
                out.write(vi.CHROM + "\t" + str(vi.start+1) +
                          "\t" + "\t".join(stats)+"\n")
                vi = next(vcfin)
                vo = next(vcfout)
            elif vi.start > vo.start:
                vo = next(vcfout)
            elif vi.start < vo.start:
                vi = next(vcfin)
            else:
                pass
    except StopIteration:
        sys.stdout.write("reach the end of vcf!\n")

    res = []
    labels = []
    sc = 0
    for k, v in dct.items():
        sc += v
        if k != "equal":
            res.append(v)
            labels.append(k)
        out.write(f"{k}\t{v}\n")

    # res = [round(i / sc, 6) for i in res]
    res.append(sum(res))
    labels.append("TotalDisc")

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    x = np.arange(len(labels)) + 1
    rects = ax.bar(x, res, label=f"#all samples genotypes({sc})")
    # ax.ticklabel_format(style="plain")
    ax.set_xlabel('Genotype Types')
    ax.set_ylabel('Discordance Counts')
    ax.set_title('Discordance Counts between original and imputed vcf')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()
    ax.bar_label(rects)

    ax2.set_title('Violine plot of Format/DP')
    ax2.set_xlabel('Genotype Types')
    ax2.set_ylabel('Format/DP')
    ax2.violinplot(dct2.values())
    # avoid userwarning  FixedLocator
    ax2.set_xticks(np.arange(len(dct2.keys())) + 1)
    ax2.set_xticklabels(dct2.keys())
    fig.tight_layout()
    plt.savefig(f"{theOutPref}.png", dpi=300)


def main():

    parser = argparse.ArgumentParser(
        description="calculate the GT discordance rate between two vcfs.")
    parser.add_argument("vcf1", metavar="VCF1",
                        help="vcf1 with GT or GP tag")
    parser.add_argument("vcf2", metavar="VCF2",
                        help="vcf2 with GT or GP tag")
    parser.add_argument("out", metavar="OUT",
                        help="output prefix")
    parser.add_argument("-chr", metavar="STRING",
                        help="chromosome to use")
    args = parser.parse_args()

    calc_beagle_glgt_discord_rate(
        args.vcf1, args.vcf2, args.out, args.chr)


if __name__ == '__main__':
    main()
