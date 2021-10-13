#!/usr/bin/env python3
import sys
import gzip
import argparse

def _parse_sites(theBeagleSitesFile):

    d = {}
    with open(theBeagleSitesFile) as f:
        for row in f:
            tmp = row.rstrip().split()
            # pos, ref, alt
            d[tmp[0]] = [tmp[1], tmp[2], tmp[3]]
    return d


def _output_sample_file(header, theOutPrefix):
    '''
    header of phased.gz: I id sample1 sample2 ...
    '''

    samples = header.rstrip().split()[2:]
    out = open(f"{theOutPrefix}.samples", "w")
    header_of_sample = "ID_1 ID_2 missing\n0 0 0\n"
    out.write(header_of_sample)
    for s in samples[::2]:
        out.write(f"{s} {s} 0\n")


def run(theChr, theBeaglePhasedFile, theBeagleSitesFile, theOutPrefix):
    '''
    output: hap.gz/samples format of shapeit2 but coustomized for bcftools convert
    see:
        http://www.htslib.org/doc/1.1/bcftools.html#convert
        https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#hapsample
    '''

    sites = _parse_sites(theBeagleSitesFile)

    out = gzip.open(f"{theOutPrefix}.hap.gz", "wt")
    # get sample from header of phased.gz
    line = 0
    with gzip.open(theBeaglePhasedFile, "rt") as f:
        for row in f:
            line += 1
            if line == 1:
                _output_sample_file(row, theOutPrefix)
            else:
                tmp = row.rstrip().split()
                snp = tmp[1]
                if sites.get(snp):
                    d = {sites[snp][1]: '0', sites[snp][2]: '1'}
                    gts = [d[g] for g in tmp[2:]]
                    out.write(f"{theChr}:{sites[snp][0]}_{sites[snp][1]}_{sites[snp][2]} {snp} {sites[snp][0]} {sites[snp][1]} {sites[snp][2]} {' '.join(gts)}\n")
                else:
                    sys.stderr.write("the sites file doesn't match the phased file!")
                    sys.exit(1)



def main():

    parser = argparse.ArgumentParser(
        description="calculate the GT discordance rate between two vcfs.")
    parser.add_argument("chr", metavar="STRING",
                        help="chromosome to use")
    parser.add_argument("phased", metavar="phased.gz",
                        help="phased.gz file by beagle v3")
    parser.add_argument("sites", metavar="sites",
                        help="sites file by beagle v3")
    parser.add_argument("out", metavar="OUT",
                        help="output prefix")
    args = parser.parse_args()

    run(args.chr, args.phased, args.sites, args.out)


if __name__ == '__main__':
    main()
