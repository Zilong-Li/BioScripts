#!/bin/bash

set -e

usage() { echo "Usage: $0 [-b <bamlist>] [-d <depth>] [-n <cores>] [-o <outdir>]" 1>&2; exit 1; }

dry_run=false

while getopts ":b:d:n:o:D" o; do
    case "${o}" in
        b)
            bamlst=${OPTARG}
            ;;
        d)
            deps=${OPTARG}
            ;;
        n)
            cores=${OPTARG}
            ;;
        o)
            OUT=${OPTARG}
            ;;
        D)
            dry_run=true
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${bamlst}" ] || [ -z "${deps}" ] || [ -z "${cores}" ] || [ -z "${OUT}" ]; then
    usage
fi

if "${dry_run}"; then
    cmd="--dry-run"
else
    cmd=""
fi


# OUT=results/downsample
# parallel --dry-run

# deps=(0.5 1.0) # target depth not fraction

for dep in ${deps[@]};do
  outdir=$OUT/${dep}x
  mkdir -p $outdir
  for bam in `cat $bamlst`;do 
    out=$outdir/`basename $bam`.bam
    echo $out >>$outdir/bams.list
    echo "frac=\`samtools view -h $bam chr21 | samtools depth -@ 2 -a - | awk '{s+=\$3;} END{print $dep*NR/s}'\`;samtools view -h -s \$frac -o $out $bam; samtools index $out"
  done | parallel $cmd -j $cores -k "sh -c {}"
done

wait
echo "downsampling finished"

