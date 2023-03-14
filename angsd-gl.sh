#!/bin/bash

set -e

usage() { echo "Usage: $0 [-f <fam>] [-b <bamlist>] [-s <sites>] [-o <outdir>]" 1>&2; exit 1; }

while getopts ":f:b:s:o:" o; do
    case "${o}" in
        b)
            bamlst=${OPTARG}
            ;;
        f)
            fam=${OPTARG}
            ;;
        s)
            sites=${OPTARG}
            ;;
        o)
            outdir=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${fam}" ] || [ -z "${bamlst}" ] || [ -z "${sites}" ] || [ -z "${outdir}" ]; then
    usage
fi

mkdir -p $outdir
newbam=$outdir/order.bamlst
R -s -e "
fam=read.table('$fam',h=F)
ids=fam[,1]
stopifnot(length(unique(ids))==dim(fam)[1])
bam=read.table('$bamlst',h=F)
bamid=sapply(strsplit(sapply(strsplit(bam[,1], '/'), '[[', 4), '.', fix=T), '[[', 1)
rownames(bam) <- bamid
cat(bam[ids,], sep='\n')
" > $newbam

angsd="angsd"

# sites=/emc/zilong/phaseless/tmp/$prefix/all.sites
# awk '{print $1"\t"$4"\t"$5"\t"$6}' ${bfile}.bim >${sites}
# sleep 3
# ${angsd} sites index ${sites}


chrs=($(cut -f1 ${sites}|uniq))
echo $chrs

for chr in "${chrs[@]}";do
{
  outgl=$outdir/$chr/genolike
  mkdir -p `dirname $outgl`
  ${angsd} -bam ${newbam} -sites ${sites} -GL 1 -out $outgl -doGlf 2 -doMajorMinor 3 -r $chr
} &
done
wait
echo "angsd -GL done"

for chr in "${chrs[@]}";do
  zcat $outdir/$chr/genolike.beagle.gz
done | awk 'NR>1 && /^marker/ {next} 1' | gzip -c >$outdir/genolike.beagle.gz
echo "concat all gz files done"
