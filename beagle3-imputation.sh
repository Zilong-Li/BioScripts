#!/usr/bin/env bash

# trace date
PS4='+ $(date "+%x %T %Z")\011 '
# terminate the script immediately if there is a error
set -e

# programs
bcftools=bcftools
vcftools=vcftools
beagle3=beagle.jar
gprobs2beagle=gprobs2beagle.jar
beagle2vcf=beagle2vcf.jar
beagle2hap=beagle_phased_to_hap_sample.py
gtdiscord=calc_imputed_gt_discord.py

runBeagle3() {
    VCF=$1
    OUT=$2
    MAF=$3
    echo "start running beagle3 by chr"
    # chrs=$(bcftools index -s $VCF | cut -f1 | tail -n+$1 | head -n $2)
    chrs=$(${bcftools} index -s $VCF | cut -f1) # may change chr name patterns to your own
    for chrom in $chrs; do
        {
            # convert PL tag to GL in beagle format
            indir=$OUT/gl
            mkdir -p $indir
            bname=$indir/$(basename $VCF).$chrom
            ${bcftools} view -r $chrom -q $MAF:minor -O u $VCF |
                ${bcftools} +tag2tag -Ov -o ${bname}.vcf -- -r --pl-to-gl &&
                ${vcftools} --vcf ${bname}.vcf --out ${bname} --BEAGLE-GL --chr $chrom &&
                rm -f ${bname}.vcf && echo 'convert to BEAGLE-GL done'
            # run imputation
            outdir=$OUT/imputed
            mkdir -p $outdir
            out=$outdir/$chrom
            input=${bname}.BEAGLE.GL
            java -Xss5m -Xmx20g -Djava.io.tmpdir=$outdir -jar $beagle3 like=$input omitprefix=true out=$out
            zcat $out.BEAGLE.GL.gprobs.gz |
                java -jar ${gprobs2beagle} 0.9 -1 |
                gzip -c >$out.bgl.gz &&
                zcat $out.BEAGLE.GL.gprobs.gz |
                awk 'NR>1{split($1,a,":");print $1,a[2],$2,$3}' >$out.bgl.sites &&
                java -jar ${beagle2vcf} $chrom $out.bgl.sites $out.bgl.gz -1 | bgzip -c >$out.vcf.gz &&
                ${bcftools} index -f $out.vcf.gz &&
                python3 ${gtdiscord} -chr $chrom $VCF $out.vcf.gz $out.sum.disc.count &&
                echo "convert to vcf done"
            ### convert phased.gz to hap then to vcf
            python3 ${beagle2hap} $chrom $out.BEAGLE.GL.phased.gz $out.bgl.sites $out &&
                ${bcftools} convert --hapsample2vcf $out |
                ${bcftools} annotate -I +'%CHROM:%POS' -Oz -o $out.phased.vcf.gz &&
                ${bcftools} index -f $out.phased.vcf.gz
            echo "$chrom done"
        } &
    done
    wait

    echo "all jobs done by chroms"
    echo "start concating files"

    out=$outdir/all
    ${bcftools} concat --threads 10 -Ob -o $out.bcf $(for i in $chrs; do echo $outdir/$i.vcf.gz; done) &&
        ${bcftools} index -f $out.bcf
    ${bcftools} concat --threads 10 $(for i in $chrs; do echo $outdir/$i.phased.vcf.gz; done) |
        ${bcftools} annotate -I +'%CHROM:%POS' -Ob -o $out.phased.bcf --threads 10 &&
        ${bcftools} index -f $out.phased.bcf
    for i in $chrs; do
        cat $outdir/$i.BEAGLE.GL.r2
    done >$out.r2

    echo "beagle3 imputation done"
}

showhelp() {
    # `cat << EOF` This means that cat should stop reading when EOF is detected
    cat <<EOF
Usage: beagle3-imputation.sh [options]
Pipeline of genotype refinement for median depth sequencing data using beagle3

-h,          Display help
-i,          Input VCF/BCF file
-o,          Output folder
-f,          MAF filters before imputation
EOF
    exit 1
}

while getopts ":h:i:f:o:" OPT; do
    case "${OPT}" in
        h)
            showhelp
            ;;
        i)
            VCF=${OPTARG}
            ;;
        o)
            OUT=${OPTARG}
            ;;
        f)
            MAF=${OPTARG}
            ;;
        *)
            showhelp
            ;;
    esac
done
shift $((OPTIND - 1))

if [ -z "${VCF}" ] || [ -z "${OUT}" ] || [ -z "${MAF}" ]; then
    showhelp
fi

runBeagle3 $VCF $OUT $MAF && echo "runBeagle3 done"
