#!/bin/bash

# trace date
PS4='+ $(date "+%x %T %Z")\011 '
# terminate the script immediately if there is a error
set -e

species=giraff
version=v3
dir=/home/zilong/project/african1k/$species
outdir=$dir/output/imputation/$version
mkdir -p $outdir
threads=40
# vcf=/home/krishang/projects/DNA/africa1kg/sampleQC/giraf/RothschildsGiraffe_GCam_withbaq_variable_sites_diallelic_king01_maf005.vcf.gz
# vcf=/home/zilong/project/african1k/buffalo/output/calls/all.buffalo.maf0.05.biallelic.sites.vcf.gz
# convert PL to GL first
# bcftools +tag2tag $vcfori -Oz -o $vcf --threads $threads -- -r --pl-to-gl && bcftools index -f $vcf && echo done

beagle41=/home/zilong/local/bin/beagle.27Jan18.7e1.jar
beagle3=/home/genis/software/beagle.jar

showhelp() {
# `cat << EOF` This means that cat should stop reading when EOF is detected
cat <<EOF
Usage: ./call-genotype-lcs.sh [command]
Pipeline of genotype calling and imputation for median depth sequencing data

-h,          Display help

beagle3,     Run beagle 3

postqc3,     Run Post QCs after imputation by Beagle 3

EOF
}

runBeagle4()
{
    vcf=/home/krishang/projects/DNA/africa1kg/imputation_datasets/giraf2/RothschildsGiraffe_GCam_variable_sites_nomultiallelics_noindels_nodups_maf005.vcf.gz
    # show verbose log
    set -x
    echo "start running beagle4 by chroms"

    for chrom in {14..1};do
        out=$outdir/$chrom.gl.imputed

        java -Xss5m -Xmx60g -jar $beagle41 gl=$vcf out=$out chrom=$chrom nthreads=$threads
        echo "convert $out.vcf.gz to tabixed vcf"
        gzip -dc $out.vcf.gz | bgzip -c >$out.t.vcf.gz && mv -f $out.t.vcf.gz $out.vcf.gz && tabix -f $out.vcf.gz

        echo "running: calc_imputed_glgt_discord.py -invcf $vcf -outvcf $out.vcf.gz -out $out.sum.disc.count -chrom $chrom"

        calc_imputed_gt_discord.py $vcf $out.vcf.gz $out.sum.disc.count $chrom
        echo -e "ID\tDR2\tINFO\tAF" >$out.info
        bcftools +impute-info $out.vcf.gz | bcftools query -f '%CHROM:%POS:%REF:%ALT\t%INFO/DR2\t%INFO/INFO\t%INFO/AF\n' >>$out.info

        # plot
        plot-r2-vs-maf.R $out.info $out.info

    done
}

runBeagle3() {
    # show verbose log
    set -x
    VCF=/home/krishang/projects/DNA/africa1kg/imputation_datasets/wildebeest/Wildebeest_wildebeest_variable_sites_nomultiallelics_noindels_maf005.bcf.gz # vcf path
    OUT=`pwd` # outputdir
    mkdir -p $OUT/vcf $OUT/imputed
    echo "start running beagle3 by chr"
    chrs=$(bcftools index -s $VCF|cut -f1)   # may change chr name patterns to your own
    for chrom in $chrs;do
    {
        # convert PL tag to GL in beagle format 
        indir=$OUT/vcf
        bname=$indir/`basename $VCF`.$chrom
        bcftools view -r $chrom -O u $VCF | bcftools +tag2tag -Ov -o ${bname}.vcf -- -r --pl-to-gl && vcftools --vcf ${bname}.vcf --out ${bname} --BEAGLE-GL --chr $chrom;

        # run imputation
        outdir=$OUT/imputed
        out=$outdir/$chrom
        input=${bname}.BEAGLE.GL
        java -Xss5m -Xmx40g -jar $beagle3 seed=2 like=$input out=$out omitprefix=true
        zcat $out.BEAGLE.GL.gprobs.gz | java -jar /home/zilong/local/bin/gprobs2beagle.jar 0.9 -1 | gzip -c >$out.bgl.gz && \
        zcat $out.BEAGLE.GL.gprobs.gz | awk 'NR>1{split($1,a,":");print $1,a[2],$2,$3}' >$out.bgl.sites && \
        java -jar /home/zilong/local/bin/beagle2vcf.jar $chrom $out.bgl.sites $out.bgl.gz -1 |bgzip -c >$out.vcf.gz && \
        bcftools index -f $out.vcf.gz && calc_imputed_gt_discord.py $vcf $out.vcf.gz $out.sum.disc.count $chrom && \

        # convert phased.gz to hap then to vcf
        beagle_phased_to_hap_sample.py $chrom $out.BEAGLE.GL.phased.gz $out.bgl.sites $out && bcftools convert --hapsample2vcf $out |bcftools annotate -I +'%CHROM:%POS'-Oz -o $out.phased.vcf.gz && bcftools index -f $out.phased.vcf.gz
        echo "$chrom done"
    } &
    done
    wait

    echo "all jobs done by chroms"

    echo "start concating files"
    out=$outdir/all
    bcftools concat --threads 10 -Ob -o $out.bcf `for i in $chrs; do echo $outdir/$i.vcf.gz;done` && bcftools index -f $out.bcf
    bcftools concat --threads 10 `for i in $chrs; do echo $outdir/$i.phased.vcf.gz;done` | bcftools annotate -I +'%CHROM:%POS' -Ob -o $out.phased.bcf --threads 10  && bcftools index -f $out.phased.bcf
    for i in $chrs;do
        cat $outdir/$i.BEAGLE.GL.r2 
    done >$out.r2

    echo "beagle3 imputation done"
}

runBeagle3PostQC() {
    # show verbose log
    set -x

    vcfin=$outdir/all.bcf
    vcfphased=$outdir/all.phased.bcf
    out=$outdir/all
    stats=$out.imputed.sites.stats
    # get maf of sites
    bcftools +fill-tags $vcfin -- -t MAF |bcftools query -f "%ID\t%MAF\n" >$out.maf
    # merge r2 and maf together
    echo -e "ID\tR2\tMAF" >$stats
    paste $out.r2 <(awk '{print $2}' $out.maf) >>$stats
    # apply filters MAF>0.05 && R2 >0.95
    # remove NaN sites in R2 file
    # awk 'NR>1 && $3>0.05 && $2>0.95 && $2!="NaN"{print $1}' $stats >$out.maf0.05.r2.0.95.sites
    awk 'NR>1 && $3>0.05 && $2>0.99 && $2!="NaN"{print $1}' $stats >$out.maf0.05.r2.0.99.sites
    bcftools view -i ID=@$out.maf0.05.r2.0.99.sites -Ob -o $out.maf0.05.r2.0.99.bcf --threads 10 $vcfin && bcftools index -f $out.maf0.05.r2.0.99.bcf
    # # for phased vcf
    bcftools view -i ID=@$out.maf0.05.r2.0.99.sites -Ob -o $out.phased.maf0.05.r2.0.99.bcf --threads 10 $vcfphased && bcftools index -f $out.phased.maf0.05.r2.0.99.bcf
}

case ${@:$OPTIND:1} in
    beagle3)
        runBeagle3
        ;;
    postqc3)
        runBeagle3PostQC
        ;;
    "")
        showhelp
        exit 0
        ;;
esac

# for people likes -h
while getopts 'h' OPT; do
    case $OPT in
        h)  showhelp;;
    esac
done
