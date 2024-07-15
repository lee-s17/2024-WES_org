#!/bin/bash

input="/Volumes/mum2023/z6/hgdp_raw"
output="/Volumes/mum2023/z6/Documents/workingdir/WES/e112"

##call gatk
gatk="/Volumes/mum2023/z6/Documents/workingdir/softwares/gatk-4.6.0.0/./gatk"
plink2="/Volumes/mum2023/z6/Documents/workingdir/softwares/./plink2"

#######Input: select variants in ORG
for file in  $(ls $input/*.vcf.gz); do
        bn="$(basename "$file" | cut -c24-)"
        $gatk SelectVariants -O $output/hgdporg_$bn  -V $file --L /Volumes/mum2023/z6/Documents/workingdir/WES/bedfile/e112.list
done

for k in {1..22}; do
  echo $input/hgdporg_chr$k.vcf.gz >> $output/count.txt
  zgrep -v "#" $input/hgdporg_chr$k.vcf.gz | wc -l >> $output/count.txt
done
concatenate org variants
for k in {1..22}; do 
 echo $output/hgdporg_chr$k.vcf.gz >> $output/list
  done

######choose ORG variants
oute110="/Volumes/mum2023/z6/Documents/workingdir/WES/e112"
ORGlist="/Volumes/mum2023/z6/Documents/workingdir/WES/bedfile/e112.list"


passsnp="/Volumes/mum2023/z6/Documents/workingdir/WES/postVC/postVQSR/annpass-OA_snp.vcf.gz"
passindel="/Volumes/mum2023/z6/Documents/workingdir/WES/postVC/postVQSR/annpass-OA_indel.vcf.gz"
passbn="$(basename $passsnp | cut -d "_" -f 1)"

$gatk SelectVariants -V $passsnp  -O $oute110/ORG-"$passbn"_snp.vcf.gz --L $ORGlist
$gatk SelectVariants -V $passindel  -O $oute110/ORG-"$passbn"_indel.vcf.gz --L $ORGlist


######annotate vcf files
pop="/Volumes/mum2023/z6/Documents/workingdir/WES/postVC/postVQSR/bcf-pop.txt"
region="/Volumes/mum2023/z6/Documents/workingdir/WES/postVC/postVQSR/bcf-region.txt"
for file in $(ls ORG-annpass-OA_indel.vcf.gz ORG-annpass-OA_snp.vcf.gz)
do
bn="$(basename $file | cut -d "." -f 1)"
bcftools view $file -S $pop -Ou | bcftools +fill-tags -o bcf-"$bn"  -- -S $region -t AN,AC,AF
echo "CHROM POS ID REF ALT AN AF AC AN_JH AF_JH AC_JH AN_TM AF_TM AC_TM" | sed 's/ /\t/g' >> bcfann-"$bn".txt
bcftools query --format  '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AN\t%INFO/AF\t%INFO/AC\t%INFO/AN_JH\t%INFO/AF_JH\t%INFO/AC_JH\t%INFO/AN_TM\t%INFO/AF_TM\t%INFO/AC_TM\n' bcf-"$bn" >> bcfann-"$bn".txt
done
for file in $(ls hgdporg-auto.vcf.gz hgdporg_chrX.vcf.gz); do
        bcftools view -S /Volumes/mum2023/z6/Documents/workingdir/WES/postVC/e110post/hgdp/hgdpsample-king $file -o bcf-"$(basename $file | cut -d "." -f 1)" --force-samples
        bcftools +fill-tags bcf-"$(basename $file | cut -d "." -f 1)" -o bcfann-"$(basename $file | cut -d "." -f 1)"  -- -S /Volumes/mum2023/z6/Documents/workingdir/hgdpdec20/2022hgdp/vcf_org/kingfilter-afrnon-sample-region -t AN,AC,AF  
        echo "CHROM POS ID REF ALT AA_chimp AN AF AC AN_AFR AF_AFR AC_AFR AN_NONAFR AF_NONAFR AC_NONAFR" | sed 's/ /\t/g' >> bcfannfreq-"$(basename $file | cut -d "." -f 1)" 
        bcftools query --format  \
        '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AA_chimp\t%INFO/AN\t%INFO/AF\t%INFO/AC\t%INFO/AN_AFR\t%INFO/AF_AFR\t%INFO/AC_AFR\t%INFO/AN_NONAFR\t%INFO/AF_NONAFR\t%INFO/AC_NONAFR\n' \
        bcfann-"$(basename $file | cut -d "." -f 1)" >> bcfannfreq-"$(basename $file | cut -d "." -f 1)";
done
plink geno/hwe, fst for OA ORG
$plink2 --set-all-var-ids @:# --hwe 10e-6 --geno 0.05 --remove ../postVC/postVQSR/king1table.king.cutoff.out.id \
    --out ORG-annpass-genohweking_snp --make-pgen --vcf ORG-annpass-OA_snp.vcf.gz \
    --update-sex ../postVC/e110post/sexf.txt
$plink2 --set-all-var-ids @:# --hwe 10e-6 --geno 0.05 --remove ../postVC/postVQSR/king1table.king.cutoff.out.id \
    --out ORG-annpass-genohweking_indel --make-bed --vcf ORG-annpass-OA_indel.vcf.gz \
    --update-sex ../postVC/e110post/sexf.txt
$plink2  --pfile ORG-annpass-genohweking_snp --fst PHENO method=hudson report-variants \
 --out ORG-genohweking_snp --update-sex ../postVC/e110post/sexf.txt --remove ../postVC/postVQSR/king1table.king.cutoff.out.id
 $plink2  --bfile ORG-annpass-genohweking_indel --fst PHENO method=hudson report-variants \
 --out ORG-genohweking_indel --update-sex ../postVC/e110post/sexf.txt --remove ../postVC/postVQSR/king1table.king.cutoff.out.id

king="../postVC/e110post/hgdp/hgdpsample-king"
$plink2 --set-all-var-ids @:# --hwe 10e-6  'midp'  --keep $king \
        --out hgdporg-auto-genohwe --make-pgen \
        --vcf hgdporg-auto.vcf.gz
$plink2 --set-all-var-ids @:# --hwe 10e-6  --keep $king \
        --out hgdporg-chrx-genohwe --make-pgen --update-sex sexhgdp.txt \
        --vcf hgdporg_chrx.vcf.gz


