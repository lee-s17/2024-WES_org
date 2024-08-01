#!/bin/bash

maindir="/Volumes/mum2023/z6/Documents/workingdir"
output="$maindir/WES/e112"

##call gatk
gatk="$maindir/softwares/gatk-4.6.0.0/./gatk"
plink2="$maindir/softwares/./plink2"
plink="$maindir/softwares/./plink"

# #######Input: select variants in ORG
hgdpin="/Volumes/mum2023/z6/hgdp_raw"
for file in  $(ls $hgdpin/*.vcf.gz); do
        bn="$(basename "$file" | cut -c24-)"
        $gatk SelectVariants -O $output/hgdporg_$bn  -V $file --L $maindir/WES/bedfile/e112.list
done

for k in {1..22}; do
  echo $hgdpin/hgdporg_chr$k.vcf.gz >> $output/count.txt
  zgrep -v "#" $hgdpin/hgdporg_chr$k.vcf.gz | wc -l >> $output/count.txt
done
concatenate org variants
for k in {1..22}; do 
 echo $output/hgdporg_chr$k.vcf.gz >> $output/list
  done

######choose ORG variants
oute110="$maindir/WES/e112"
ORGlist="$maindir/WES/bedfile/e112.list"

passsnp="$maindir/WES/postVC/postVQSR/annpass-OA_snp.vcf.gz"
passindel="$maindir/WES/postVC/postVQSR/annpass-OA_indel.vcf.gz"
passbn="$(basename $passsnp | cut -d "_" -f 1)"

$gatk SelectVariants -V $passsnp  -O $output/ORG-"$passbn"_snp.vcf.gz --L $ORGlist
$gatk SelectVariants -V $passindel  -O $output/ORG-"$passbn"_indel.vcf.gz --L $ORGlist


######Annotate vcf files
pop="$maindir/WES/postVC/postVQSR/bcf-pop.txt"
region="$maindir/WES/postVC/postVQSR/bcf-region.txt"
for file in $(ls ORG-annpass-OA_indel.vcf.gz ORG-annpass-OA_snp.vcf.gz)
do
bn="$(basename $file | cut -d "." -f 1)"
bcftools view $file -S $pop -Ou | bcftools +fill-tags -o bcf-"$bn"  -- -S $region -t AN,AC,AF
echo "CHROM POS ID REF ALT AN AF AC AN_JH AF_JH AC_JH AN_TM AF_TM AC_TM" | sed 's/ /\t/g' >> bcfann-"$bn".txt
bcftools query --format  '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AN\t%INFO/AF\t%INFO/AC\t%INFO/AN_JH\t%INFO/AF_JH\t%INFO/AC_JH\t%INFO/AN_TM\t%INFO/AF_TM\t%INFO/AC_TM\n' bcf-"$bn" >> bcfann-"$bn".txt
done
for file in $(ls hgdporg-auto.vcf.gz hgdporg_chrX.vcf.gz); do
        bcftools view -S $maindir/WES/postVC/e110post/hgdp/hgdpsample-king $file -o bcf-"$(basename $file | cut -d "." -f 1)" --force-samples
        bcftools +fill-tags bcf-"$(basename $file | cut -d "." -f 1)" -o bcfann-"$(basename $file | cut -d "." -f 1)"  -- -S $maindir/hgdpdec20/2022hgdp/vcf_org/kingfilter-afrnon-sample-region -t AN,AC,AF  
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


king="$maindir/WES/postVC/e110post/hgdp/hgdpsample-king"

$plink2 --set-all-var-ids @:# --hardy  'midp'  --keep $king \
        --out $output/hgdporg-auto-homhet --make-pgen --update-sex $output/sexhgdp.txt \
        --vcf $output/hgdporg-auto.vcf.gz
$plink2 --set-all-var-ids @:# --hardy  'midp'  --keep $king \
        --out $output/hgdporg-chrx-homhet --make-pgen --update-sex $output/sexhgdp.txt \
        --vcf $output/hgdporg_chrx.vcf.gz

# #plink geno/hwe, fst for OA ORG
$plink2 --set-all-var-ids @:# --hardy  'midp'  --remove $maindir/WES/postVC/postVQSR/king1table.king.cutoff.out.id \
        --out $output/OA-homhet-snp --make-pgen --vcf $output/ORG-annpass-OA_snp.vcf.gz \
        --update-sex $maindir/WES/postVC/e110post/sexf.txt
$plink2 --set-all-var-ids @:# --hardy   'midp'   --remove $maindir/WES/postVC/postVQSR/king1table.king.cutoff.out.id \
        --out $output/OA-homhet-indel --make-pgen --vcf $output/ORG-annpass-OA_indel.vcf.gz \
        --update-sex $maindir/postVC/e110post/sexf.txt

$gatk SelectVariants -V $passsnp  -O king_annpass-OA_snp.vcf.gz --sample-name $output/bcf-pop.args
$gatk SelectVariants -V $passindel  -O king_annpass-OA_indel.vcf.gz --sample-name $output/bcf-pop.args

varld="$output/LD/varID"
afr="$output/LD/hgdpafr98.txt"
ea="$output/LD/hgdpea220.txt"
bcftools view -r chr1:248000000-248800000  $hgdpin/hgdp_wgs.20190516.full.chr1.vcf.gz -Oz -o $output/LD/chr1:248000000-248800000.vcf.gz

$plink  --vcf  $output/LD/chr1:248000000-248800000.vcf.gz --keep $afr --r2 --ld-snp-list $varld --chr chr1 --ld-window-kb 200 --ld-window-r2 0.1 \
	--out $output/LD/afr-chr1 --update-sex $output/LD/hgdpsex.txt --set-missing-var-ids @:#\$1,\$2
$plink  --vcf  $output/LD/chr1:248000000-248800000.vcf.gz --keep $ea --r2 --ld-snp-list $varld --chr chr1 --ld-window-kb 200 --ld-window-r2 0.1 \
	--out $output/LD/ea-chr1 --update-sex $output/LD/hgdpsex.txt --set-missing-var-ids @:#\$1,\$2

$plink  --vcf  $output/king_annpass-OA-all.vcf.gz --r2 --ld-snp rs9330305  --chr chr1 --ld-window-kb 200 --ld-window-r2 0.1 \
	--out $output/LD/oa-chr1 --update-sex $output/LD/ld-oasex.txt --set-missing-var-ids @:#\$1,\$2
$plink  --vcf  $output/king_annpass-OA-all.vcf.gz --r2 --ld-snp-list varID  --chr chr9 --ld-window-kb 200 --ld-window-r2 0.1 \
 	--out $output/LD/oa-chr9 --update-sex $output/LD/ld-oasex.txt --set-missing-var-ids @:#\$1,\$2


# bcftools view -r chr11:4000000-6000000  $hgdpin/hgdp_wgs.20190516.full.chr11.vcf.gz -Oz -o $output/LD/chr11:4000000-6000000.vcf.gz
LDblock="$maindir/softwares/LDBlockShow/bin"
$LDblock/./LDBlockShow -InVCF $output/king_annpass-OA-all.vcf.gz -OutPut $output/LD/oa-chr9block -BlockType 1 -Region  chr9:122615000-122670000 -OutPdf -SeleVar 2
$LDblock/./ShowLDSVG -InPreFix  $output/LD/oa-chr9block   -OutPut $output/LD/oa-chr9blockshow   -ShowNum -PointSize 3 -OutPng -ResizeH 8000 -SpeSNPName $output/LD/rsid.txt

