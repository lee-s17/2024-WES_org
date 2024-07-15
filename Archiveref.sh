#1-metrics_postfilt_pca
#!/bin/bash
set -e  #exit if any command fails

#####1collectmetrics.script
#setdir
softdir="/home/z6/Downloads/softwares"
rawdata="/media/z6/DATADRIVE1/temuanjahai_wes/m3-backup/szemei/"
maindir="/home/z6/Documents/workingdir/WES/"

#tools
plink2="$softdir/./plink2"
gatk="$softdir/gatk-4.2.3.0/./gatk"

snp="$rawdata/joint_calling/results/gatk_VQSR_snp/OA.vcf.gz"
indel="$rawdata/joint_calling/results/gatk_VQSR_indel/OA.vcf.gz"
grchref="$rawdata/preprocessing/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
dbsnp138="$rawdata/gatk_bundle/Homo_sapiens_assembly38.dbsnp138.vcf"

out="$rawdata/postVC/metrics"

#before VQSR filtering
$gatk CollectVariantCallingMetrics -I $rawdata/joint_calling/results/gatk_selectvariants_snp/OA.vcf.gz \
         --DBSNP $dbsnp138 -O $out/metrics_snpunfilt-dbsnp138 &
$gatk CollectVariantCallingMetrics -I $rawdata/joint_calling/results/gatk_selectvariants_indel/OA.vcf.gz \
         --DBSNP $dbsnp138 -O $out/metrics_indelunfilt-dbsnp138 

#after VQSR filtering
$gatk CollectVariantCallingMetrics -I $snp --DBSNP $dbsnp138 -O $out/metrics_snpVQSR-dbsnp138 &
$gatk CollectVariantCallingMetrics -I $indel --DBSNP $dbsnp138 -O $out/metrics_indelVQSR-dbsnp138


#####2postfilt.script
output="$maindir/WES/postVC/postVQSR"

dbsnp155="$maindir/gatk_bundle/Human_GRCh38.dbSNP155.vcf.gz"
ref="$rawdata/preprocessing/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

source ~/anaconda3/etc/profile.d/conda.sh
conda activate vcf
#collect variants passed VQSR
bcftools view -f PASS $snp | bgzip -c > $output/pass-OA_snp.vcf.gz
bcftools view -f PASS $indel | bgzip -c > $output/pass-OA_indel.vcf.gz
tabix -p vcf $output/pass-OA_snp.vcf.gz
tabix -p vcf $output/pass-OA_indel.vcf.gz
#annotate using dbsnp155
gatk VariantAnnotator -R $ref -V $output/pass-OA_snp.vcf.gz --dbsnp $dbsnp155 -O $output/annpass-OA_snp.vcf.gz
gatk VariantAnnotator -R $ref -V $output/pass-OA_indel.vcf.gz --dbsnp $dbsnp155 -O $output/annpass-OA_indel.vcf.gz

#main file for downstream
passsnp="$output/annpass-OA_snp.vcf.gz"
passindel="$output/annpass-OA_indel.vcf.gz"
passbn="$(basename $passsnp | cut -d "_" -f 1)"

#biallelic
$gatk SelectVariants -V $passsnp  -O $output/"$passbn"_biallsnp.vcf.gz -select-type-to-include SNP --restrict-alleles-to BIALLELIC
$gatk SelectVariants -V $passindel  -O $output/"$passbn"_bindel.vcf.gz -select-type-to-include INDEL --restrict-alleles-to BIALLELIC

#filter by hwe/geno; familial relationship (non-autosomal will be excluded by king and plink pca tools themselves)
$plink2 --set-all-var-ids @:# --hwe 10e-6 --geno 0.05 --out $output/"$passbn"-genohwebiallsnp --make-bed \
        --vcf $output/"$passbn"_biallsnp.vcf.gz
$plink2 --king-table-filter 0.177 --bfile $output/"$passbn"-genohwebiallsnp --out $output/king1table --make-king-table
$plink2 --king-cutoff 0.177 --bfile $output/"$passbn"-genohwebiallsnp --out $output/king1table

#define file for king-removed individuals
king="$output/king1table.king.cutoff.out.id"

#filter geno, hwe, king- make bed
$plink2 --set-all-var-ids @:# --hwe 10e-6 --geno 0.05 --remove $king \
        --out $output/"$passbn"-genohweking-biallsnp --make-bed \
        --vcf $output/"$passbn"_biallsnp.vcf.gz
$plink2 --set-all-var-ids @:# --hwe 10e-6 --geno 0.05 --remove $king \
        --out $output/"$passbn"-genohweking-bindel --make-bed \
        --vcf $output/"$passbn"_bindel.vcf.gz



#####3PCA-script
mkdir $maindir/postVC/pca
cd $maindir/postVC/pca

#use bcftools to select variants in OA (pca-OA-genohweking-biallsnp.list) from references
gatk SelectVariants --L pca-OA-genohweking-biallsnp.list -V /media/z6/DATADRIVE1/VCF-NAS/hg38_vcf/full/oa-biallasnps-auto.vcf.gz  -O pca-OAWGS.vcf.gz
gatk SelectVariants --L pca-OA-genohweking-biallsnp.list -V /media/z6/DATADRIVE1/rawhgdp/hgdpwhole_biallsnps/biallsnps_autosomes.vcf.gz  -O pca-hgdp.vcf.gz
bcftools annotate -x  FORMAT/AD pca-hgdp.vcf.gz | bgzip -c > pca-annhgdp.vcf.gz 
tabix -p vcf pca-annhgdp.vcf.gz 
gatk SelectVariants --L pca-OA-genohweking-biallsnp.list -V $maindir/postVC/pca/Andamanese/Andamanese-hg38.vcf.gz  -O pca-Andaman.vcf.gz

gatk SelectVariants --L pca-OA-genohweking-biallsnp.list -V $maindir/postVC/postVQSR/annpass-OA_biallsnp.vcf.gz  -O pca-OAWES.vcf.gz
#merge ref and OA
bcftools merge -Oz -o merged-allwAndaman.vcf.gz  pca-annhgdp.vcf.gz pca-OAWGS.vcf.gz pca-OAWES.vcf.gz pca-Andaman.vcf.gz
tabix -p vcf merged-allwAndaman.vcf.gz
gatk SelectVariants -V merged-allwAndaman.vcf.gz  -O merged-allwAndaman_biallsnp.vcf.gz -select-type-to-include SNP --restrict-alleles-to BIALLELIC

$plink2 --vcf merged-allwAndaman.vcf.gz  --set-all-var-ids @:#  --out merged-allwAndaman-genohweking-mafld-biallsnp --hwe 10e-6 --geno 0.05 \
        --remove /media/z6/DATADRIVE1/rawhgdp/hgdpwhole_biallsnps/kingthird.king.cutoff.out.id --maf 0.01 --indep-pairwise 50 5 0.2
$plink2 -vcf merged-allwAndaman_biallsnp.vcf.gz --set-all-var-ids @:#  -remove /media/z6/DATADRIVE1/rawhgdp/hgdpwhole_biallsnps/kingthird.king.cutoff.out.id \
        --extract merged-allwAndaman-genohweking-mafld-biallsnp.prune.in --out pruned-merged-allwAndaman-genohweking-mafld-biallsnp --make-bed

#num=$(wc -l < pruned-merged-genohweking-mafld-biallsnp.fam)
$plink2 --out pca-merged-allwAndaman-genohweking-mafld-biallsnp --pca $num --set-all-var-ids @:#  --bfile pruned-merged-allwAndaman-genohweking-mafld-biallsnp



#2-ORG_vep
######2postfilt.script-select SNP in ORG
set -e  #exit if any command fails
#tools
plink2="/home/z6/Documents/workingdir/softwares/./plink2"
gatk="/home/z6/Documents/workingdir/softwares/gatk-4.4.0.0/./gatk"

#intersect Agilent WES and 1000GP strict mask
cd /home/z6/Documents/workingdir/WES/bedfile
bedtools intersect  -a Agilentv6r2_hg38/S07604514_Regions.bed -b 20160622_genome_mask_GRCh38/StrictMask/20160622.allChr.mask.bed > intersect-agilent_strict1000gp.bed
#intersect Agilent WES-1000GP Strict Mask and ORG
cut -f 7,8,9,2,7,12,15 e112_masterlist.txt | awk '{print "chr"$2 "\t"  $3 "\t" $4 "\t" $1 "\t" $5 "\t" $6 }' | sort -V -k 1 > e112.bed
bedtools intersect -a e112.bed -b intersect-agilent_strict1000gp.bed -wao > e112_agilent-strict1000gp.txt
#combine total coverage for each ORG
cut -f 4,10 e112_agilent-strict1000gp.txt | awk 'BEGIN{FS=OFS="\t"}{a[$1]+=$2}END{for(i in a) print i,a[i]}'


######choose ORG variants
oute110="/home/z6/Documents/workingdir/WES/postVC/e110post"
ORGlist="/home/z6/Documents/workingdir/WES/bedfile/e110/e110final.list"

source ~/anaconda3/etc/profile.d/conda.sh
conda activate gatk

passsnp="/home/z6/Documents/workingdir/WES/postVC/postVQSR/annpass-OA_snp.vcf.gz"
passindel="/home/z6/Documents/workingdir/WES/postVC/postVQSR/annpass-OA_indel.vcf.gz"
passbn="$(basename $passsnp | cut -d "_" -f 1)"

$gatk SelectVariants -V $passsnp  -O $oute110/ORG-"$passbn"_snp.vcf.gz --L $ORGlist
$gatk SelectVariants -V $passindel  -O $oute110/ORG-"$passbn"_indel.vcf.gz --L $ORGlist
#biallelic
$gatk SelectVariants -V $oute110/ORG-"$passbn"_snp.vcf.gz  -O $oute110/ORG-"$passbn"_biallsnp.vcf.gz -select-type-to-include SNP --restrict-alleles-to BIALLELIC
#no multiallelic indels
$gatk SelectVariants -V $oute110/ORG-"$passbn"_indel.vcf.gz  -O $oute110/ORG-"$passbn"_bindel.vcf.gz -select-type-to-include INDEL --restrict-alleles-to BIALLELIC

#plink filter geno/ hwe
$plink2 --set-all-var-ids @:# --hwe 10e-6 --geno 0.05 \
	--remove /home/z6/Documents/workingdir/WES/postVC/postVQSR/king1table.king.cutoff.out.id \
    --out $oute110/ORG-annpass-genohweking_biallsnp --make-bed --vcf $oute110/ORG-"$passbn"_biallsnp.vcf.gz \
    --update-sex /home/z6/Documents/workingdir/WES/postVC/e110post/sexf.txt
$plink2 --set-all-var-ids @:# --hwe 10e-6 --geno 0.05 \
	--remove /home/z6/Documents/workingdir/WES/postVC/postVQSR/king1table.king.cutoff.out.id \
    --out $oute110/ORG-annpass-genohweking_indel --make-bed --vcf $oute110/ORG-"$passbn"_bindel.vcf.gz \
    --update-sex /home/z6/Documents/workingdir/WES/postVC/e110post/sexf.txt
conda deactivate


######annotate vcf files
pop="/home/z6/Documents/workingdir/WES/postVC/postVQSR/bcf-pop.txt"
region="/home/z6/Documents/workingdir/WES/postVC/postVQSR/bcf-region.txt"
source /home/z6/mambaforge/etc/profile.d/conda.sh 
conda activate bcftools

for file in $(ls $oute110/ORG-"$passbn"_snp.vcf.gz $oute110/ORG-"$passbn"_indel.vcf.gz)
do
bn="$(basename $file | cut -d "." -f 1)"
bcftools view $file -S $pop -Ou | bcftools +fill-tags -o $oute110/bcf-"$bn"  -- -S $region -t AN,AC,AF
echo "CHROM POS ID REF ALT AN AF AC AN_JH AF_JH AC_JH AN_TM AF_TM AC_TM" | sed 's/ /\t/g' >> $oute110/bcfann-"$bn".txt
bcftools query --format  '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AN\t%INFO/AF\t%INFO/AC\t%INFO/AN_JH\t%INFO/AF_JH\t%INFO/AC_JH\t%INFO/AN_TM\t%INFO/AF_TM\t%INFO/AC_TM\n' $oute110/bcf-"$bn" >> $oute110/bcfann-"$bn".txt
done

conda deactivate


######vep
for vcf in $(ls $oute110/ORG-"$passbn"_snp.vcf.gz $oute110/ORG-"$passbn"_indel.vcf.gz)
do 
        bn="$(basename $vcf)";

/home/z6/Documents/workingdir/softwares/ensembl-vep/./vep --cache --offline --format vcf --tab --force_overwrite --mane \
--check_existing \
--polyphen b --sift b \
--dir_plugins  /home/z6/.vep/Plugins/ \
--fasta  /home/z6/.vep/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz --input_file  $vcf \
--output_file $oute110/vep-$bn  --everything \
--plugin  AncestralAllele,/home/z6/Documents/workingdir/softwares/ensembl-vep/cache/homo_sapiens_ancestor_GRCh38.fa.gz \
--plugin CADD,/home/z6/Documents/workingdir/softwares/ensembl-vep/cache/whole_genome_SNVs.tsv.gz,/home/z6/Documents/workingdir/softwares/ensembl-vep/cache/gnomad.genomes.r3.0.indel.tsv.gz;

        grep "#" $oute110/vep-$bn | cut --complement -f 8-12,26-34,39-43,50-69,76-80 | tail -1 >> $oute110/vep-cut-$bn;
        grep -v "#" $oute110/vep-$bn | cut --complement -f 8-12,26-34,39-43,50-69,76-80 >> $oute110/vep-cut-$bn;
done
#Match ENSG
awk -F'\t' 'NR==FNR{c[$2]++;next};c[$4]' /home/z6/Documents/workingdir/WES/bedfile/e110/org-e110-master.list $oute110/vep-cut-ORG-annpass-OA_indel.vcf.gz > $oute110/matchensg-vep-cut-ORG-annpass-OA_indel
awk -F'\t' 'NR==FNR{c[$2]++;next};c[$4]' /home/z6/Documents/workingdir/WES/bedfile/e110/org-e110-master.list $oute110/vep-cut-ORG-annpass-OA_snp.vcf.gz > $oute110/matchensg-vep-cut-ORG-annpass-OA_snp



#2A-hgdp
#!/bin/bash

input="/media/z6/DATADRIVE1/rawhgdp"
output="/home/z6/Documents/workingdir/WES/postVC/e110post/hgdp"

##call gatk
plink2="/home/z6/Documents/workingdir/softwares/./plink2"
gatk="/home/z6/Documents/workingdir/softwares/gatk-4.4.0.0/./gatk"
source ~/anaconda3/etc/profile.d/conda.sh
conda activate gatk

#######Input: select variants in ORG
for file in  $(ls $input*.vcf.gz); do
        bn="$(basename "$file" | cut -c24-)"
        $gatk SelectVariants -O $output/org_$bn  -V $file --L /home/z6/Documents/workingdir/WES/bedfile/e110/e110final.list
done

for k in {1..22}; do
  echo $input/org_chr$k.vcf.gz >> $output/count.txt
  zgrep -v "#" $input/org_chr$k.vcf.gz | wc -l >> $output/count.txt
done
#concatenate org variants
for k in {1..22}; do 
  echo $output/org_chr$k.vcf.gz >> $output/list
  done 

conda deactivate
source /home/z6/mambaforge/etc/profile.d/conda.sh 
conda activate bcftools

str=$(less $output/list)
bcftools concat  $str -Oz -o $output/org-auto.vcf.gz
tabix -p vcf $output/org-auto.vcf.gz



#######Select variant in ORG based on types
        $gatk SelectVariants -V $output/org-auto.vcf.gz -O $output/org-snps_auto.vcf.gz --select-type-to-include SNP
        $gatk SelectVariants -V $output/org-auto.vcf.gz -O $output/org-biallsnps_auto.vcf.gz --select-type-to-include SNP --restrict-alleles-to BIALLELIC
        $gatk SelectVariants -V $output/org-auto.vcf.gz -O $output/org-bindels_auto.vcf.gz --select-type-to-include INDEL --restrict-alleles-to BIALLELIC
        $gatk SelectVariants -V $output/org-auto.vcf.gz -O $output/org-indels_auto.vcf.gz --select-type-to-include INDEL 
        $gatk SelectVariants -V $output/org-auto.vcf.gz -O $output/org-mixed_auto.vcf.gz --select-type-to-include MIXED
        
for t in $(ls $output/org-*_auto.vcf.gz); do 
  echo $t >> $output/count.txt
  zgrep -v "#" $t | wc -l >> $output/count.txt
done

cd $output
for file in $(ls org-auto.vcf.gz org_chrX.vcf.gz); do
        bcftools view -S hgdpsample-king $file -o bcf-hgdp"$(basename $file | cut -d "." -f 1)" --force-samples
        bcftools +fill-tags bcf-hgdp"$(basename $file | cut -d "." -f 1)" -o bcfann-hgdp"$(basename $file | cut -d "." -f 1)"  -- -S /home/z6/Documents/workingdir/hgdpdec20/2022hgdp/vcf_org/kingfilter-afrnon-sample-region -t AN,AC,AF  
        echo "CHROM POS ID REF ALT AA_chimp AN AF AC AN_AFR AF_AFR AC_AFR AN_NONAFR AF_NONAFR AC_NONAFR" | sed 's/ /\t/g' >> bcfannfreq-hgdp"$(basename $file | cut -d "." -f 1)" 
        bcftools query --format  \
        '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AA_chimp\t%INFO/AN\t%INFO/AF\t%INFO/AC\t%INFO/AN_AFR\t%INFO/AF_AFR\t%INFO/AC_AFR\t%INFO/AN_NONAFR\t%INFO/AF_NONAFR\t%INFO/AC_NONAFR\n' \
        bcfann-hgdp"$(basename $file | cut -d "." -f 1)" >> bcfannfreq-hgdp"$(basename $file | cut -d "." -f 1)";
done





