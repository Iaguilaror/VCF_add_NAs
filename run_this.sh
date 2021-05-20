mkdir -p results
bcftools view -h data/out.recode.vcf > results/header.tmp
bcftools view -H data/out.recode.vcf > results/body.tmp
# Luego se corre el script en R
Rscript --vanilla rscript.R
cat results/header.tmp results/body_with_missingGT.tmp | bcftools +fill-tags > results/out.recode.addedmissingGT.vcf
bgzip results/out.recode.addedmissingGT.vcf
tabix -p vcf results/out.recode.addedmissingGT.vcf.gz
rm results/*.tmp
