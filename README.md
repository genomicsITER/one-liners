# one-liners
Frequently used commands in Bioinformatics
<br>
<br>
#  Select region from 1KGP chr.vcf.gz, extract, compress, and tabix
region="3:10539482-11710471"
<br>
infile="ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
<br>
outfile="1KGP.chr3.region"
<br>
tabix -h ${infile} ${region} > ${outfile}.vcf
<br>
bgzip -c ${outfile}.vcf > ${outfile}.vcf.gz
<br>
tabix -p vcf ${outfile}.vcf.gz
<br>
