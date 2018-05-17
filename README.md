# one-liners
Frequently used commands in Bioinformatics
<br>
<br>
##3  Select region from 1KGP chr.vcf.gz, extract, compress, and tabix
----------
region="3:10000-11000"
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
# Select a subset of samples (list here contains GBR individuals from 1KGP), and tabix
----------
outfile="1KGP.chr3"
lista="pop_GBR"
suffix="GBR"
vcf-subset -c ${lista} ${outfile}.vcf.gz | bgzip -c > ${outfile}.${suffix}.vcf.gz
tabix -p vcf ${outfile}.${suffix}.vcf.gz
<br>
<br>
# Remove "chr" string in variants in a VCF
----------
awk '{gsub(/^chr/,""); print}' infile.vcf > infile.no_chr.vcf
<br>
<br>
# Add "chr" string in variants in VCF
----------
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' infile.no_chr.vcf > infile.vcf
<br>
<br>
*Sort karyotipically a VCF (version 1)*
----------
sort -K1,1 -k2,2n infile.vcf > outfile.vcf
<br>
<br>
*Sort karyotypically a VCF (version 2). Use '-V', natural sort of (version) numbers within text:*
----------
sort -V -k1,1 -k2,2n infile.vcf > outfile.vcf
<br>
<br>

