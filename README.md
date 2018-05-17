# one-liners
**Frequently used commands in Bioinformatics**
<br>
<br>
<hr>
**Select region from 1KGP chr.vcf.gz, extract, compress, and tabix**
region="3:10000-11000"
infile="ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
outfile="1KGP.chr3.region"
tabix -h ${infile} ${region} > ${outfile}.vcf
bgzip -c ${outfile}.vcf > ${outfile}.vcf.gz
tabix -p vcf ${outfile}.vcf.gz
<br>
<br>
<hr>
**Select a subset of samples (list here contains GBR individuals from 1KGP), and tabix**
outfile="1KGP.chr3"
lista="pop_GBR"
suffix="GBR"
vcf-subset -c ${lista} ${outfile}.vcf.gz | bgzip -c > ${outfile}.${suffix}.vcf.gz
tabix -p vcf ${outfile}.${suffix}.vcf.gz
<br>
<br>
<hr>
**Remove "chr" string in variants in a VCF**
awk '{gsub(/^chr/,""); print}' infile.vcf > infile.no_chr.vcf
<br>
<br>
<hr>
**Add "chr" string in variants in VCF**
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' infile.no_chr.vcf > infile.vcf
<br>
<br>
<hr>
**Sort karyotipically a VCF (version 1)**
sort -K1,1 -k2,2n infile.vcf > outfile.vcf
<br>
<br>
<hr>
**Sort karyotypically a VCF (version 2). Use '-V', natural sort of (version) numbers within text**
sort -V -k1,1 -k2,2n infile.vcf > outfile.vcf


