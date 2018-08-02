# one-liners, frequently used commands in Bioinformatics

**Select region from 1KGP chr.vcf.gz, extract, compress, and tabix**<br>
region="3:10000-11000"<br>
infile="ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"<br>
outfile="1KGP.chr3.region"<br>
tabix -h ${infile} ${region} > ${outfile}.vcf<br>
bgzip -c ${outfile}.vcf > ${outfile}.vcf.gz<br>
tabix -p vcf ${outfile}.vcf.gz
<br>
<br>
**Select a subset of samples (list here contains GBR individuals from 1KGP), and tabix**<br>
outfile="1KGP.chr3"<br>
lista="pop_GBR"<br>
suffix="GBR"<br>
vcf-subset -c ${lista} ${outfile}.vcf.gz | bgzip -c > ${outfile}.${suffix}.vcf.gz<br>
tabix -p vcf ${outfile}.${suffix}.vcf.gz<br>
<br>
<br>
**Remove "chr" string in variants in a VCF**<br>
awk '{gsub(/^chr/,""); print}' infile.vcf > infile.no_chr.vcf<br>
<br>
<br>
**Add "chr" string in variants in VCF**<br>
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' infile.no_chr.vcf > infile.vcf<br>
<br>
<br>
**Sort karyotipically a VCF (version 1)**<br>
##Notice that header must be removed at some point <br>
grep '^#' in.vcf > out.vcf && grep -v '^#' in.vcf | sort -V -k1,1 -k2,2n >> out.vcf
<br>
<br>
**Sort karyotypically a VCF (version 2). Use '-V', natural sort of (version) numbers within text**<br>
sort -V -k1,1 -k2,2n infile.vcf > outfile.vcf<br>
<br>
<br>
**Replace spaces with a single tab**<br>
sed 's/ \+/\t/g' infile > outfile<br>
<br>
<br>
**Compute BAM coverage with BEDtools**<br>
bedtools genomecov -ibam infile.bam -bg > coverage.txt<br>
<br>
<br>
**Find duplicated lines in a VCF matching the whole line**<br>
awk ' !uniq[$0]++ ' infile.vcf<br>
<br>
<br>
**Find duplicated lines in a VCF matching files 1, 2, and 5**<br>
awk ' !uniq[$1 FS $2 FS $5]++ ' infile.vcf<br>
<br>
<br>
**Find a line by a field on it, delete it, and save the result**<br>
string=abc<br>
grep -v $string infile > outfile<br>
<br>
<br>

