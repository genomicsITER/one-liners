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
**Count the number of mapped reads for each read in PE reads**<br>
samtools view -F 0x40 infile.bam | cut -f1 | sort | uniq | wc -l<br>
#Left read: view the BAM content filtering by the provided SAM flags, cut the first column, sort the data on that column, keep only the uniq data on that column, and count the number of lines.<br>
samtools view -f 0x40 -F 0x4 infile.bam | cut -f1 | sort | uniq | wc -l<br>
#Right read<br>
samtools view -f 0x80 -F 0x4 infile.bam | cut -f1 | sort | uniq  | wc -l<br>
<br>Note: for flags information, see page 5 of https://samtools.github.io/hts-specs/SAMv1.pdf
<br>
<br>
**Replace white spaces with tabs**<br>
awk -v OFS="\t" '$1=$1' infile > outfile
<br>
<br>
**Add rs from INFO field (avsnp150) to ID column in a VCF**<br>
cat infile | awk '
BEGIN { FS=OFS="\t" }
{
if ($3 == ".") {
$3=match($8, /avsnp150=((rs[0-9]+)|.)/, arr)
$3=arr[0]
sub(/avsnp150=/, "", $3)
}
print
}' > outfile
<br>
<br>
**Remove duplicated lines in a file while keeping the original order**<br>
awk '!visited[$0]++' infile > deduplicated_infile
<br>
<br>Note: see https://iridakos.com/how-to/2019/05/16/remove-duplicate-lines-preserving-order-linux.html
<br>
<br>Or equivalently (using cat, sort and cut: cat showing numbers, sort using the second column and keep the first ocurrence of number, sort using the number, and cut the second column; clever!):
<br>
cat -n infile | sort -uk2 | sort -nk1 | cut -f2-
<br>
<br>
**Number of files per extension type in the current directory**<br>
nfiletypes () { find . -maxdepth 1 -type f | sed 's/.*\.//' | sort | uniq -c | sed 's/^ *//g' | sed 's/ /\t/g'; }
<br>
<br>Note: run this code and then write "nfiletypes" at the prompt and will see the count of files per extension at the current directory.
<br>
<br>
**Parse file with AWK, sum column values in each line, and shows the result**<br>
awk -F'[\t]' 'BEGIN{sum=0; OFS="\t"} { for (i=1;i<=NF;i++) a[i]+=$i } END{ for (i in a) print a[i] }' infile
<br>
<br>

**Count genotypes in a VCF imputed at Michigan Imputation Server**br>
Based on script by tommycarstensen<br>
See: https://gatkforums.broadinstitute.org/gatk/discussion/5692/count-genotypes-per-position-in-a-multi-sample-vcf<br>
<br>
zgrep -v "^#" infile.vcf.gz | awk '{<br>
    unknown=0; homref=0; het=0; homalt=0; for(i=10;i<=NF;i++) {<br>
    split($i,a,":"); split(a[1],GT,"[/|]");<br>
    if(GT[1]=="."&&GT[2]==".") {unknown++}<br>
    else if(GT[1]==0&&GT[2]==0) {homref++}<br>
    else if(GT[1]==GT[2]) {homalt++}<br>
    else {het++}};<br>
    print $1,$2,$3":"$4"-"$5,$4,$5,unknown,homref,het,homalt,het+homalt}' > tmp1<br>
<br>
##Replace spaces with tab<br>
awk -v OFS="\t" '$1=$1' tmp1 > tmp2<br>
<br>
##Prepare a header<br>
echo -e "#CHR\tPOS\tID\tREF\tALT\tnumGT.Unknowns\tnumGT.HomRef\tnumGT.Het\tnumGT.HomAlt\tnumGT.(Het+HomAlt)" > header<br>
cat header tmp2 > infile.variant-genotypes.counts<br>
<br>
##Clean the house<br>
rm tmp1 tmp2<br>
<br>
##Extract variants that have a genotype count equal or higher than a genotype count threshold<br>
awk -F'[ ]' '{ if ($10 >= 5) print $3 }' infile.variant-genotypes.counts > variant-list<br>
<br>
<br>

