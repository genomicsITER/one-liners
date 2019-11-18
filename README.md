# one-liners, frequently used commands in Bioinformatics

**Select region from 1KGP chr.vcf.gz, extract, compress, and tabix**
```Bash
#!/bin/bash
region="3:10000-11000"
infile="ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
outfile="1KGP.chr3.region"
tabix -h ${infile} ${region} > ${outfile}.vcf
bgzip -c ${outfile}.vcf > ${outfile}.vcf.gz
tabix -p vcf ${outfile}.vcf.gz
# End of script
```


**Remove "chr" string in variants in a VCF**
```Bash
#!/bin/bash
awk '{gsub(/^chr/,""); print}' infile.vcf > infile.no_chr.vcf
# End of script
```


**Add "chr" string in variants in VCF**
```Bash
#!/bin/bash
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' infile.no_chr.vcf > infile.vcf
# End of script
```


**Sort karyotipically a VCF (version 1)**
```Bash
#!/bin/bash
##Notice that header must be removed at some point
grep '^#' in.vcf > out.vcf && grep -v '^#' in.vcf | sort -V -k1,1 -k2,2n >> out.vcf
# End of script
```


**Sort karyotypically a VCF (version 2). Use '-V', natural sort of (version) numbers within text**
```Bash
#!/bin/bash
sort -V -k1,1 -k2,2n infile.vcf > outfile.vcf
# End of script
```


**Replace spaces with a single tab**
```Bash
#!/bin/bash
sed 's/ \+/\t/g' infile > outfile
# End of script
```


**Compute BAM coverage with BEDtools**
```Bash
#!/bin/bash
bedtools genomecov -ibam infile.bam -bg > coverage.txt
# End of script
```


**Find duplicated lines in a VCF matching the whole line**
```Bash
#!/bin/bash
awk ' !uniq[$0]++ ' infile.vcf
# End of script
```

**Find duplicated lines in a VCF matching files 1, 2, and 5**
```Bash
#!/bin/bash
awk ' !uniq[$1 FS $2 FS $5]++ ' infile.vcf
# End of script
```

**Find a line by a field on it, delete it, and save the result**
```Bash
#!/bin/bash
string=abc<br>
grep -v $string infile > outfile
# End of script
```


**Count the number of mapped reads for each read in PE reads**
```Bash
#!/bin/bash
samtools view -F 0x40 infile.bam | cut -f1 | sort | uniq | wc -l

#Left read: view the BAM content filtering by the provided SAM flags, cut the first column, sort the data on that column, keep only the uniq data on that column, and count the number of lines.
samtools view -f 0x40 -F 0x4 infile.bam | cut -f1 | sort | uniq | wc -l

#Right read<br>
samtools view -f 0x80 -F 0x4 infile.bam | cut -f1 | sort | uniq  | wc -l

Note: for flags information, see page 5 of https://samtools.github.io/hts-specs/SAMv1.pdf
# End of script
```


**Replace white spaces with tabs**
```Bash
#!/bin/bash
awk -v OFS="\t" '$1=$1' infile > outfile
# End of script
```


**Add rs from INFO field (avsnp150) to ID column in a VCF**
```Bash
#!/bin/bash
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
# End of script
```


**Remove duplicated lines in a file while keeping the original order**
```Bash
#!/bin/bash
awk '!visited[$0]++' infile > deduplicated_infile

#Note: see https://iridakos.com/how-to/2019/05/16/remove-duplicate-lines-preserving-order-linux.html
# Or equivalently (using cat, sort and cut: cat showing numbers, sort using the second column and keep the first ocurrence of number, sort using the number, and cut the second column; clever!):

cat -n infile | sort -uk2 | sort -nk1 | cut -f2-
# End of script
```


**Number of files per extension type in the current directory**
```Bash
#!/bin/bash
nfiletypes () { find . -maxdepth 1 -type f | sed 's/.*\.//' | sort | uniq -c | sed 's/^ *//g' | sed 's/ /\t/g'; }

# Note: run this code and then write "nfiletypes" at the prompt and will see the count of files per extension at the current directory.
# End of script
```


**Parse file with AWK, sum column values in each line, and shows the result**
```Bash
#!/bin/bash
awk -F'[\t]' 'BEGIN{sum=0; OFS="\t"} { for (i=1;i<=NF;i++) a[i]+=$i } END{ for (i in a) print a[i] }' infile
<br>
<br>

**Count genotypes in a VCF imputed at Michigan Imputation Server**br>
Based on script by tommycarstensen<br>
See: https://gatkforums.broadinstitute.org/gatk/discussion/5692/count-genotypes-per-position-in-a-multi-sample-vcf<br>
In each line of your infile, split the line by ":" and the split again by genotype delimiters (unphased, "/"; phased, "|"). Then, split again the genotype field and assign its values (REF, delimiter, and ALT) to GT awk-variable. Then, parse the GT variable according to the expected genotype values (./. for missing or untyped genotype; 0/0 for hom-ref; ref=alt for hom-alt; and heterozygous for the rest of conditions: 1/0 or 0/1). Finally, show the count of genotypes and, possible, the count of variant genotypes (Het+HomAlt).<br>

```Bash
#!/bin/bash
zgrep -v "^#" infile.vcf.gz | awk '{
    unknown=0; homref=0; het=0; homalt=0; for(i=10;i<=NF;i++) {
    split($i,a,":"); split(a[1],GT,"[/|]");
    if(GT[1]=="."&&GT[2]==".") {unknown++}
    else if(GT[1]==0&&GT[2]==0) {homref++}
    else if(GT[1]==GT[2]) {homalt++}
    else {het++}};
    print $1,$2,$3":"$4"-"$5,$4,$5,unknown,homref,het,homalt,het+homalt}' > tmp1

##Replace spaces with tab
awk -v OFS="\t" '$1=$1' tmp1 > tmp2

##Prepare a header
echo -e "#CHR\tPOS\tID\tREF\tALT\tnumGT.Unknowns\tnumGT.HomRef\tnumGT.Het \
     \tnumGT.HomAlt\tnumGT.(Het+HomAlt)" > header<br>
cat header tmp2 > infile.variant-genotypes.counts

##Extract variants that have a genotype count equal or higher than a genotype count threshold
awk -F'[ ]' '{ if ($10 >= 5) print $3 }' infile.variant-genotypes.counts > variant-list

# End of script
```

