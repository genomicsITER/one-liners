# one-liners, frequently used commands in Bioinformatics

**Markdown References**

+ https://guides.github.com/features/mastering-markdown/

+ https://help.github.com/en/github/writing-on-github/basic-writing-and-formatting-syntax


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


**Sort karyotypically a VCF (version 3): using vcf-sort (from vcftools)**
```Bash
#!/bin/bash

cat infile.vcf | vcf-sort --chromosomal-order > infile.sorted.vcf

# End of script
```


**Sort karyotypically a VCF (version 4): using PICARD**
```Bash
#!/bin/bash

java -jar picard.jar SortVcf I=unsorted.infile.vcf O=sorted.infile.vcf

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

string=abc
grep -v $string infile > outfile

# End of script
```


**Count the number of mapped reads for each read in PE reads**
```Bash
#!/bin/bash

samtools view -F 0x40 infile.bam | cut -f1 | sort | uniq | wc -l

#Left read: view the BAM content filtering by the provided SAM flags, cut the first column, 
# sort the data on that column, keep only the uniq data on that column, and count the number of lines.
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
# Or equivalently (using cat, sort and cut: cat showing numbers, sort using the second column and 
# keep the first ocurrence of number, sort using the number, and cut the second column; clever!):

cat -n infile | sort -uk2 | sort -nk1 | cut -f2-

# End of script
```


**Number of files per extension type in the current directory**
```Bash
#!/bin/bash

nfiletypes () { find . -maxdepth 1 -type f | sed 's/.*\.//' | sort | uniq -c \
    | sed 's/^ *//g' | sed 's/ /\t/g'; }

# Note: run this code and then write "nfiletypes" at the prompt and will see the count of 
# files per extension at the current directory.

# End of script
```


**Parse file with AWK, sum column values in each line, and shows the result**
```Bash
#!/bin/bash

awk -F'[\t]' 'BEGIN{sum=0; OFS="\t"} { for (i=1;i<=NF;i++) a[i]+=$i } END{ for (i in a) print a[i] }' infile

# End of script
```


**Count genotypes in a VCF imputed at Michigan Imputation Server**

Based on script by tommycarstensen

See: https://gatkforums.broadinstitute.org/gatk/discussion/5692/count-genotypes-per-position-in-a-multi-sample-vcf

In each line of your infile, split the line by ":" and the split again by genotype delimiters 
(unphased, "/"; phased, "|"). Then, split again the genotype field and assign its values 
(REF, delimiter, and ALT) to GT awk-variable. Then, parse the GT variable according to the 
expected genotype values (./. for missing or untyped genotype; 0/0 for hom-ref; ref=alt for 
hom-alt; and heterozygous for the rest of conditions: 1/0 or 0/1). Finally, show the count 
of genotypes and, possible, the count of variant genotypes (Het+HomAlt).

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
    
# End of script
```


**Replace spaces with tab**
```Bash
#!/bin/bash

awk -v OFS="\t" '$1=$1' tmp1 > tmp2

# End of script
```


**Prepare a header**
```Bash
#!/bin/bash

echo -e "#CHR\tPOS\tID\tREF\tALT\tnumGT.Unknowns\tnumGT.HomRef\tnumGT.Het \
     \tnumGT.HomAlt\tnumGT.(Het+HomAlt)" > header<br>
cat header tmp2 > infile.variant-genotypes.counts

# End of script
```


**Extract variants that have a genotype count equal or higher than a genotype count threshold**
```Bash
#!/bin/bash

awk -F'[ ]' '{ if ($10 >= 5) print $3 }' infile.variant-genotypes.counts > variant-list

# End of script
```

**Count the lines of each file in a dir and get the sum of all lines**

Credits: http://stackoverflow.com/questions/13727917/ddg#13728131

```Bash
#!/bin/bash

find ${indir} -type f -name "*.selected-files" -print0 | wc -l --files0-from=-

# End of script
```


**Count the total number of lines in the selected files of a dir**

Credits: http://stackoverflow.com/questions/13727917/ddg#13728131

```Bash
#!/bin/bash

find ${indir} -type f -name "*.selected-files" -exec cat {} + | wc -l

# End of script
```


**Grab the body of a file excluding the header**
```Bash
#!/bin/bash

tail -n +2 ${infile} > file-without-header

# End of script
```


**Count number of lines in a file**
```Bash
#!/bin/bash

wc -l ${infile}

# End of script
```


**Count number of columns in a file**
```Bash
#!/bin/bash

head -n 1 ${infile} | awk '{print NF}'
#Or
awk '{print NF; exit}' ${infile}
# End of script
```


**Replace many spaces with single tabs, and specially the leading spaces of a PLINK '*.frq' file**
```Bash
#!/bin/bash

sed 's/ \+/\t/g' ${infile} | sed -e 's/^[\t]*//' >${infile}.no-trailing-spaces.tabs

# End of script
```


**Recode (0/0, 0/1, 1/0, 1/1) genotypes into (0,1,2) genotypes**
```Bash
#!/bin/bash

sed 's/0|0/0/g' test | sed 's/0|1/1/g' | sed 's/1|0/1/g' | sed 's/1|1/2/g

# End of script
```


**Convert a VCF into a table of variants combining PERL and GATK4 VariantsToTable**

Credits: https://gatkforums.broadinstitute.org/gatk/profile/ericminikel && https://gatkforums.broadinstitute.org/gatk/profile/dobenshain

```Bash
#!/bin/bash

# Grab INFO tags 
zcat infile.vcf.gz | perl -ne '/^.*INFO.*ID=([a-zA-Z0-9_]+),/ && print "-F $1 \\ "' | uniq
# Example of output (allele frequencies from gnomAD):
-F AC -F AF -F AF_afr -F AF_amr -F AF_asj \

# Grab FORMAT tags
zcat infile.vcf.gz | perl -ne '/^.*FORMAT.*ID=([a-zA-Z0-9_]+),/ && print "-GF $1 \\ "' | uniq
# Example of output:
-GF AD \ -GF DP \ -GF GQ \ -GF GT \ -GF MIN_DP \ -GF PGT \ -GF PID \ -GF PL \ -GF RGQ \ -GF SB \

# Complex example using GATK4 (notice that INFO tags are written after "-F" and FORMAT tags are written after "GF").
# GATK3.8 does not produce a proper header in the output and does not allow for a filtering control as GATK4.

module load java-jre/1.8.0_77
module load gatk/4.1.4.0

gatk VariantsToTable \
--variant ${infile} \
--show-filtered true \
--showHidden true \
--verbosity INFO \
-F CHROM \
-F POS \
-F ID \
-F REF \
-F ALT \
-F QUAL \
-F FILTER \
-F AAChange.refGene -F AC -F AF -F AF_afr -F AF_amr -F AF_asj -F AF_eas -F AF_female -F AF_fin -F AF_male -F AF_nfe -F AF_oth -F AF_popmax \
-F AF_raw -F AF_sas -F AN -F BA1 -F BP1 -F BP2 -F BP3 -F BP4 -F BP5 -F BP6 -F BP7 -F BS1 -F BS2 -F BS3 -F BS4 -F BaseQRankSum -F CADD_phred \
-F CADD_raw -F CADD_raw_rankscore -F CLNALLELEID -F CLNDISDB -F CLNDN -F CLNREVSTAT -F CLNSIG -F ClippingRankSum -F DANN_rankscore \
-F DANN_score -F DB -F DP -F DS -F Eigen-PC-raw -F Eigen-raw -F Eigen_coding_or_noncoding -F ExAC_AFR -F ExAC_ALL -F ExAC_AMR -F ExAC_EAS \
-F ExAC_FIN -F ExAC_NFE -F ExAC_OTH -F ExAC_SAS -F ExcessHet -F ExonicFunc.refGene -F FATHMM_converted_rankscore -F FATHMM_pred -F FATHMM_score \
-F FS -F Func.refGene -F GERP++_RS -F GERP++_RS_rankscore -F GTEx_V6_gene -F GTEx_V6_tissue -F Gene.refGene -F GeneDetail.refGene \
-F GenoCanyon_score -F GenoCanyon_score_rankscore -F HaplotypeScore -F InbreedingCoeff -F InterVar_automated -F Interpro_domain \
-F LRT_converted_rankscore -F LRT_pred -F LRT_score -F M-CAP_pred -F M-CAP_rankscore -F M-CAP_score -F MLEAC -F MLEAF -F MQ -F MQRankSum \
-F MetaLR_pred -F MetaLR_rankscore -F MetaLR_score -F MetaSVM_pred -F MetaSVM_rankscore -F MetaSVM_score -F MutationAssessor_pred \
-F MutationAssessor_score -F MutationAssessor_score_rankscore -F MutationTaster_converted_rankscore -F MutationTaster_pred \
-F MutationTaster_score -F NEGATIVE_TRAIN_SITE -F PM1 -F PM2 -F PM3 -F PM4 -F PM5 -F PM6 -F POSITIVE_TRAIN_SITE -F PP1 -F PP2 -F PP3 \
-F PP4 -F PP5 -F PROVEAN_converted_rankscore -F PROVEAN_pred -F PROVEAN_score -F PS1 -F PS2 -F PS3 -F PS4 -F PVS1 -F Polyphen2_HDIV_pred \
-F Polyphen2_HDIV_rankscore -F Polyphen2_HDIV_score -F Polyphen2_HVAR_pred -F Polyphen2_HVAR_rankscore -F Polyphen2_HVAR_score -F QD \
-F RAW_MQ -F ReadPosRankSum -F SIFT_converted_rankscore -F SIFT_pred -F SIFT_score -F SNPEFF_AMINO_ACID_CHANGE -F SNPEFF_CODON_CHANGE \
-F SNPEFF_EFFECT -F SNPEFF_EXON_ID -F SNPEFF_FUNCTIONAL_CLASS -F SNPEFF_GENE_BIOTYPE -F SNPEFF_GENE_NAME -F SNPEFF_IMPACT \
-F SNPEFF_TRANSCRIPT_ID -F SOR -F SiPhy_29way_logOdds -F SiPhy_29way_logOdds_rankscore -F VEST3_rankscore -F VEST3_score -F VQSLOD \
-F avsnp150 -F controls_AF_popmax -F culprit -F cytoBand -F fathmm-MKL_coding_pred -F fathmm-MKL_coding_rankscore -F fathmm-MKL_coding_score \
-F integrated_confidence_value -F integrated_fitCons_score -F integrated_fitCons_score_rankscore -F non_cancer_AF_popmax -F non_neuro_AF_popmax \
-F non_topmed_AF_popmax -F phastCons100way_vertebrate -F phastCons100way_vertebrate_rankscore -F phastCons20way_mammalian \
-F phastCons20way_mammalian_rankscore -F phyloP100way_vertebrate -F phyloP100way_vertebrate_rankscore -F phyloP20way_mammalian \
-F phyloP20way_mammalian_rankscore \
-GF GT \
-GF AD \
-GF DP \
-GF GQ \
-GF PL \
--output ${infile}.Variants-to-Table.txt

# End of script
```


**RefSeq gene annotation with ANNOVAR from a VCF file**

Credits: ANNOVAR, https://annovar.openbioinformatics.org/en/latest/user-guide/gene/

The output generates two files starting from "infile.vcf": 'infile.ANNOVAR.gene_annotation.variant_function' and 'infile.ANNOVAR.gene_annotation.exonic_variant_function'

Output "infile.ANNOVAR.gene_annotation.variant_function". The first column in the file output tells whether the variant hit exons or hit intergenic regions, or hit introns, or hit a non-coding RNA genes. If the variant is exonic/intronic/ncRNA, the second column gives the gene name (if multiple genes are hit, comma will be added between gene names); if not, the second column will give the two neighboring genes and the distance to these neighboring genes.

Output "infile.ANNOVAR.gene_annotation.exonic_variant_function". The second output file contains the amino acid changes as a result of the exonic variant.
The exact format of the output below may change slightly between different versions of ANNOVAR.

```Bash
#!/bin/bash

module load annovar/18.04.16

infile="${workdir}/infile.vcf"
outfile="${workdir}/infile.ANNOVAR.out.avinput"

convert2annovar.pl -format vcf4 ${infile} -outfile ${outfile} -allsample -withfreq -include

infile="${workdir}/infile.ANNOVAR.out.avinput"
outfile="${workdir}/infile.ANNOVAR.gene_annotation"
humandb="..ANNOVAR/humandb/"

annotate_variation.pl -out ${outfile} -build hg19 ${infile} ${humandb}

# Outputs two files: 'infile.ANNOVAR.gene_annotation.variant_function' and 'infile.ANNOVAR.gene_annotation.exonic_variant_function'

# End of script
```


**Split multiallelic variants (MNPs) into biallelic variants**

Source: http://samtools.github.io/bcftools/bcftools.html#norm

SNPs and indels will be merged separately into two records.

```Bash
#!/bin/bash

module load bcftools

bcftools norm -m -both infile.vcf.gz -Oz -o outfile.vcf.gz

# End of script
```


**Annotate the ID field of each variant in a VCF file using dbSNP database**

Source: http://samtools.github.io/bcftools/bcftools.html#annotate

```Bash
#!/bin/bash

module load htslib
module load bcftools

# dbSNP database index tbi file should be place within the same folder
dbsnp="../dbSNP/GCF_000001405.25.gz"

indir=".."
infile="${indir}/infile.vcf"
outfile="${infile}.dbsnp.vcf"

# First, compress and tabix the VCF file
bgzip -c ${infile} > ${infile}.gz
tabix -p vcf ${infile}.gz

# Then annotate the ID field
bcftools annotate -c ID -a ${dbsnp} -o ${outfile} ${infile}.gz

# End of script
```


**Annotate the ID field of each variant in a VCF file using dbSNP database**

Source: http://samtools.github.io/bcftools/bcftools.html#annotate

```Bash
#!/bin/bash

# Get the union of variants from two files
awk 'FNR==NR{a[$1];next}$1 on a' ${infile1} ${infile2} > ${union}

# End of script
```


**Parse a VCF file and count the number of HomRef, HomAlt, and Het genotypes from dose data**

In each line of a VCF file from Michigan Imputation Server (FORMAT field use to be GT:DS:GP), split the fields from the 10th column by "." (GT:DS:etc.). Then, split the "DS" but without using any delimiter... in other words, do not split. And then take the second field of the array, and assign it to DS variable. Then filter by DS value. From the sums "het+homalt" and "hethomref", we keep the minumum to filter out the variant afterwards.

```Bash
#!/bin/bash

zgrep -v "^#" ${indir}/${infile} \
| awk 'BEGIN{OFS="\t"} {
 homref=0; het=0; homalt=0; for(i=10;i<=NF;i++) {
  split($i,a,":"); split(a[2],DS);
  if(DS[1]<=0.5) {homref++}
  else if(DS[1]>1.5) {homalt++}
  else {het++}};
 if(het+homalt<het+homref) print $1,$2,$3,$4,$5,homref,het,homalt,het+homalt
  else if(het+homalt>=het+homref) print $1,$2,$3,$4,$5,homref,het,homalt,het+homref}' > ${outdir}/${outfile}
  
# End of script
```


**Transposing the row-vector of individuals into a column-vector**

```Bash
#!/bin/bash

awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str"\t"a[i,j];
        }
        print j"->"str" MLDOSE" # Prepare a matrix with columns 1, 2, and 3 following probABEL specifications
    }
}' $output/individuals > $output/t_individuals

# End of script
```


**Change the sample name in a BAM header using SAMtools**

```Bash
#!/bin/bash

# Define inbam and outbam files
samtools view -H ${inbam} | sed "s/SM:[^\t]*/SM:new-sample-name/g" | samtools reheader - ${bamfile} > ${outbam}

# End of script
```


**Sort a file using a certain numerical column while keeping the header at the top of the output**

```Bash
#!/bin/bash

# The header is the first line at the top of the file
# Column3 in the input file is numeric (i.e. physical position) and will be used as the index to order the dataset
awk 'NR == 1; NR > 1 {print $0 | "sort -k 3"}' ${infile} > ${outfile}

# End of script
```


**Convert a multisample columnar FASTA file into a multisample single-line FASTA file**

```Bash
#!/bin/bash

infile="columnar.fasta"
outfile="sinleline.fasta"

awk '{ 
if (NR==1) { print $0 }
else
    {
    if ($0 ~ /^>/)
        { print "\n"$0 }
    else
        { printf $0 }
    }
}' ${infile} > ${outfile}

# End of script
```


**Extract a region of interest from a BAM file with SAMtools and index**

```Bash
#!/bin/bash

# Use at HPC:
module load xz/5.2.5/gcc htslib/1.12/gcc samtools/1.12/gcc

#Define region
region="chr2:163140072-163160072"
indir="path-to-input-dir"
infile="input-file"
outfile="${infile}.${region}"
#Extract Region of Interest
samtools view -b -h ${indir}/${infile}.bam ${region} > ${indir}/${outfile}.bam
#Index
samtools index ${indir}/${outfile}.bam

# End of script
```


**Compress files without directory structure**

```Bash
#!/bin/bash

zip -r -j *.zip

# End of script
```


**Shorten the current directory path on terminal**

```Bash
#!/bin/bash

#See (1): https://www.gnu.org/software/bash/manual/html_node/Controlling-the-Prompt.html
#See (2): https://unix.stackexchange.com/questions/381113/how-do-i-shorten-the-current-directory-path-shown-on-terminal
#Open the shell, write this, and press ENTER:

PS1='\u:\W\$ '

#Done!
# End of script
```


**Check if directory/file exist/does not exist**

```Bash
#!/bin/bash

##### DIRECTORIES

#Check whether a directory exists or not.
indir="your-directory"

#Check if a directory exists
if [ -d "${indir}" ]; then
    echo "${indir} directory exists."
fi

#Check if a directory does not exist
if [ ! -d "${indir}" ]; then
    echo "${indir} directory does not exist."
fi

##### FILES

#Check whether file exists or not.
infile="your-file"

#Check if a file exists
if [ -f "${infile}" ]; then
    echo "${infile} exists."
fi

#Check is a file does not exist
if [ ! -f "${infile}"]; then
    echo "${infile} does not exist."
fi

#Check if a file is not empty
if [ -s "${infile}" ]; then
    echo "${infile} is not empty."
fi

#Together
if [ -f "${infile}" ]; then
    echo "${infile} exists."
else 
    echo "${infile} does not exist."
fi

#Check if multiple files exist
#Use '-a' or '&&' together with [[

#With '-a'
if [ -f "${infile1}" -a -f "${infile2}" ]; then
    echo "${infile1} and ${infile2} exist."
fi

#With '&&'
if [[ -f "${infile1}" && -f "${infile2}" ]]; then
    echo "${infile1} and ${infile2} exist."
fi

# End of script
```


**Determine strand orientation and allele switch in a VCF**

```Bash
#!/bin/bash

Source: http://samtools.github.io/bcftools/howtos/plugin.fixref.html

Use BCFtools plugins.

reference="path-to-reference-fasta"
infile="infile.vcf.gz"

#Stats to assess the number of REF allele mismatches and the number of non-biallelic sites
bcftools +fixref ${infile} -- -f ${reference}

#Check the reference allele mismatches
bcftools norm --check-ref e -f ${reference} ${infile} -Ou -o /dev/null

# End of script
```


**Get the file with the smallest number of lines

```Bash
#!/bin/bash

wc -l file1 > count1
wc -l file2 > count2
wc -l file3 > count3

if [ ${count1} -lt ${count2} ] && [ ${count1} -lt ${count3} ]
then
    cod=1
    count=${count1}
    file="${infile1}"
    echo "File number: ${cod}"
    echo "Number of lines: "${count}
elif [ ${count2} -lt ${count1} ] && [ ${count2} -lt ${count3} ]
then
    cod=2
    count=${count2}
    file="${infile2}"
    echo "File number: ${cod}"
    echo "Number of lines: "${count}
else
    cod=3
    count=${count3}
    file="${infile3}"
    echo "File number: ${cod}"
    echo "Number of lines: "${count}
fi

# End of script
```


**Get the shared positions by three lists**

```Bash
#!/bin/bash

#cat list1
#1
#2
#3
#4

#cat list2
#1
#2
#3
#4
#5

#cat list3
#2
#3
#4

# Cat list of positions
# Then sort them
# Count the number of duplicated lines
# Trim any leading and trailing white spaces
# Grab only the lines starting by 3 (3 lines duplicated)
# Keep only the selected lines and save into 'list123' file

cat list1 list2 list3 | sort -k1,1n | uniq -c | awk '{$1=$1;print}' | grep "^3" | awk '{ print $2 }' > list123

#cat list1 list2 list3 | sort -k1,1n | uniq -c
#      2 1
#      3 2
#      3 3
#      2 4
#      2 5
#      1 6
#      1 7
#      1 8

#cat list1 list2 list3 | sort -k1,1n | uniq -c | awk '{$1=$1;print}'
#2 1
#3 2
#3 3
#2 4
#2 5
#1 6
#1 7
#1 8

#cat list1233
#2
#3

# End of script
```


**Reredirect .out and .err outputs to a log file**

```Bash
#!/bin/bash

log="your-log-file"

./your-script.sh 2>&1 | tee -a ${log}
'your-line-commands-here' 2>&1 | tee -a ${log}

# End of script
```

**Filter a VCF file using certain values in the FILTER tag**

```Bash
#!/bin/bash

bcftools view -i "%FILTER='PASS' | %FILTER='.'" infile.vcf.gz -Oz -o outfile.vcf.gz 

# End of script
```


**Check that a VCF is sorted (via indexing)**

```Bash
#!/bin/bash

bgzip -c infile.vcf > infile.vcf.gz
bcftools index -t -f infile.vcf.gz

#If the latter command fails, the file is not sorted.

# End of script
```


**String Functions in GNU**

+ https://www.gnu.org/software/gawk/manual/html_node/String-Functions.html


**Sort numerically a list of files within a directory and save the list into a new file**

Source: https://stackoverflow.com/questions/13360925/sort-files-numerically-in-bash

```Bash
#!/bin/bash

#Search for and sort '*.chr22:a-b.Fst' files, where 'a' and 'b' are numbers (e.g., typically, a chromosome region)
#Use '-v1' parameter:
# -v: natural sort of (version) numbers within text
# -1: list one file per line

#List files numerically and save the list into ${list}
ls -v1 "*.Fst" > ${list}

# End of script
```


**Replace './.' or unknown genotypes with '0/0' or homozygous-reference genotypes**

```Bash
#!/bin/bash


sed 's/\.\/\./0\/0/g' infile.vcf > infile.genotypes-replaced.vcf

# End of script
```
