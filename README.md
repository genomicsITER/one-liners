<a name="inicio"></a>

<!-- ------------------ HEADER ------------------ -->
<!-- Developed and maintained by Genomics Division
<!-- of the Institute of Technology an Renewable Energy (ITER)
<!-- Tenerife, Canary Islands, SPAIN
<!-- See the "Contact us" section to collaborate with us to growth
<!-- this repository. ;=)

<!-- ------------------ SECTION ------------------ -->
<p align="left">
  <a href="https://github.com/genomicsITER/one-liners" title="Instituto Tecnológico y de Energ&iacute;as Renovables (ITER) / Institute of Technology and Renewable Energy (ITER)">
    <img src="https://github.com/genomicsITER/one-liners/blob/master/images/ITER_logo.png" width="auto" /> 
      </a>
</p>

# One-liners, frequently used commands in Bioinformatics

**Markdown References**

+ https://guides.github.com/features/mastering-markdown/

+ https://help.github.com/en/github/writing-on-github/basic-writing-and-formatting-syntax


**Table of Contents**

<details>
<summary>Click to display the content of this site:</summary>
<ul>
<li><a href="#code1">Select region from 1KGP chr.vcf.gz, extract, compress, and tabix</li></a>
<li><a href="#code2">Remove "chr" string in variants in a VCF</li></a>
<li><a href="#code3">Add "chr" string in variants in VCF</li></a>
<li><a href="#code4">Sort karyotipically a VCF (version 1)</li></a>
<li><a href="#code5">Sort karyotypically a VCF (version 2). Use '-V', natural sort of (version) numbers within text</li></a>
<li><a href="#code6">Sort karyotypically a VCF (version 3): using vcf-sort (from vcftools)</li></a>
<li><a href="#code7">Sort karyotypically a VCF (version 4): using PICARD</li></a>
<li><a href="#code8">Replace spaces with a single tab</li></a>
<li><a href="#code9">Compute BAM coverage with BEDtools</li></a>
<li><a href="#code10">Find duplicated lines in a VCF matching the whole line</li></a>
<li><a href="#code11">Find duplicated lines in a VCF matching files 1, 2, and 5</li></a>
<li><a href="#code12">Find a line by a field on it, delete it, and save the result</li></a>
<li><a href="#code13">Count the number of mapped reads for each read in PE reads</li></a>
<li><a href="#code14">Replace white spaces with tabs</li></a>
<li><a href="#code15">Add rs from INFO field (avsnp150) to ID column in a VCF</li></a>
<li><a href="#code16">Remove duplicated lines in a file while keeping the original order</li></a>
<li><a href="#code17">Number of files per extension type in the current directory</li></a>
<li><a href="#code18">Parse file with AWK, sum column values in each line, and shows the result</li></a>
<li><a href="#code19">Count genotypes in a VCF imputed at Michigan Imputation Server</li></a>
<li><a href="#code20">Replace spaces with tab</li></a>
<li><a href="#code21">Prepare a header</li></a>
<li><a href="#code22">Extract variants that have a genotype count equal or higher than a genotype count threshold</li></a>
<li><a href="#code23">Count the lines of each file in a dir and get the sum of all lines</li></a>
<li><a href="#code24">Count the total number of lines in the selected files of a dir</li></a>
<li><a href="#code25">Grab the body of a file excluding the header</li></a>
<li><a href="#code26">Count number of lines in a file</li></a>
<li><a href="#code27">Count number of columns in a file</li></a>
<li><a href="#code28">Replace many spaces with single tabs, and specially the leading spaces of a PLINK '*.frq' file</li></a>
<li><a href="#code29">Recode (0/0, 0/1, 1/0, 1/1) genotypes into (0,1,2) genotypes</li></a>
<li><a href="#code30">Convert a VCF into a table of variants combining PERL and GATK4 VariantsToTable</li></a>
<li><a href="#code31">RefSeq gene annotation with ANNOVAR from a VCF file</li></a>
<li><a href="#code32">Split multiallelic variants (MNPs) into biallelic variants</li></a>
<li><a href="#code33">Annotate the ID field of each variant in a VCF file using dbSNP database</li></a>
<li><a href="#code34">Annotate the ID field of each variant in a VCF file using dbSNP database</li></a>
<li><a href="#code35">Parse a VCF file and count the number of HomRef, HomAlt, and Het genotypes from dose data</li></a>
<li><a href="#code36">Transposing the row-vector of individuals into a column-vector</li></a>
<li><a href="#code37">Change the sample name in a BAM header using SAMtools</li></a>
<li><a href="#code38">Sort a file using a certain numerical column while keeping the header at the top of the output</li></a>
<li><a href="#code39">Convert a multisample columnar FASTA file into a multisample single-line FASTA file</li></a>
<li><a href="#code40">Extract a region of interest from a BAM file with SAMtools and index</li></a>
<li><a href="#code41">Compress files without directory structure</li></a>
<li><a href="#code42">Shorten the current directory path on terminal</li></a>
<li><a href="#code43">Check if directory/file exist/does not exist</li></a>
<li><a href="#code44">Determine strand orientation and allele switch in a VCF</li></a>
<li><a href="#code45">Get the file with the smallest number of lines</li></a>
<li><a href="#code46">Get the shared positions by three lists</li></a>
<li><a href="#code47">Reredirect .out and .err outputs to a log file</li></a>
<li><a href="#code48">Filter a VCF file using certain values in the FILTER tag</li></a>
<li><a href="#code49">Check that a VCF is sorted (via indexing)</li></a>
<li><a href="#code50">String Functions in GNU</li></a>
<li><a href="#code51">Sort numerically a list of files within a directory and save the list into a new file</li></a>
<li><a href="#code52">Replace './.' or unknown genotypes with '0/0' or homozygous-reference genotypes</li></a>
<li><a href="#code53">Replace metadata ids in a FASTA file with new ids</li></a>
<li><a href="#code54">Avoid the '--' group-separator in a grep output providing multiple matches</li></a>
<li><a href="#code55">Iterate over a list of elements within a tuple</li></a>
<li><a href="#code56">Annotate fields (AF, AC, MAF, etc.) in a VCF</li></a>
<li><a href="#code57">Concatenate VCF files into a single-all-chrs VCF file</li></a>
<li><a href="#code58">Monitor processes in your HPC (squeue, sacct, etc.)</li></a>
<li><a href="#code59">List file recursively by last modification time</li></a>
<li><a href="#code60">Remove extension from a file name</li></a>
<li><a href="#code61">Create a template <i>'Empy Document'</i> in order to show this option in the contextual menu in Ubuntu</li></a>
<li><a href="#code62">Set an infinite history length of bashrc in Ubuntu using HISTSIZE and HISTFILESIZE in bash</li></a>
<li><a href="#code63">Set a string to lowercase or uppercase in AWK</li></a>
<li><a href="#code64">Mount a network NTFS volume into a local folder from the shell</li></a>
<li><a href="#code65">Configure a workstation folder mounting on a remote NAS volume from the shell</li></a>
<li><a href="#code66">Grab the header of a file excluding the last line</li></a>
<li><a href="#code67">Break MNPs and exclude spanning deletions (SD) from a VCF/gVCF file using BCFtools</li></a>
<li><a href="#code68">Split a file with a list in a number of sublists using BASH</li></a>
<li><a href="#code69">Add new line ('\n') at the end of a file using BASH</li></a>
<li><a href="#code70">Keep only SNPs from a VCF using BCFtools</li></a>
<li><a href="#code71">Grab the list of variants from a VCF using BCFtools</li></a>
<li><a href="#code72">Merge a list of VCF files using BCFtools</li></a>
<li><a href="#code73">List the content of a file, line by line, showing line numbers with BASH</li></a>
<li><a href="#code74">Add a number of leading zeroes ('0') to folder names within a directory with BASH</li></a>
</details>

<hr>

<a name="code1"></a>

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


<hr>

<a name="code2"></a>

**Remove "chr" string in variants in a VCF**
```Bash
#!/bin/bash

awk '{gsub(/^chr/,""); print}' infile.vcf > infile.no_chr.vcf

# End of script
```


<hr>

<a name="code3"></a>

**Add "chr" string in variants in VCF**
```Bash
#!/bin/bash

awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' infile.no_chr.vcf > infile.vcf

# End of script
```


<hr>

<a name="code4"></a>

**Sort karyotipically a VCF (version 1)**
```Bash
#!/bin/bash

##Notice that header must be removed at some point
grep '^#' in.vcf > out.vcf && grep -v '^#' in.vcf | sort -V -k1,1 -k2,2n >> out.vcf

# End of script
```


<hr>

<a name="code5"></a>

**Sort karyotypically a VCF (version 2). Use '-V', natural sort of (version) numbers within text**
```Bash
#!/bin/bash

sort -V -k1,1 -k2,2n infile.vcf > outfile.vcf

# End of script
```


<hr>

<a name="code6"></a>

**Sort karyotypically a VCF (version 3): using vcf-sort (from vcftools)**
```Bash
#!/bin/bash

cat infile.vcf | vcf-sort --chromosomal-order > infile.sorted.vcf

# End of script
```


<hr>

<a name="code7"></a>

**Sort karyotypically a VCF (version 4): using PICARD**
```Bash
#!/bin/bash

java -jar picard.jar SortVcf I=unsorted.infile.vcf O=sorted.infile.vcf

# End of script
```


<hr>

<a name="code8"></a>

**Replace spaces with a single tab**
```Bash
#!/bin/bash

sed 's/ \+/\t/g' infile > outfile

# End of script
```


<hr>

<a name="code9"></a>

**Compute BAM coverage with BEDtools**
```Bash
#!/bin/bash

bedtools genomecov -ibam infile.bam -bg > coverage.txt

# End of script
```


<hr>

<a name="code10"></a>

**Find duplicated lines in a VCF matching the whole line**
```Bash
#!/bin/bash

awk ' !uniq[$0]++ ' infile.vcf

# End of script
```


<hr>

<a name="code11"></a>

**Find duplicated lines in a VCF matching files 1, 2, and 5**
```Bash
#!/bin/bash

awk ' !uniq[$1 FS $2 FS $5]++ ' infile.vcf

# End of script
```


<hr>

<a name="code12"></a>

**Find a line by a field on it, delete it, and save the result**
```Bash
#!/bin/bash

string=abc
grep -v $string infile > outfile

# End of script
```


<hr>

<a name="code13"></a>

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


<hr>

<a name="code14"></a>

**Replace white spaces with tabs**
```Bash
#!/bin/bash

awk -v OFS="\t" '$1=$1' infile > outfile

# End of script
```


<hr>

<a name="code15"></a>

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


<hr>

<a name="code16"></a>

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


<hr>

<a name="code17"></a>

**Number of files per extension type in the current directory**
```Bash
#!/bin/bash

nfiletypes () { find . -maxdepth 1 -type f | sed 's/.*\.//' | sort | uniq -c \
    | sed 's/^ *//g' | sed 's/ /\t/g'; }

# Note: run this code and then write "nfiletypes" at the prompt and will see the count of 
# files per extension at the current directory.

# End of script
```


<hr>

<a name="code18"></a>

**Parse file with AWK, sum column values in each line, and shows the result**
```Bash
#!/bin/bash

awk -F'[\t]' 'BEGIN{sum=0; OFS="\t"} { for (i=1;i<=NF;i++) a[i]+=$i } END{ for (i in a) print a[i] }' infile

# End of script
```


<hr>

<a name="code19"></a>

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


<hr>

<a name="code20"></a>

**Replace spaces with tab**
```Bash
#!/bin/bash

awk -v OFS="\t" '$1=$1' tmp1 > tmp2

# End of script
```


<hr>

<a name="code21"></a>

**Prepare a header**
```Bash
#!/bin/bash

echo -e "#CHR\tPOS\tID\tREF\tALT\tnumGT.Unknowns\tnumGT.HomRef\tnumGT.Het \
     \tnumGT.HomAlt\tnumGT.(Het+HomAlt)" > header<br>
cat header tmp2 > infile.variant-genotypes.counts

# End of script
```


<hr>

<a name="code22"></a>

**Extract variants that have a genotype count equal or higher than a genotype count threshold**
```Bash
#!/bin/bash

awk -F'[ ]' '{ if ($10 >= 5) print $3 }' infile.variant-genotypes.counts > variant-list

# End of script
```


<hr>

<a name="code23"></a>

**Count the lines of each file in a dir and get the sum of all lines**

Credits: http://stackoverflow.com/questions/13727917/ddg#13728131

```Bash
#!/bin/bash

find ${indir} -type f -name "*.selected-files" -print0 | wc -l --files0-from=-

# End of script
```


<hr>

<a name="code24"></a>

**Count the total number of lines in the selected files of a dir**

Credits: http://stackoverflow.com/questions/13727917/ddg#13728131

```Bash
#!/bin/bash

find ${indir} -type f -name "*.selected-files" -exec cat {} + | wc -l

# End of script
```


<hr>

<a name="code25"></a>

**Grab the body of a file excluding the header**
```Bash
#!/bin/bash

tail -n +2 ${infile} > file-without-header

# End of script
```


<hr>

<a name="code26"></a>

**Count number of lines in a file**
```Bash
#!/bin/bash

wc -l ${infile}

# End of script
```


<hr>

<a name="code27"></a>

**Count number of columns in a file**
```Bash
#!/bin/bash

head -n 1 ${infile} | awk '{print NF}'
#Or
awk '{print NF; exit}' ${infile}
# End of script
```


<hr>

<a name="code28"></a>

**Replace many spaces with single tabs, and specially the leading spaces of a PLINK '*.frq' file**
```Bash
#!/bin/bash

sed 's/ \+/\t/g' ${infile} | sed -e 's/^[\t]*//' >${infile}.no-trailing-spaces.tabs

# End of script
```


<hr>

<a name="code29"></a>

**Recode (0/0, 0/1, 1/0, 1/1) genotypes into (0,1,2) genotypes**
```Bash
#!/bin/bash

sed 's/0|0/0/g' test | sed 's/0|1/1/g' | sed 's/1|0/1/g' | sed 's/1|1/2/g

# End of script
```


<hr>

<a name="code30"></a>

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


<hr>

<a name="code31"></a>

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


<hr>

<a name="code32"></a>

**Split multiallelic variants (MNPs) into biallelic variants**

Source: http://samtools.github.io/bcftools/bcftools.html#norm

SNPs and indels will be merged separately into two records.

```Bash
#!/bin/bash

module load bcftools

bcftools norm -m -both infile.vcf.gz -Oz -o outfile.vcf.gz

# End of script
```


<hr>

<a name="code33"></a>

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


<hr>

<a name="code34"></a>

**Annotate the ID field of each variant in a VCF file using dbSNP database**

Source: http://samtools.github.io/bcftools/bcftools.html#annotate

```Bash
#!/bin/bash

# Get the union of variants from two files
awk 'FNR==NR{a[$1];next}$1 on a' ${infile1} ${infile2} > ${union}

# End of script
```


<hr>

<a name="code35"></a>

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


<hr>

<a name="code36"></a>

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


<hr>

<a name="code37"></a>

**Change the sample name in a BAM header using SAMtools**

```Bash
#!/bin/bash

# Define inbam and outbam files
samtools view -H ${inbam} | sed "s/SM:[^\t]*/SM:new-sample-name/g" | samtools reheader - ${bamfile} > ${outbam}

# End of script
```


<hr>

<a name="code38"></a>

**Sort a file using a certain numerical column while keeping the header at the top of the output**

```Bash
#!/bin/bash

# The header is the first line at the top of the file
# Column3 in the input file is numeric (i.e. physical position) and will be used as the index to order the dataset
awk 'NR == 1; NR > 1 {print $0 | "sort -k 3"}' ${infile} > ${outfile}

# End of script
```


<hr>

<a name="code39"></a>

**Convert a multisample columnar FASTA file into a multisample single-line FASTA file**

```Bash
#!/bin/bash

infile="columnar.fasta"
outfile="singleline.fasta"

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


<hr>

<a name="code40"></a>

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


<hr>

<a name="code41"></a>

**Compress files without directory structure**

```Bash
#!/bin/bash

zip -r -j *.zip

# End of script
```


<hr>

<a name="code42"></a>

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


<hr>

<a name="code43"></a>

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


<hr>

<a name="code44"></a>

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


<hr>

<a name="code45"></a>

**Get the file with the smallest number of lines**

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


<hr>

<a name="code46"></a>

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


<hr>

<a name="code47"></a>

**Reredirect .out and .err outputs to a log file**

```Bash
#!/bin/bash

log="your-log-file"

./your-script.sh 2>&1 | tee -a ${log}
'your-line-commands-here' 2>&1 | tee -a ${log}

# End of script
```


<hr>

<a name="code48"></a>

**Filter a VCF file using certain values in the FILTER tag**

```Bash
#!/bin/bash

bcftools view -i "%FILTER='PASS' | %FILTER='.'" infile.vcf.gz -Oz -o outfile.vcf.gz 

# End of script
```


<hr>

<a name="code49"></a>

**Check that a VCF is sorted (via indexing)**

```Bash
#!/bin/bash

bgzip -c infile.vcf > infile.vcf.gz
bcftools index -t -f infile.vcf.gz

#If the latter command fails, the file is not sorted.

# End of script
```


<hr>

<a name="code50"></a>

**String Functions in GNU**

+ https://www.gnu.org/software/gawk/manual/html_node/String-Functions.html



<hr>

<a name="code51"></a>

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


<hr>

<a name="code52"></a>

**Replace './.' or unknown genotypes with '0/0' or homozygous-reference genotypes**

```Bash
#!/bin/bash

sed 's/\.\/\./0\/0/g' infile.vcf > infile.with-replaced-genotypes.vcf

#Alternative using a BCFtools plugin (replace '.' with '0')
bcftools +setGT infile.vcf -Ov -o infile.with-replaced-genotypes.vcf -- -t . -n 0

# End of script
```


<hr>

<a name="code53"></a>

**Replace metadata ids in a FASTA file with new ids**

Source: https://unix.stackexchange.com/questions/652523/replacing-the-seq-ids-of-fasta-file-based-on-the-new-ids-from-a-list

```Bash
#!/bin/bash

#Files:
#'map' = a two-column file with old and new ids
#'old-ids.fasta' = a FASTA file with old-ids
#'new-ids.fasta' = a new FASTA file with new ids replacing old-ids

#'map' file is a two column file as follows:
#old-id1 [tab] new-id1
#old-id2 [tab] new-id2
#old-id3 [tab] new-id3
#... [tab] ...

#Example of 'map' file:
#HUNSC_ITER_150265155	hCoV-19/Spain/CN-HUNSC_ITER_150265155/2022
#HUNSC_ITER_900188949	hCoV-19/Spain/CN-HUNSC_ITER_900188949/2022
#HUNSC_ITER_180111581	hCoV-19/Spain/CN-HUNSC_ITER_180111581/2022
#HUNSC_ITER_402054223	hCoV-19/Spain/CN-HUNSC_ITER_402054223/2022
#...    ...

#Example of 'old-ids.fasta' file:
#>HUNSC_ITER_150265155
#TCTTGTAGATCTGTTCTCTAAACGAACTTTAAAAT...
#>HUNSC_ITER_900188949
#TTGTAGATCTGTTCTCTAAACGAACTTTAAAATCT...
#>HUNSC_ITER_180111581
#TTGTAGATCTGTTCTCTAAACGAACTTTAAAATCT...

awk -F'\t' '
    NR==FNR { map[">"$1] = ">"$2; next }
    $0 in map { $0 = map[$0] }
    { print }
' map old-ids.fasta > new-ids.fasta

#Example of output in 'new-ids.fasta' file:
#>hCoV-19/Spain/CN-HUNSC_ITER_150265155/2022
#TCTTGTAGATCTGTTCTCTAAACGAACTTTAAAAT...
#>hCoV-19/Spain/CN-HUNSC_ITER_900188949/2022
#TTGTAGATCTGTTCTCTAAACGAACTTTAAAATCT...
#>hCoV-19/Spain/CN-HUNSC_ITER_180111581/2022
#TTGTAGATCTGTTCTCTAAACGAACTTTAAAATCT...

# End of script
```


<hr>

<a name="code54"></a>

**Avoid the '--' group-separator in a grep output providing multiple matches**

```Bash
#!/bin/bash

# According to 'man grep':
#-A NUM, --after-context=NUM
#        Print NUM  lines  of  trailing  context  after  matching  lines.
#        Places   a  line  containing  a  group  separator  (--)  between
#        contiguous groups of matches.  With the  -o  or  --only-matching
#        option, this has no effect and a warning is given.
    
#Grep a multisample linear FASTA file witn '-A' and get rid of the annoying
#'--' lines when multiple matches are found
grep -A 1 --no-group-separator infile.fasta > outfile.fasta

#Or alternatively clean the output avoind the '--' lines
my-grep-command | sed '/^--$/d'

# End of script
```


<hr>

<a name="code55"></a>

**Iterate over a list of elements within a tuple**

Source: https://stackoverflow.com/questions/9713104/loop-over-tuples-in-bash

```Bash
#!/bin/bash

#Example:
#Instead of nesting two for-loops (one for population and one for number of individuals 
#within the selected population), define an array of tuples with two elements each.

for tuple in "CEU 100" "IBS 121" "GBR 106"
do
# convert the "tuple" into the param args $1 $2...
set -- ${tuple} 
pop=${1}
n=${2}
done

#Alternatively, define explicitely the separator between tuple values using IFS
for tuple in CEU,100 IBS,121 GBR,106
do 
IFS=',' read -a strarr <<< "${tuple}"
pop=${strarr[0]}
n=${strarr[1]}
done

# End of script
```


<hr>

<a name="code56"></a>

**Annotate fields (AF, AC, MAF, etc.) in a VCF**

Source: https://www.biostars.org/p/180894/

```Bash
#!/bin/bash

infile="in.vcf"
outfile="out.vcf"

bcftools +fill-tags ${infile}  -- -t AF,AC,MAF > ${outfile}

# End of script
```

<hr>

<a name="code57"></a>

**Concatenate VCF files into a single-all-chrs VCF file**

```Bash
#!/bin/bash

#Prepare a list of individuals VCF.gz files sorted numerically
ls -v1 *.vcf.gz > vcf.gz.list
cat vcf.gz.list

#Concatenate VCF files into a single-all-chrs VCF file
bcftools concat -f vcf.gz.list -Oz -o all-chrs-VCF.vcf.gz 

# End of script
```

<hr>

<a name="code58"></a>

**Monitor processes in your HPC (squeue, sacct, etc.)**

```Bash
#!/bin/bash

#Use the shell and connect to your HPC running SLURM
#Example1: monitor partitions and jobs
watch squeue

#Example2: monitor running jobs and subjobs
watch sacct

# End of script
```

<hr>

<a name="code59"></a>

**List file recursively by last modification time**

```Bash
#!/bin/bash

ls -Rt

# End of script
```

<hr>

<a name="code60"></a>

**Remove extension from a file name**

```Bash
#!/bin/bash

filename="basename.ext"

n=${filename%.*}

echo ${n}`

# End of script
```

<hr>

<a name="code61"></a>

**Create a template <i>'Empy Document'</i> in order to show this option in the contextual menu in Ubuntu**

```Bash
#!/bin/bash

#Open a shell and run this:
touch ~/Templates/Empty\ Document

#Then go to any folder and right-click to see that the contextual
#item 'New Document' > 'Empty Document' has been created.

# End of script
```

<hr>

<a name="code62"></a>

**Set an infinite history length of bashrc in Ubuntu using HISTSIZE and HISTFILESIZE in bash**

```Bash
#!/bin/bash

#Edit the '.bashrc' file in your /home/user directory
#Find the HISTSIZE and HISTFILESIZE sections and replace with these values:

# Numeric values less than zero result in every command being saved on the history list (there is no limit). 
# Old value: HISTSIZE=1000
HISTSIZE=-1

# Non-numeric values and numeric values less than zero inhibit truncation.
# Old value: HISTFILESIZE=2000
HISTFILESIZE=-1

# End of script
```
<hr>

<a name="code63"></a>

**Set a string to lowercase or uppercase in AWK**

```Bash
#!/bin/bash

#To lowercase field in position 1
awk -F'[\t]' '{ tolower($1) }' infile > outfile

#To uppercase field in position 2
awk -F'[\t]' '{ tolower($2) }' infile > outfile

# End of script
```
<hr>

<a name="code64"></a>

**Mount a network NTFS volume into a local folder from the shell**

```Bash
#!/bin/bash

sudo mount.cifs -v //192.168.xx.yy/<remote-folder> \
/mnt/cifs1 -o sec=ntlmv2,vers=1.0,domain=<your-domain>,username=<username>

# End of script
```
<hr>

<a name="code65"></a>

**Configure a workstation folder mounting on a remote NAS volume from the shell**

```Bash
#!/bin/bash

sudo mount -t cifs \
-o user=<username>,password=<password>,iocharset=utf8,noperm \
//192.168.xx.yy/<remote-folder> \
/mnt/<local-folder>

# End of script
```

<hr>

<a name="code66"></a>

**Grab the header of a file excluding the last line**

```Bash
#!/bin/bash

#Example: keep the header of a VCF file except the last line (i.e., the line with CHR, POS, ID,
#REF, ALT, and the rest of fields).

head -n -1 ${infile} > file-without-the-last-line

# End of script
```

<hr>

<a name="code67"></a>

**Break MNPs and exclude spanning deletions (SD) from a VCF/gVCF file using BCFtools**

```Bash
#!/bin/bash

#Define your files
infile="file.vcf.gz"
outfile="file.noMNP.noSD.vcf.gz"

# For each single VCF or gVCF, split MNPS and avoid spanning-deletions (marked as ALT=*)(SD), and index
bcftools norm -m -both ${infile} -Ou | \
bcftools view -e 'ALT="*"' -Oz -o ${outfile}
bcftools index -f -t ${outfile}

# End of script
```

<hr>

<a name="code68"></a>

**Split a file with a list in a number of sublists using BASH**

```Bash
#!/bin/bash

# Split a 'list' file (i.e., with sample codes in rows) into 10 new sublists
# 'sublist' is the name preffix for new files

split -d list sublist -n 10

# End of script
```

<hr>

<a name="code69"></a>

**Add new line ('\n') at the end of a file using BASH**

```Bash
#!/bin/bash

sed -i -z 's/$/\n/g' infile > outfile

# End of script
```

<hr>

<a name="code70"></a>

**Keep only SNPs from a VCF using BCFtools**

```Bash
#!/bin/bash

infile="in.vcf.gz"
outfile="out.onlySNPs.vcf.gz"

bcftools view --types snps -m 2 -M 2 ${infile} -Oz -o ${outfile}

# End of script
```

<hr>

<a name="code71"></a>

**Grab the list of variants from a VCF using BCFtools**

```Bash
#!/bin/bash

infile="in.vcf.gz"
outfile="variants.list"

bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' ${infile} > ${outfile}

# End of script
```

<hr>

<a name="code72"></a>

**Merge a list of VCF files using BCFtools**

```Bash
#!/bin/bash

list="list.of.vcf.files"
outfile="merged.vcf.gz"

#Merge
bcftools merge -m none -l ${list} -Oz -o ${outfile}
#Index
bcftools index -f -t ${outfile}

# End of script
```

<hr>

<a name="code73"></a>

**List the content of a file, line by line, showing line numbers with BASH**

```Bash
#!/bin/bash
infile="filename-here"

#Option1
awk '{print NR  "> " $s}' ${infile}

#Option2
#See 'nl --help'  for options
nl -b a ${infile}

#Option3
#See 'pr --help' for options
pr -t -n ${infile}

# End of script
```

<hr>

<a name="code74"></a>

**Add a number of leading zeroes ('0') to folder names within a directory with BASH**

```Bash
#!/bin/bash

#Rename folders in a directory to write leading zeroes to folder names
#BEFORE (list of folders)
#1_folderA
#2_folderB
#...
#11_folderK
#12_folderL
#...

#BEFORE (list of renamed folders)
#001_folderA
#002_folderB
#...
#011_folderK
#012_folderL
#...

#The number of leading zeroes is controled by '%03d' (for three zeroes)

indir="my-dir"
cd ${indir}

for d in *; do
  if [ -d "$d" ] && [[ "$d" =~ ^([0-9]+)_(.*)$ ]]; then
    num="${BASH_REMATCH[1]}"
    rest="${BASH_REMATCH[2]}"
    new="$(printf "%03d" "$num")_$rest"
    if [ "$d" != "$new" ]; then
      mv -v -- "$d" "$new"
    fi
  fi
done

# End of script
```

<hr>
