# RES-Scanner2
RES-Scanner 2 is a software package for genome-wide identification and annotation of RNA-editing sites for any species with matched RNA-Seq and DNA-Seq data. Compared with last version of perl, we add the function of identification for hyper editing sites, optimize some logic problems, and reduce the labor intensity to run. We also make this pipeline adapt to many different inputs as possible.
## Requirements
Python >= 3.3 and [Pysam](https://pysam.readthedocs.io/en/latest/api.html).  
[BWA](http://bio-bwa.sourceforge.net/), [SOAPnuke](https://github.com/BGI-flexlab/SOAPnuke), [JAVA](https://en.wikipedia.org/wiki/Java_(programming_language)), [PILON](https://github.com/broadinstitute/pilon/wiki) and [BLAT](https://genome.ucsc.edu/FAQ/FAQblat.html) are needed by some functions.  
## Examples
### 1. Paired-end DNA-seq and fr strand-specific RNA-seq data of one sample.
Input files:  
* dna_1.fq, dna_2.fq: paired-end DNA-seq reads
* rna_1.fq, rna_2.fq: fr strand-specific RNA-seq reads
* config: config file in json format, used in step1 and should be set manually.
* ref.fa: reference genome fasta file

Notes:
* Commands in step 2~3 could be run in parallel.
* Input files are in the testdata directory.
* Ires is a shell script in the bin directory of the package, which calls the same name of python script in the library directory. So you could use python directly.

##### Step1: preprocess input data to generate bam file.
```
ires preprocess -o preprocess/ --bwa <path/to/bwa> --soapnuke <path/to/soapnuke> --java <path/to/java> --pilon <path/to/pilon.jar> --DNA-fq1 dna_1.fq --DNA-fq2 dna_2.fq --RNA-fq1 rna_1.fq --RNA-fq2 rna_2.fq --config config --ss T ref.fa
```
##### Step2: transform bam result to single base format.
```
ires sam2base --trim 0,0 preprocess/ref.fa.fix_snps.fa preprocess/DNA_aln.best.bam
ires sam2base --trim 6,6 preprocess/ref.fa.fix_snps.fa preprocess/RNA_aln.negative.bam
ires sam2base --trim 6,6 preprocess/ref.fa.fix_snps.fa preprocess/RNA_aln.positive.bam
```
##### Step3: identify editing sites in genome wide.
```
ires scanner -o scanner/ --DNA-depth 10 --strand - --rate 2 --blat <path/to/blat> preprocess/ref.fa.fix_snps.fa preprocess/DNA_aln.best.bam.sb.gz preprocess/RNA_aln.negative.bam.sb.gz preprocess/RNA_aln.negative.bam
ires scanner -o scanner/ --DNA-depth 10 --strand + --rate 2 --blat <path/to/blat> preprocess/ref.fa.fix_snps.fa preprocess/DNA_aln.best.bam.sb.gz preprocess/RNA_aln.positive.bam.sb.gz preprocess/RNA_aln.positive.bam
```
##### Step4: identify hyper-editing sites using unmapped reads in step1.
```
ires hyper -o hyper/ preprocess/ref.fa.fix_snps.fa preprocess/aln/RNA_aln.bam preprocess/DNA_aln.best.bam.sb.gz preprocess/soapnuke/RNA/rna_1.fq_1.clean.fq.gz preprocess/soapnuke/RNA/rna_2.fq_2.clean.fq.gz
```
