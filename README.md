# Bio Data Zoo

This repo contains example data in various genomics file formats. It is intended for bioinformatics tool developers to make testing software easier. It includes examples of valid file formats, edge cases, and invalid formats.

## Browse

**TODO**

## Formats included

|Format|Extensions|
|--|--|
|FASTA|.fa, .fa.gz|
|FASTQ|.fastq, .fastq.gz|
|BAM|.bam, .bam.bai, .bam.csi, .sam, .sam.gz, .sam.gz.csi, .sam.gz.tbi|
|VCF|.vcf, .vcf.gz, .vcf.gz.csi, .vcf.gz.tbi, .bcf, .bcf.csi|
|BED|.bed, .bed.gz, .bed.gz.csi, .bed.gz.tbi|
|CRAM|**TODO**|
|GFF|**TODO**|


## Data Source

|Path|Source|Preview file|Download file|
|--|--|--|--|
| `fastq/good/basic_R1.fastq` | `s3://1000genomes` | [Preview on 42basepairs](https://42basepairs.com/browse/s3/1000genomes/phase3/data/NA12878/sequence_read?file=ERR001268_1.filt.fastq.gz&preview=) | [Download file](https://42basepairs.com/download/s3/1000genomes/phase3/data/NA12878/sequence_read/ERR001268_1.filt.fastq.gz) |
| `bam/good/basic.bam` | `s3://1000genomes` | [Preview on 42basepairs](https://42basepairs.com/browse/s3/1000genomes/phase3/data/NA12878/alignment?file=NA12878.chrom11.ILLUMINA.bwa.CEU.low_coverage.20121211.bam&preview=) | [Download file](https://42basepairs.com/download/s3/1000genomes/phase3/data/NA12878/alignment/NA12878.chrom11.ILLUMINA.bwa.CEU.low_coverage.20121211.bam) |
| `vcf/good/basic_multisample.vcf` | `s3://human-pangenomics` | [Preview on 42basepairs](https://42basepairs.com/browse/s3/1000genomes/technical/working/20140708_previous_phase3/chrXY_v1a?file=ALL.chrY.phase3_integrated_v1a.20130502.genotypes.vcf.gz&preview=) | [Download file]("https://42basepairs.com/download/s3/1000genomes/technical/working/20140708_previous_phase3/chrXY_v1a/ALL.chrY.phase3_integrated_v1a.20130502.genotypes.vcf.gz") |
| `vcf/good/basic.vcf` | `s3://human-pangenomics` | [Preview on 42basepairs](https://42basepairs.com/browse/s3/1000genomes/technical/working/20140708_previous_phase3/chrXY_v1a?file=ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz&preview=) | [Download file](https://42basepairs.com/download/s3/1000genomes/technical/working/20140708_previous_phase3/chrXY_v1a/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz) |
| `bed/good/basic.bed` | `s3://human-pangenomics` | [Preview on 42basepairs](https://42basepairs.com/browse/s3/human-pangenomics/pangenomes/freeze/freeze1/minigraph?file=hprc-v1.0-minigraph-chm13.bb.bed.gz&preview=) | [Download file](https://42basepairs.com/download/s3/human-pangenomics/pangenomes/freeze/freeze1/minigraph/hprc-v1.0-minigraph-chm13.bb.bed.gz) |


## Contributing

**TODO**

## Dev (skip this section if only using the data)

### Prerequisites to develop on this repo:

* `samtools`
* `bedtools`
* `bcftools`
* `tabix`
* `bgzip`
* `seqtk`
* `gshuf` (`brew install coreutils` on Mac)
* `zcat`

### File names

Basic example files are called `basic*`, and files derived from it must not have a file name that start with that pattern.
