# Bio Data Zoo

This repo contains example data in various genomics file formats. It is intended for bioinformatics tool developers to make testing software easier. It includes examples of valid file formats, edge cases, and invalid formats.

## Browse

**TODO**

## Formats included

|Format|Extensions|
|--|--|
|BAM|.bam, .bam.bai, .bam.csi, .sam, .sam.gz, .sam.gz.csi, .sam.gz.tbi|
|BED|.bed, .bed.gz, .bed.gz.csi, .bed.gz.tbi|
|FASTA|.fa, .fa.gz|
|FASTQ|.fastq, .fastq.gz|
|CRAM|**TODO**|
|VCF|**TODO**|
|GFF|**TODO**|


## Data Source

|Path|Source|Preview file|Download file|
|--|--|--|--|
| `fastq/good/basic_R1.fastq` | `s3://1000genomes` | [Preview on 42basepairs](https://42basepairs.com/browse/s3/1000genomes/phase3/data/NA12878/sequence_read?file=ERR001268_1.filt.fastq.gz&preview=) | [Download file](https://42basepairs.com/download/s3/1000genomes/phase3/data/NA12878/sequence_read/ERR001268_1.filt.fastq.gz) |
| `bam/good/basic.bam` | `s3://1000genomes` | [Preview on 42basepairs](https://42basepairs.com/browse/s3/1000genomes/phase3/data/NA12878/alignment?file=NA12878.chrom11.ILLUMINA.bwa.CEU.low_coverage.20121211.bam&preview=) | [Download file](https://42basepairs.com/download/s3/1000genomes/phase3/data/NA12878/alignment/NA12878.chrom11.ILLUMINA.bwa.CEU.low_coverage.20121211.bam) |
| `bed/good/basic.bed` | `s3://human-pangenomics` | [Preview on 42basepairs](https://42basepairs.com/browse/s3/human-pangenomics/pangenomes/freeze/freeze1/minigraph?file=hprc-v1.0-minigraph-chm13.bb.bed.gz&preview=) | [Download file](https://42basepairs.com/download/s3/human-pangenomics/pangenomes/freeze/freeze1/minigraph/hprc-v1.0-minigraph-chm13.bb.bed.gz) |


## Contributing

**TODO**

## Dev (skip this section if only using the data)

### Prerequisites to develop on this repo:

* `samtools`
* `bedtools`
* `seqtk`
* `bgzip`
* `tabix`
* `gshuf` (`brew install coreutils` on Mac)

### File names

Basic example files are called `basic*`, and files derived from it must not have a file name that start with that pattern.
