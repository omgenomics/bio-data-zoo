URL_FASTQ = "https://42basepairs.com/download/s3/1000genomes/phase3/data/NA12878/sequence_read/ERR001268_1.filt.fastq.gz"
URL_BAM = "https://42basepairs.com/download/s3/1000genomes/phase3/data/NA12878/alignment/NA12878.chrom11.ILLUMINA.bwa.CEU.low_coverage.20121211.bam"
URL_BED = "https://42basepairs.com/download/s3/human-pangenomics/pangenomes/freeze/freeze1/minigraph/hprc-v1.0-minigraph-chm13.bb.bed.gz"
URL_VCF = "https://42basepairs.com/download/s3/1000genomes/technical/working/20140708_previous_phase3/chrXY_v1a/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz"
URL_VCF_MULTISAMPLE = "https://42basepairs.com/download/s3/1000genomes/technical/working/20140708_previous_phase3/chrXY_v1a/ALL.chrY.phase3_integrated_v1a.20130502.genotypes.vcf.gz"

# Non-file targets
.PHONY: all clean generate download

# Targets
clean:
	@rm -r ./data

generate:
	@find ./data -type f ! -name 'basic*' ! -name README.md -exec rm {} \;
	@./src/generate_fasta.sh && \
		./src/generate_fastq.sh && \
		./src/generate_bam.sh && \
		./src/generate_vcf.sh && \
		./src/generate_bed.sh

init:
	@mkdir -p ./data/{fasta,fastq,bam,vcf,bed}/{good,bad}
	@echo "Generating FASTA files..." && \
		echo ">sequence1\nAATTCTCATTACTGTATCACAGCAAGTTGTATTTACAACAAAAATCCAAA\n>sequence2\nGCCTACCAGAAAACGTTGTATTTTGGCAAAGTTCAAAAAGTCAGTCCAGA\n>sequence3\nGTATAATTCACAGAGTTTCATGTGGTTGTTGTTGACTCTACATATTGTCT" > "./data/fasta/good/basic_dna.fa" && \
		echo ">sequence1\nATCTACGATCGAGCTACT\n>sequence2\nATC----ATCGACCCACT" > "./data/fasta/good/basic_aligned.fa" && \
		echo ">sequence1\nADHWNARNNAKFWVYSHGPLWGIMHSHFPAGLAQGKNLHEIIPSMKQCIRPEWVDYCHMF\n>sequence2\nISTTGEGMSHFVQNWVPLVWGFAVHYAQVTLFRDTRNGGYEVSVEWLGLYVSQLDASWNI\n>sequence3\nVKMHIEVVRPIWEHSQNIHFAQLTDNPAAKACDGFAPVTMKKTCGTDTIHCYHTYHACWR" > "./data/fasta/good/basic_protein.fa"

	@echo "Downloading FASTQ files..." && \
		curl -s "$(URL_FASTQ)" | seqtk seq -l 0 | head -n 12 > "./data/fastq/good/basic_R1.fastq" && \
		curl -s "$(subst _1,_2,$(URL_FASTQ))" | seqtk seq -l 0 | head -n 12 > "./data/fastq/good/basic_R2.fastq"

	@echo "Downloading BAM files..." && \
		samtools view -b -h "$(URL_BAM)" 11:82365011-82366010 > "./data/bam/good/basic.bam" && \
		samtools view -b -h "$(URL_BAM)" 11:128989445-128990444 11:82365011-82366010 > "./data/bam/good/basic_unsorted.bam" && \
		samtools view -h "./data/bam/good/basic.bam" > "./data/bam/good/basic.sam"

	@echo "Downloading VCF files..." && \
		bcftools view --no-version "$(URL_VCF)" | head -n 300 > "./data/vcf/good/basic.vcf" && \
			bcftools view --no-version "$(URL_VCF_MULTISAMPLE)" | head -n 150 > "./data/vcf/good/basic_multisample.vcf" && \
			bcftools view --no-version "./data/vcf/good/basic.vcf" -o "./data/vcf/good/basic.bcf" && \
			bcftools view --no-version "./data/vcf/good/basic_multisample.vcf" -o "./data/vcf/good/basic_multisample.bcf"

	@echo "Downloading BED file..." && \
		curl -s "$(URL_BED)" | zcat | head -n 10 > "./data/bed/good/basic.bed"
