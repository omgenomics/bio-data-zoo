#!/bin/bash

./src/generate_fastq.sh && \
    ./src/generate_bam.sh && \
    ./src/generate_bed.sh
