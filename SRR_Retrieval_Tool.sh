#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

SRRLIST=$SCRIPT_DIR/../data/srr_list

cd $SCRIPT_DIR/../data

cat $SRRLIST | while read line; do


    if [[ $line == "SRR"* ]]; then

        # retrieve fastq files from SRA
        wget -O $line.fastq.gz "https://www.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=$line&clipped=1"

        gzip -d $line.fastq.gz
        # deinterleave fastq file
        paste - - - - - - - - < $line.fastq \
    | tee >(cut -f 1-4 | tr "\t" "\n" > $line"_r1".fastq) \
    | cut -f 5-8 | tr "\t" "\n" > $line"_r2".fastq

    fi

done
