#!/bin/bash

for CORES in 1 2 4 8 16 32; do
    (
        env BLASTDB=/home/rutger.vos/data/ncbi/nt/ BLASTDB_LMDB_MAP_SIZE=1024G time blastn \
            -task megablast \
            -evalue 1e-5 \
            -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore staxids" \
            -num_threads $CORES \
            -max_target_seqs 10 \
            -word_size 28 \
            -dust "20 64 1" \
            -soft_masking true \
            -db /home/rutger.vos/data/ncbi/nt/nt \
            -query /tmp/tmpo173d0or.fasta \
            -taxids 50557 \
            -out /tmp/tmpo173d0or.fasta.tsv 2>> time.log
    )
done
