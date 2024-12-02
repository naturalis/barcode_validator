#!/bin/bash

for MEM in 67108864 134217728 268435456 536870912 1073741824; do
    MEM_GB=$((MEM/1048576))
    echo "Testing with ${MEM_GB}GB limit" >> time_mem.log
    (
	export BLASTDB_LMDB_MAP_SIZE=${MEM_GB}G
        ulimit -v $MEM
        env BLASTDB=/home/rutger.vos/data/ncbi/nt/ time blastn \
            -task megablast \
            -evalue 1e-5 \
            -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore staxids" \
            -num_threads 32 \
            -max_target_seqs 10 \
            -word_size 28 \
            -dust "20 64 1" \
            -soft_masking true \
            -db /home/rutger.vos/data/ncbi/nt/nt \
            -query /tmp/tmpo173d0or.fasta \
            -taxids 50557 \
            -out "/tmp/tmpo173d0or.fasta_${MEM_GB}GB.tsv" 2>> time_mem.log
    )
    sleep 5
done
