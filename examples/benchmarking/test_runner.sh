#!/bin/bash

# Base input file path
INPUT="data/BGE00146_MGE-BGE_r_1_s_100.fas"

# Function to format core count for config filename
format_cores() {
    echo "${1}_test_config.yml"
}

# Launch jobs for different core counts
for CORES in 1 2 4 8 16 32; do
    CONFIG_FILE="config/$(format_cores $CORES)"
    OUTPUT_TSV="${INPUT}_${CORES}.tsv"
    OUTPUT_LOG="${INPUT}_${CORES}.log"
    
    # Use parentheses to group commands and handle redirections properly
    (
        nohup bash -c "/usr/bin/time -v python barcode_validator -c $CONFIG_FILE -f $INPUT > $OUTPUT_TSV 2>> $OUTPUT_LOG"
    ) > /dev/null 2>&1 &
    
    echo "Launched job with $CORES cores"
    sleep 1  # Small delay to prevent potential race conditions
done

# Optional: wait for all background jobs to complete
# wait

echo "All jobs launched"
