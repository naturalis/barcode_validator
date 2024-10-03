#!/bin/bash

# Save me in /usr/local/bin/run-barcode-validator.sh so that I can be referenced from barcode-validator.service
# Check all the paths in config.yml and map any relevant paths into the container. Also ensure that the num_threads
# and BLAST RAM variable in the config.yml are compatible with the settings here:
docker run \
  --name barcode_validator \
  --restart unless-stopped \
  --cpus 56 \
  --memory 1300G \
  --memory-swap 1300G \
  -v /path/to/your/config.yml:/config/config.yml \
  -v /home/rutger.vos/data/ncbi:/home/rutger.vos/data/ncbi \
  barcode_validator