# Use a base image with Miniconda pre-installed
FROM continuumio/miniconda3:latest

# Set working directory
WORKDIR /app

# Copy environment.yml and requirements.txt
COPY environment.yml requirements.txt ./

# Create conda environment
RUN conda env create -f environment.yml

# Activate conda environment
SHELL ["conda", "run", "-n", "barcode_validator", "/bin/bash", "-c"]

# Copy program folders
COPY barcode_validator ./barcode_validator

# Add barcode_validator to PYTHONPATH
ENV PYTHONPATH="${PYTHONPATH}:/app/barcode_validator"

# Create a directory for the config
RUN mkdir /config

# Set the entry point with a default config location that can be overridden
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "barcode_validator", "python", "barcode_validator/daemon.py", "-c"]
CMD ["/config/config.yml"]