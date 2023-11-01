#!/bin/bash

# This is just a simple wrapper script around boldigger_cline
# that sets up certain defaults.

# The credentials should be in an .env file analogous to .env.example
source .env

# Marker can be coi, its, rbcl, case sensitive
MARKER=coi

# Input FASTA file is first argument
INFILE=$1

# Output folder is second
OUTFOLDER=$2

# Batch size 
BATCH=1

# Modify as per: 
# https://github.com/DominikBuchner/BOLDigger/tree/master#use-the-bold-identification-engine-for-coi-its-and-rbcl--matk
if [[ $MARKER == 'coi' ]]; then
	BATCH=50
elif [[ $MARKER == 'its' ]]; then
	BATCH=10
elif [[ $MARKER == 'rbcl' ]]; then
	BATCH=5
fi

# Do the incantation
boldigger-cline ie_${MARKER} $USER $PASS $INFILE $OUTFOLDER $BATCH
