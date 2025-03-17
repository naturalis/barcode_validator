| File list | Description |
| --- | --- |
XE-4013_barcodes_out.fasta |	1,431 barcode consensus sequences (hitpicked by fasta_compare) from 1,804 sample XE-4013 datasetc snakemake workflow, run with no contaminant sequences, 6 r and s MGE parameter combinations, using both 'fastp' and 'standard' MGE snakemake workflow |
XE-4013_barcodes_out.submitted.fasta |	1,357 'project' barcode consensus sequences from XE-4013_barcodes_out.fasta |
XE-4013_barcodes_out.unsubmitted.fasta |	74 'non-project' barcode consensus sequences from XE-4013_barcodes_out.fasta | 
chiro_demux_barcode_out.fasta |	106 barcode consensus sequences  (hitpicked by fasta_compare) from 279 sample Chironomid dataset, run with no contaminant sequences, 6 r and s MGE parameter combinations, using both 'fastp' and 'standard' MGE snakemake workflow |
mge_fastp_r13s100_nocontam.fasta |	Barcode consensus sequences from 570 sample benchmarking dataset, run with no contaminant sequences, r 1.3 and s 100 MGE parameters, using the 'fastp' (i.e. fastp-merge version of MGE snakemake workflow) |
mge_fastp_r13s50_nocontam.fasta |	Barcode consensus sequences from 570 sample benchmarking dataset, run with no contaminant sequences, r 1.3 and s 50 MGE parameters, using the 'fastp' (i.e. fastp-merge version of MGE snakemake workflow) |
mge_fastp_r15s100_nocontam.fasta |	Barcode consensus sequences from 570 sample benchmarking dataset, run with no contaminant sequences, r 1.5 and s 100 MGE parameters, using the 'fastp' (i.e. fastp-merge version of MGE snakemake workflow) |
mge_fastp_r15s50_nocontam.fasta |	Barcode consensus sequences from 570 sample benchmarking dataset, run with no contaminant sequences, r 1.5 and s 50 MGE parameters, using the 'fastp' (i.e. fastp-merge version of MGE snakemake workflow) |
mge_fastp_r1s100_nocontam.fasta |	Barcode consensus sequences from 570 sample benchmarking dataset, run with no contaminant sequences, r 1 and s 100 MGE parameters, using the 'fastp' (i.e. fastp-merge version of MGE snakemake workflow) |
mge_fastp_r1s50_nocontam.fasta |	Barcode consensus sequences from 570 sample benchmarking dataset, run with no contaminant sequences, r 1 and s 50 MGE parameters, using the 'fastp' (i.e. fastp-merge version of MGE snakemake workflow) |
mge_standard_r13s100_nocontam.fasta |	Barcode consensus sequences from 570 sample benchmarking dataset, run with no contaminant sequences, r 1.3 and s 100 MGE parameters, using the 'standard' (i.e. concat version of MGE snakemake workflow) |
mge_standard_r13s50_nocontam.fasta |	Barcode consensus sequences from 570 sample benchmarking dataset, run with no contaminant sequences, r 1.3 and s 50 MGE parameters, using the 'standard' (i.e. concat version of MGE snakemake workflow) |
mge_standard_r15s100_nocontam.fasta |	Barcode consensus sequences from 570 sample benchmarking dataset, run with no contaminant sequences, r 1.5 and s 100 MGE parameters, using the 'standard' (i.e. concat version of MGE snakemake workflow) |
mge_standard_r15s50_nocontam.fasta |	Barcode consensus sequences from 570 sample benchmarking dataset, run with no contaminant sequences, r 1.5 and s 50 MGE parameters, using the 'standard' (i.e. concat version of MGE snakemake workflow) |
mge_standard_r1s100_nocontam.fasta |	Barcode consensus sequences from 570 sample benchmarking dataset, run with no contaminant sequences, r 1 and s 100 MGE parameters, using the 'standard' (i.e. concat version of MGE snakemake workflow) |
mge_standard_r1s50_nocontam.fasta |	Barcode consensus sequences from 570 sample benchmarking dataset, run with no contaminant sequences, r 1 and s 50 MGE parameters, using the 'standard' (i.e. concat version of MGE snakemake workflow) |
concat.sh | Script to create 'concatenated_untriaged.fasta' |
concat_tsv.py | Creates 'concatenated_untriaged.tsv' |
concatenated_untriaged.fasta | Contains 6,840 barcode consensus sequences, I think compiled from 'mge_fastp|standard_*' benchmarking dataset MGE runs? |
concatenated_untriaged.tsv | Combined TSV summary statistics file for 'mge_fastp|standard_*' benchmarking dataset MGE runs |

- For each FASTA file containing barcode consensus seqeunces, a corresponding barcode validator log and TSV file output by barcode validator exist.
- For each 'mge_fastp|standard_*' FASTA file, a corresponding MGE run YAML file and CSV file output by MGE snakemake workflow containing summary stats also exist.
