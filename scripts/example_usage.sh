#!/bin/bash

# Example usage of the histotab-format CLI
# This script shows how to run the CLI with the TCGA LUAD data

# Set the base data directory
DATA_DIR="/home/val/workspaces/histotab/data/raw/TCGA_LUAD"

# Run the CLI command
python histotab/data/cli_format_tcga_tab.py \
    "$DATA_DIR/luad_tcga_pan_can_atlas_2018/data_clinical_patient.txt" \
    "$DATA_DIR/luad_tcga_firehose_legacy/data_clinical_patient.txt" \
    "$DATA_DIR/luad_tcga_gdc/data_clinical_patient.txt" \
    "$DATA_DIR/TCGA_LUAD_PanCancer_gene_oncoprint.tsv" \
    "/home/val/workspaces/histotab/data/processed/tabular_data/tcga_luad_pancancer_combined_clinical_genes_cli.csv" \
    --morphology-col-pancancer "ICD_O_3_HISTOLOGY" \
    --morphology-col-legacy "ICD_O_3_HISTOLOGY" \
    --morphology-col-gdc "MORPHOLOGY"
