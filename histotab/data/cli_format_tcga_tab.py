"""
CLI for formatting TCGA LUAD data combining clinical and gene mutation data.
"""

from pathlib import Path

import click
import pandas as pd


def create_histology_mappings():
    """Create the histology and pattern mappings."""
    histology_map = {
        "8140/3": "Adenocarcinoma, NOS",
        "8255/3": "Adenocarcinoma with mixed subtypes",
        "8260/3": "Papillary adenocarcinoma, NOS",
        "8550/3": "Acinar cell carcinoma",
        "8480/3": "Mucinous adenocarcinoma",
        "8310/3": "Clear cell adenocarcinoma, NOS",
        "8252/3": "Bronchiolo-alveolar carcinoma, non-mucinous",
        "8253/3": "Invasive mucinous adenocarcinoma",
        "8230/3": "Solid carcinoma, NOS",
        "8507/3": "Invasive micropapillary carcinoma",
        "8250/3": "Lepidic adenocarcinoma",
        "8490/3": "Signet ring cell carcinoma",
    }

    pattern_map = {
        "Lepidic adenocarcinoma": "Lepidic",
        "Bronchiolo-alveolar carcinoma, non-mucinous": "Lepidic",  # older term
        "Acinar cell carcinoma": "Acinar",
        "Papillary adenocarcinoma, NOS": "Papillary",
        "Solid carcinoma, NOS": "Solid",
        "Invasive micropapillary carcinoma": "Micropapillary",
        "Mucinous adenocarcinoma": "To drop",  # not part of 5 canonical patterns
        "Invasive mucinous adenocarcinoma": "To drop",  # not part of 5 canonical patterns
        "Clear cell adenocarcinoma, NOS": "To drop",  # not part of 5 canonical patterns
        "Signet ring cell carcinoma": "To drop",  # not part of 5 canonical patterns
        "Adenocarcinoma with mixed subtypes": "Mixed",  # optionally drop or keep as its own group
        "Adenocarcinoma, NOS": None,  # too vague
    }

    return histology_map, pattern_map


def process_clinical_data(
    clinical_file: Path, morphology_col: str, histology_map: dict, pattern_map: dict
):
    """Process clinical data by adding histologic subtype and LUAD pattern."""
    df_clinical = pd.read_csv(clinical_file, sep="\t", comment="#", skip_blank_lines=True)

    # Add histologic subtype and LUAD pattern
    df_clinical["HISTOLOGIC_SUBTYPE"] = df_clinical[morphology_col].map(histology_map)
    df_clinical["LUAD_PATTERN"] = df_clinical["HISTOLOGIC_SUBTYPE"].map(pattern_map)

    # Filter out patterns marked as "To drop"
    df_clinical = df_clinical[~df_clinical["LUAD_PATTERN"].isin(["To drop"])]

    return df_clinical


def process_gene_data(gene_file: Path):
    """Process gene mutation data from oncoprint file."""
    df = pd.read_csv(gene_file, sep="\t", index_col=0)

    # Filter only mutation or CNA rows
    df_filtered = df[df["track_type"].isin(["MUTATIONS", "CNA"])]

    # Drop the 'track_type' column
    df_filtered = df_filtered.drop(columns="track_type")

    # Group by gene name and mark as 'MUT' if any value exists for a patient
    df_binary = (
        df_filtered.groupby(df_filtered.index).apply(lambda gene_df: gene_df.notnull().any()).T
    )

    # Convert boolean to 'MUT'/'WT'
    df_binary = df_binary.replace({True: "MUT", False: "WT"})

    # Rename index to patient_id
    df_binary.index.name = "patient_id"

    gene_columns = df_binary.columns.tolist()

    return df_binary, gene_columns


def merge_and_complete_data(df_clinical_main, df_clinical_legacy, df_clinical_gdc, df_binary):
    """Merge clinical and gene data, then complete missing values from other sources."""

    # Set index to PATIENT_ID for all clinical dataframes
    df_clinical_main = df_clinical_main.set_index("PATIENT_ID")
    df_clinical_legacy = df_clinical_legacy.set_index("PATIENT_ID")
    df_clinical_gdc = df_clinical_gdc.set_index("PATIENT_ID")

    # Merge main clinical data with gene data
    df_merged = df_clinical_main.join(df_binary, how="inner")

    # Define relevant clinical columns
    relevant_columns = [
        "LUAD_PATTERN",
        "HISTOLOGIC_SUBTYPE",
        "AGE",
        "SEX",
        "OS_STATUS",
        "OS_MONTHS",
        "DSS_STATUS",
        "DSS_MONTHS",
        "DFS_STATUS",
        "DFS_MONTHS",
        "PFS_STATUS",
        "PFS_MONTHS",
    ]

    # Check which columns exist in the merged dataframe
    available_clinical_columns = [col for col in relevant_columns if col in df_merged.columns]

    # Combine with mutation data (gene columns come from df_binary)
    final_columns = available_clinical_columns + list(df_binary.columns)
    df_final = df_merged[final_columns]

    # Complete missing values from legacy and GDC datasets
    click.echo("Checking for missing values to complete from legacy and GDC datasets...")

    # For each clinical column, try to fill missing values from other datasets
    for col in available_clinical_columns:
        if col in df_final.columns:
            missing_mask = df_final[col].isna()
            if missing_mask.any():
                click.echo(f"Found {missing_mask.sum()} missing values in column '{col}'")
                if col in df_clinical_gdc.columns:
                    # Filter GDC dataframe to include only matching PATIENT_IDs
                    gdc_values = df_clinical_gdc.loc[
                        df_clinical_gdc.index.intersection(df_final.index), col
                    ]
                    filled_from_gdc = df_final.loc[missing_mask, col].fillna(gdc_values)
                    df_final.loc[missing_mask, col] = filled_from_gdc
                    new_missing = df_final[col].isna().sum()
                    filled_count = missing_mask.sum() - new_missing
                    if filled_count > 0:
                        click.echo(f"  Filled {filled_count} values from GDC dataset")

                # Try to fill remaining missing values from Legacy dataset
                remaining_missing_mask = df_final[col].isna()
                if remaining_missing_mask.any() and col in df_clinical_legacy.columns:
                    # Filter Legacy dataframe to include only matching PATIENT_IDs
                    legacy_values = df_clinical_legacy.loc[
                        df_clinical_legacy.index.intersection(df_final.index), col
                    ]
                    filled_from_legacy = df_final.loc[remaining_missing_mask, col].fillna(
                        legacy_values
                    )
                    df_final.loc[remaining_missing_mask, col] = filled_from_legacy
                    final_missing = df_final[col].isna().sum()
                    filled_count = new_missing - final_missing
                    if filled_count > 0:
                        click.echo(f"  Filled {filled_count} values from Legacy dataset")

                # Check if there are still missing values
                final_missing = df_final[col].isna().sum()
                if final_missing > 0:
                    click.echo(
                        f"  Still {final_missing} missing values in '{col}' after completion"
                    )

    return df_final


@click.command()
@click.argument("clinical_pancancer", type=click.Path(exists=True, path_type=Path))
@click.argument("clinical_legacy", type=click.Path(exists=True, path_type=Path))
@click.argument("clinical_gdc", type=click.Path(exists=True, path_type=Path))
@click.argument("gene_oncoprint", type=click.Path(exists=True, path_type=Path))
@click.argument("output", type=click.Path(path_type=Path))
@click.option(
    "--morphology-col-pancancer",
    default="ICD_O_3_HISTOLOGY",
    help="Morphology column name for PanCancer data",
)
@click.option(
    "--morphology-col-legacy",
    default="ICD_O_3_HISTOLOGY",
    help="Morphology column name for legacy data",
)
@click.option(
    "--morphology-col-gdc", default="MORPHOLOGY", help="Morphology column name for GDC data"
)
def format_data(
    clinical_pancancer,
    clinical_legacy,
    clinical_gdc,
    gene_oncoprint,
    output,
    morphology_col_pancancer,
    morphology_col_legacy,
    morphology_col_gdc,
):
    """
    Format TCGA LUAD data by combining clinical and gene mutation data.

    This command processes clinical data from three sources (PanCancer, legacy, GDC)
    and gene mutation data, then creates a combined dataset with missing value completion.

    Arguments:
        CLINICAL_PANCANCER: Path to PanCancer clinical data file
        CLINICAL_LEGACY: Path to legacy clinical data file
        CLINICAL_GDC: Path to GDC clinical data file
        GENE_ONCOPRINT: Path to gene oncoprint TSV file
        OUTPUT: Output CSV file path
    """

    # Create output directory if it doesn't exist
    output.parent.mkdir(parents=True, exist_ok=True)

    click.echo("Creating histology mappings...")
    histology_map, pattern_map = create_histology_mappings()

    click.echo("Processing PanCancer clinical data...")
    df_clinical_pancancer = process_clinical_data(
        clinical_pancancer, morphology_col_pancancer, histology_map, pattern_map
    )

    click.echo("Processing legacy clinical data...")
    df_clinical_legacy = process_clinical_data(
        clinical_legacy, morphology_col_legacy, histology_map, pattern_map
    )

    click.echo("Processing GDC clinical data...")
    df_clinical_gdc = process_clinical_data(
        clinical_gdc, morphology_col_gdc, histology_map, pattern_map
    )

    click.echo("Processing gene mutation data...")
    df_binary, gene_columns = process_gene_data(gene_oncoprint)

    click.echo("Merging data and completing missing values...")
    df_final = merge_and_complete_data(
        df_clinical_pancancer, df_clinical_legacy, df_clinical_gdc, df_binary
    )

    # Reset index and save
    df_final = df_final.reset_index().rename(columns={"index": "patient_id"})

    click.echo(f"Saving final dataset to {output}...")
    df_final.to_csv(output, index=False)

    click.echo(f"✅ Successfully created combined dataset with shape: {df_final.shape}")
    click.echo(f"   Columns: {list(df_final.columns)}")

    # Show some statistics about LUAD patterns
    if "LUAD_PATTERN" in df_final.columns:
        pattern_counts = df_final["LUAD_PATTERN"].value_counts()
        click.echo("\nLUAD Pattern distribution:")
        for pattern, count in pattern_counts.items():
            click.echo(f"  {pattern}: {count}")

    no_mutation_patients = df_final[df_final[gene_columns].eq("WT").all(axis=1)]["patient_id"]

    if not no_mutation_patients.empty:
        click.echo(
            f"\nPatients with no mutations in the current panel of genes ({len(no_mutation_patients)} patients):"
        )
        for patient_id in no_mutation_patients:
            click.echo(f"  {patient_id}")
    else:
        click.echo("\nAll patients have at least one mutation in the current panel of genes.")


if __name__ == "__main__":
    format_data()
