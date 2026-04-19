#!/usr/bin/env python

import argparse
import csv
import pandas as pd
from pandas.errors import EmptyDataError, ParserError
import numpy as np
import re
import matplotlib.pyplot as plt
import plotly.graph_objects as go


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Calculate depth in gene and pseudogene regions")
    parser.add_argument("--class_tsv", required=True, help="Input read classification table (TSV)")
    parser.add_argument("--positions", required=True, help="Input positions file (TSV)")
    parser.add_argument("--report", required=True, help="Output classification report")
    parser.add_argument("--output_coverage", required=True, help="Output coverage file")
    parser.add_argument("--output_plot", required=True, help="Output plot file")
    return parser.parse_args()


def load_classification(class_tsv):
    """Load and count read classifications."""
    df = pd.read_csv(class_tsv, sep="\t")
    if df.empty or df.shape[1] < 6:
        classification_counts = pd.Series([0, 0, 0], index=["Gene", "Pseudogene", "Ambiguous"])
    elif "Classification" in df.columns:
        classification_counts = df["Classification"].value_counts().reindex(["Gene", "Pseudogene", "Ambiguous"]).fillna(0).astype(int)
    else:
        classification_counts = (
            df.iloc[:, 5].value_counts().reindex(["Gene", "Pseudogene", "Ambiguous"]).fillna(0).astype(int)
        )
    report_df = classification_counts.reset_index()
    report_df.columns = ["Classification", "Count"]
    return df, report_df


def write_empty_outputs(output_report, output_coverage, output_plot):
    pd.DataFrame(columns=["Classification", "Count"]).to_csv(output_report, index=False)
    pd.DataFrame(
        columns=[
            "Gene_Position",
            "Gene_Coverage",
            "Gene_Base",
            "Pseudogene_Base",
            "Ambiguous_Coverage",
            "Pseudogene_Position",
            "Pseudogene_Coverage",
        ]
    ).to_csv(output_coverage, sep="\t", index=False)
    plt.figure(figsize=(8, 4))
    plt.title("No coverage data available")
    plt.xlabel("Position")
    plt.ylabel("Coverage")
    plt.tight_layout()
    plt.savefig(output_plot)
    plt.close()


def calculate_coverage(df, positions):
    output_df = pd.DataFrame(
        columns=["Gene_Position", "Gene_Coverage", "Ambiguous_Coverage", "Pseudogene_Position", "Pseudogene_Coverage"]
    )
    for _, row in positions.iterrows():
        gene_position = row.iloc[0]
        pseudogene_position = row.iloc[2]
        read_bases = row.iloc[5:]
        if gene_position == "-":
            continue

        covered_reads = read_bases[read_bases != "-"].index
        classifications = df.loc[covered_reads, "Classification"]

        gene_cov = (classifications == "Gene").sum()
        amb_cov = (classifications == "Ambiguous").sum()
        pseudo_cov = (classifications == "Pseudogene").sum()

        new_row = pd.DataFrame(
            {
                "Gene_Position": [gene_position],
                "Gene_Coverage": [gene_cov],
                "Gene_Base": [row["Gene_Base"]],
                "Pseudogene_Base": [row["Pseudogene_Base"]],
                "Ambiguous_Coverage": [amb_cov],
                "Pseudogene_Position": [pseudogene_position],
                "Pseudogene_Coverage": [pseudo_cov],
            }
        )
        output_df = pd.concat([output_df, new_row], ignore_index=True)
    return output_df


def generate_report(report_df, output_report):
    report_df.to_csv(output_report, index=False)


def main():
    """
    Main function to execute the script.

    This function performs the following steps:
    1. Parses command-line arguments.
    2. Loads classification data from a TSV file and generates a report.
    3. Loads positions data from a CSV file and calculates coverage.
    4. Prepares data for plotting.
    5. Generates a stacked bar plot of coverage per position and saves it as an image.

    The function expects the following command-line arguments:
    - class_tsv: Path to the TSV file containing classification data.
    - report: Path to save the generated report.
    - positions: Path to the CSV file containing positions data.
    - output_coverage: Path to save the calculated coverage data.
    - output_plot: Path to save the generated plot image.

    The plot displays coverage for 'gene', 'ambiguous', and 'pseudogene' classifications
    for each position, with regions of interest (ROI) extracted from the class_tsv filename.
    """
    """Main function to execute the script."""
    args = parse_arguments()

    # Step 1: Load classification and generate report
    try:
        df, report_df = load_classification(args.class_tsv)
    except (EmptyDataError, ParserError, FileNotFoundError):
        write_empty_outputs(args.report, args.output_coverage, args.output_plot)
        return

    generate_report(report_df, args.report)
    if df.empty or "Read_Name" not in df.columns:
        write_empty_outputs(args.report, args.output_coverage, args.output_plot)
        return
    df.set_index("Read_Name", inplace=True)

    # Step 2: Load positions and calculate coverage
    try:
        positions_df = pd.read_csv(args.positions, sep="\t")
    except EmptyDataError:
        write_empty_outputs(args.report, args.output_coverage, args.output_plot)
        return

    if positions_df.empty:
        write_empty_outputs(args.report, args.output_coverage, args.output_plot)
        return

    output_df = calculate_coverage(df, positions_df)
    output_df.to_csv(args.output_coverage, sep="\t", index=False)

    if output_df.empty:
        write_empty_outputs(args.report, args.output_coverage, args.output_plot)
        return

    # Prepare data
    positions = output_df["Gene_Position"]
    gene_cov = output_df["Gene_Coverage"]
    amb_cov = output_df["Ambiguous_Coverage"]
    pseudo_cov = output_df["Pseudogene_Coverage"]

    x = np.arange(len(positions))
    width = 0.25

    pattern = r"(\d+:\d+-\d+)_(\d+:\d+-\d+)"
    match = re.search(pattern, args.class_tsv)
    region1, region2 = "NA", "NA"
    if match:
        region1, region2 = match.groups()

    fig = go.Figure(
        data=[
            go.Bar(name="Gene", x=positions, y=gene_cov),
            go.Bar(name="Ambiguous", x=positions, y=amb_cov),
            go.Bar(name="Pseudogene", x=positions, y=pseudo_cov),
        ]
    )

    fig.update_layout(
        barmode="overlay",
        title={
            "text": f"Coverage per Position on 'gene' ROI {region1}<br>and 'pseudogene' ROI {region2}",
            "y": 0.9,
            "x": 0.5,
            "xanchor": "center",
            "yanchor": "top",
        },
        xaxis=dict(title="Position", tickangle=45, dtick=50),
        yaxis=dict(title="Coverage", range=[0, max(gene_cov.max(), amb_cov.max(), pseudo_cov.max()) + 10]),
        legend=dict(title="Classification"),
    )

    fig.write_image(args.output_plot)


if __name__ == "__main__":
    main()