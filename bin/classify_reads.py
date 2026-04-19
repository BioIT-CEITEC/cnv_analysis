#!/usr/bin/env python

import pandas as pd
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def classify_quantile(gene_support, pseudogene_support, ambiguous_support, gene_threshold, pseudogene_threshold):
    if gene_support >= gene_threshold and gene_support > pseudogene_support and gene_support > ambiguous_support:
        return "Gene"
    elif (
        pseudogene_support >= pseudogene_threshold
        and pseudogene_support > gene_support
        and pseudogene_support > ambiguous_support
    ):
        return "Pseudogene"
    else:
        return "Ambiguous"


def classify(
    gene_support,
    pseudogene_support,
    ambiguous_support,
    gene_quantile,
    pseudogene_quantile,
    total_positions,
    threshold,
    metric,
    content_length,
    gene_mean,
    pseudogene_mean,
):
    if metric == "threshold":
        # Compute ratios and classify
        if total_positions > 0:
            gene_ratio = gene_support * (pseudogene_mean / gene_mean) / content_length
            pseudogene_ratio = pseudogene_support * (gene_mean / pseudogene_mean) / content_length
            if gene_ratio >= threshold:
                return "Gene"
            elif pseudogene_ratio >= threshold or pseudogene_support > gene_support:
                return "Pseudogene"
            else:
                return "Ambiguous"
        else:
            return "Ambiguous"

    elif metric == "majority":
        if gene_support > pseudogene_support and gene_support > ambiguous_support:
            return "Gene"
        elif pseudogene_support > gene_support and pseudogene_support > ambiguous_support:
            return "Pseudogene"
        else:
            return "Ambiguous"

    elif metric == "quantile":
        return classify_quantile(
            gene_support, pseudogene_support, ambiguous_support, gene_quantile, pseudogene_quantile
        )


def plot_classification_results(classification_df, tsv_file, output):
    plt.figure(figsize=(15, 5))
    plt.suptitle(f"Classification Results for {tsv_file}")

    plt.subplot(1, 3, 1)
    sns.histplot(classification_df["Gene_Support"], kde=True, bins=30, color="blue")
    plt.title("Gene Support Distribution")

    plt.subplot(1, 3, 2)
    sns.histplot(classification_df["Pseudogene_Support"], kde=True, bins=30, color="green")
    plt.title("Pseudogene Support Distribution")

    plt.subplot(1, 3, 3)
    sns.histplot(classification_df["Ambiguous_Support"], kde=True, bins=30, color="red")
    plt.title("Ambiguous Support Distribution")

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(output.replace(".tsv", "_histogram.png"))

    sns.pairplot(classification_df[["Gene_Support", "Pseudogene_Support", "Ambiguous_Support"]])
    plt.savefig(output.replace(".tsv", "_pairplot.png"))

    corr = classification_df[["Gene_Support", "Pseudogene_Support", "Ambiguous_Support"]].corr()
    sns.heatmap(corr, annot=True, cmap="coolwarm")
    plt.title("Correlation Matrix")
    plt.savefig(output.replace(".tsv", "_correlation.png"))

    # Add print statement
    print(f"Classification results saved to {output}")


def classify_reads(tsv_file, output, threshold=0.7, metric="threshold"):
    # Load the TSV file
    df = pd.read_csv(tsv_file, sep="\t")

    # List of all columns starting from the 6th column onward are the reads
    read_columns = df.columns[5:]

    # Initialize a results list to store classifications
    supports = []

    # Iterate over each read column
    for read in read_columns:
        # Initialize counts
        gene_support = 0
        pseudogene_support = 0
        ambiguous_support = 0
        read_bases = df[read]

        non_gap_positions = read_bases[read_bases != "-"].index
        total_positions = non_gap_positions.size
        content_length = len(non_gap_positions)
        if non_gap_positions.empty or content_length == 0 or total_positions == 0:
            continue

        first_idx = non_gap_positions[0]
        last_idx = non_gap_positions[-1]

        # For positions where Difference_Type != "M" and indices between first_idx and last_idx
        for idx in range(first_idx, last_idx + 1):
            row = df.iloc[idx]
            difference_type = row["Difference_Type"]
            if difference_type != "M":
                total_positions += 1
                read_base = row[read]
                gene_base = row["Gene_Base"]
                pseudogene_base = row["Pseudogene_Base"]

                if read_base == "-":
                    if gene_base == "-" and pseudogene_base != "-":
                        gene_support += 1
                    elif pseudogene_base == "-" and gene_base != "-":
                        pseudogene_support += 1
                    else:
                        ambiguous_support += 1
                else:
                    if gene_base == pseudogene_base:
                        ambiguous_support += 1
                    else:
                        if read_base == gene_base and read_base != pseudogene_base:
                            gene_support += 1
                        elif read_base == pseudogene_base and read_base != gene_base:
                            pseudogene_support += 1
                        else:
                            ambiguous_support += 1

        # # Get starting genomic coordinates of the read
        # if classification == "Gene" or classification == "Ambiguous":
        #     start_coordinate = df.iloc[first_idx]["Gene_Position"]
        # elif classification == "Pseudogene":
        #     start_coordinate = df.iloc[first_idx]["Pseudogene_Position"]
        # else:
        #     start_coordinate = None

        supports.append(
            {
                "Read_Name": read,
                "Gene_Support": gene_support,
                "Pseudogene_Support": pseudogene_support,
                "Ambiguous_Support": ambiguous_support,
                "Total_Positions": total_positions,
                "Gene_Start_Coordinate": df.iloc[first_idx]["Gene_Position"],
                "Pseudo_Start_Coordinate": df.iloc[first_idx]["Pseudogene_Position"],
            }
        )

    # Convert the results to a DataFrame
    support_df = pd.DataFrame(supports)
    if support_df.empty:
        # Keep a stable schema even when no reads passed filtering.
        empty_columns = [
            "Read_Name",
            "Gene_Support",
            "Pseudogene_Support",
            "Ambiguous_Support",
            "Total_Positions",
            "Gene_Start_Coordinate",
            "Pseudo_Start_Coordinate",
            "Classification",
            "Read_Start",
        ]
        pd.DataFrame(columns=empty_columns).to_csv(output, sep="\t", index=False)
        return
    print(support_df.head())
    plot_classification_results(support_df, tsv_file, output)

    gene_threshold = support_df["Gene_Support"].quantile(0.75)
    pseudogene_threshold = support_df["Pseudogene_Support"].quantile(0.75)
    gene_mean = support_df["Gene_Support"].mean()
    pseudogene_mean = support_df["Pseudogene_Support"].mean()
    classification_df = support_df.copy()
    classification_df["Classification"] = classification_df.apply(
        lambda x: classify(
            x["Gene_Support"],
            x["Pseudogene_Support"],
            x["Ambiguous_Support"],
            gene_threshold,
            pseudogene_threshold,
            x["Total_Positions"],
            threshold,
            metric,
            content_length,
            gene_mean=gene_mean,
            pseudogene_mean=pseudogene_mean,
        ),
        axis=1,
    )
    # Set Read_Start based on Classification
    classification_df["Read_Start"] = classification_df.apply(
        lambda x: (
            x["Gene_Start_Coordinate"] if x["Classification"] in ["Gene", "Ambiguous"] else x["Pseudo_Start_Coordinate"]
        ),
        axis=1,
    )
    # Drop the Gene_Start_Coordinate and Pseudogene_Start_Coordinate columns
    classification_df.drop(columns=["Gene_Start_Coordinate", "Pseudo_Start_Coordinate"], inplace=True)

    print(classification_df.head())
    # Save the classification results to the output file
    classification_df.to_csv(output, sep="\t", index=False)
    print(classification_df.head())


def main():
    parser = argparse.ArgumentParser(description="Classify reads from a TSV file.")
    parser.add_argument("--tsv", "-t", required=True, help="Path to the TSV file.")
    parser.add_argument("--output", "-o", required=True, help="Path to the output file.")
    parser.add_argument(
        "--threshold", "-th", required=False, type=float, default=0.7, help="Threshold for classification."
    )
    parser.add_argument(
        "--metric",
        "-m",
        required=False,
        type=str,
        choices=["threshold", "quantile", "majority"],
        default="threshold",
        help="Metric for classification.",
    )
    args = parser.parse_args()

    classify_reads(args.tsv, args.output, threshold=args.threshold, metric=args.metric)


if __name__ == "__main__":
    main()