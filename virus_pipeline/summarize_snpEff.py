import sys
import pandas as pd
import glob
import os
import json
import argparse
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', type=str, required=True, help='Path to input directory containing SnpEff summary files.')
    parser.add_argument('--output_dir', type=str, required=True, help='Path to output directory for summary files.')
    args = parser.parse_args(argv)

    input_dir = args.input_dir
    output_dir = args.output_dir

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # File pattern
    file_pattern = "*_snpEff_summary.genes.txt"

    # Expected columns to keep
    columns_to_keep = [
        "TranscriptId",
        "variants_impact_HIGH",
        "variants_impact_LOW",
        "variants_impact_MODERATE",
        "variants_impact_MODIFIER",
        "variants_effect_downstream_gene_variant",
        "variants_effect_frameshift_variant",
        "variants_effect_missense_variant",
        "variants_effect_splice_region_variant",
        "variants_effect_synonymous_variant"
    ]

    # Full column list for files with headers
    full_columns = [
        "GeneName", "GeneId", "TranscriptId", "BioType",
        "variants_impact_HIGH", "variants_impact_LOW", "variants_impact_MODERATE", "variants_impact_MODIFIER",
        "variants_effect_conservative_inframe_deletion",
        "variants_effect_disruptive_inframe_insertion",
        "variants_effect_downstream_gene_variant", "variants_effect_frameshift_variant",
        "variants_effect_missense_variant", "variants_effect_splice_region_variant",
        "variants_effect_stop_gained", "variants_effect_synonymous_variant", "variants_effect_upstream_gene_variant"
    ]

    # Initialize data list
    all_data = []

    # Check if directory exists
    if not os.path.exists(input_dir):
        logging.error(f"Directory '{input_dir}' does not exist.")
        sys.exit(1)

    # Get list of files
    files = glob.glob(os.path.join(input_dir, file_pattern))
    if not files:
        logging.error(f"No files found matching pattern '{file_pattern}' in '{input_dir}'.")
        logging.error(f"Files in directory: {os.listdir(input_dir)}")
        sys.exit(1)

    # Iterate over files
    for file_path in files:
        sample_name = os.path.basename(file_path).replace("_snpEff_summary.genes.txt", "")
        logging.info(f"Processing file: {file_path} (Sample: {sample_name})")
        
        try:
            # Read first few lines for debugging
            with open(file_path, 'r') as f:
                logging.debug(f"First 3 lines of {file_path}:")
                lines = []
                for i in range(3):
                    try:
                        line = next(f).strip()
                        lines.append(line)
                        logging.debug(f"Line {i+1}: {line}")
                        fields = line.split("\t")
                        logging.debug(f"Fields in line {i+1}: {len(fields)} (Expected: 13, 14, or 15)")
                        if i > 0 and len(fields) not in [13, 14, 15]:
                            logging.warning(f"Incorrect field count in line {i+1}")
                    except StopIteration:
                        break
            
            # Try reading with tab delimiter
            try:
                df = pd.read_csv(file_path, sep="\t", comment="#", header=1)
            except Exception as e:
                logging.warning(f"Tab delimiter failed: {str(e)}. Trying comma delimiter...")
                try:
                    df = pd.read_csv(file_path, sep=",", comment="#", header=1)
                except Exception as e:
                    logging.error(f"Comma delimiter failed: {str(e)}. Skipping file.")
                    continue
            
            logging.info(f"Columns in {file_path}: {list(df.columns)}")
            
            if len(df.columns) not in [13, 14, 15]:
                logging.error(f"Column count mismatch in {file_path}. Expected 13, 14, or 15, found {len(df.columns)}")
                continue
            
            expected_cols = set(full_columns)
            actual_cols = set(df.columns)
            if not any(col in actual_cols for col in columns_to_keep):
                logging.error(f"No expected columns found in {file_path}. Attempting manual column assignment.")
                if len(df.columns) in [13, 14, 15]:
                    logging.info(f"Assigning columns: {full_columns[:len(df.columns)]}")
                    df.columns = full_columns[:len(df.columns)]
                else:
                    continue
            
            for col in columns_to_keep:
                if col not in df.columns:
                    logging.warning(f"Column {col} missing in {file_path}. Filling with zeros.")
                    df[col] = 0
            
            missing_cols = [col for col in columns_to_keep if col not in df.columns]
            if missing_cols:
                logging.error(f"Missing columns in {file_path}: {missing_cols}")
                continue
            
            df = df[columns_to_keep].copy()
            df["Sample"] = sample_name
            all_data.append(df)
        except Exception as e:
            logging.error(f"Error reading {file_path}: {str(e)}")
            continue

    if not all_data:
        logging.error("No valid dataframes to concatenate. Check file contents or column names.")
        sys.exit(1)

    combined_df = pd.concat(all_data, ignore_index=True)
    summary_df = combined_df.groupby("TranscriptId").agg({
        "variants_impact_HIGH": "sum",
        "variants_impact_LOW": "sum",
        "variants_impact_MODERATE": "sum",
        "variants_impact_MODIFIER": "sum",
        "variants_effect_downstream_gene_variant": "sum",
        "variants_effect_frameshift_variant": "sum",
        "variants_effect_missense_variant": "sum",
        "variants_effect_splice_region_variant": "sum",
        "variants_effect_synonymous_variant": "sum",
        "Sample": lambda x: ";".join(x)
    }).reset_index()

    summary_df.columns = [
        "TranscriptId",
        "Total_variants_impact_HIGH",
        "Total_variants_impact_LOW",
        "Total_variants_impact_MODERATE",
        "Total_variants_impact_MODIFIER",
        "Total_variants_effect_downstream_gene_variant",
        "Total_variants_effect_frameshift_variant",
        "Total_variants_effect_missense_variant",
        "Total_variants_effect_splice_region_variant",
        "Total_variants_effect_synonymous_variant",
        "Samples"
    ]

    output_file = os.path.join(output_dir, "summary_table.csv")
    summary_df.to_csv(output_file, index=False)
    logging.info(f"Summary table saved to {output_file}")

    chart_data = {
        "labels": summary_df["TranscriptId"].tolist(),
        "datasets": [
            {
                "label": "High Impact",
                "data": summary_df["Total_variants_impact_HIGH"].tolist(),
                "backgroundColor": "rgba(255, 99, 132, 0.7)"
            },
            {
                "label": "Low Impact",
                "data": summary_df["Total_variants_impact_LOW"].tolist(),
                "backgroundColor": "rgba(54, 162, 235, 0.7)"
            },
            {
                "label": "Moderate Impact",
                "data": summary_df["Total_variants_impact_MODERATE"].tolist(),
                "backgroundColor": "rgba(255, 206, 86, 0.7)"
            },
            {
                "label": "Modifier Impact",
                "data": summary_df["Total_variants_impact_MODIFIER"].tolist(),
                "backgroundColor": "rgba(75, 192, 192, 0.7)"
            }
        ]
    }

    chart_json_file = os.path.join(output_dir, "chart_data.json")
    with open(chart_json_file, "w") as f:
        json.dump(chart_data, f, indent=2)
    logging.info(f"Chart data for external visualization saved to {chart_json_file}")

if __name__ == '__main__':
    main()