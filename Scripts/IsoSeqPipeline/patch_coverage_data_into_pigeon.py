#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

import pandas as pd


# Define target columns 
target_cols = ["min_sample_cov", "min_cov", "min_cov_pos", "sd_cov"]


## MAIN SCRIPT 
def main():

    # Set argument parser 
    parser = argparse.ArgumentParser(
        description="Patch junction coverage columns from SQANTI3 classification into Pigeon's classification file"
    )
    parser.add_argument("--pigeon-classify", required=True, help="Input Pigeon classification file (.TXT)")
    parser.add_argument("--sqanti-classify", required=True, help="Input SQANTI3 classification file (.TXT)")
    parser.add_argument("--output", required=True, help="Output path for the patched classification file (.TXT)")

    args = parser.parse_args()

    # Run the above QC parser functions 
    pigeon_classify = Path(args.pigeon_classify)
    sqanti_classify = Path(args.sqanti_classify)
    patched_output = Path(args.output)

    # Read in the classification files
    print("Reading Pigeon dataset...")
    df_pigeon = pd.read_csv(pigeon_classify, sep="\t", dtype=str, keep_default_na=False)

    print("Reading SQANTI3 dataset...")
    df_sqanti = pd.read_csv(sqanti_classify, sep="\t", dtype=str, keep_default_na=False)


    # Map junction metrics from SQANTI to Pigeon 
    print("Mapping SQANTI3 metrics directly onto Pigeon column positions...")
    
    for col in target_cols:
        if col in df_sqanti.columns:

            # Build an exact dictionary lookup based on Transcript ID
            mapping_dict = dict(zip(df_sqanti["isoform"], df_sqanti[col]))
            
            print(f"Injecting values for '{col}' without moving column...")
            # If a transcript was dropped by Pigeon, it's ignored. 
            # If it's missing in SQANTI3, it retains Pigeon's original value (which is "NA")
            df_pigeon[col] = df_pigeon["isoform"].map(mapping_dict).fillna(df_pigeon[col])
        
        else:
            print(f"Warning: Column '{col}' not found in SQANTI3 file. Skipping in-place update.")

    # Integrity verification
    print("\n--- Structural Alignment Complete ---")
    print(f"Final Table Structure: {df_pigeon.shape[0]:,} rows x {df_pigeon.shape[1]} columns")

    # Save back out, maintaining the pristine layout Pigeon expects
    df_pigeon.to_csv(patched_output, sep="\t", index=False)
    print(f"\nSuccessfully created patched Pigeon file: {patched_output}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(f"ERROR: {e}\n")
        sys.exit(1)
