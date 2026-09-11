#!/usr/bin/env python3

import argparse
import json
import re
import sys


## FUNCTIONS TO PARSE RELEVANT QC METRICS FROM ISOSEQ PREPROCESSING SUMMARY REPORTS 
def parse_ccs_report(ccs_path):
    """
    Parse key QC metrics from the IsoSeq CCS report (.TXT)
    Input: Filepath to single SMRT cell CCS report
    Returns: A dictionary of relevant count metrics 
    """
    metrics = {}

    # Use regex to parse relevant sections
    patterns = {
        "ZMWs_input": r"^ZMWs input\s*:\s*(\d+)",
        "ZMWs_pass": r"^ZMWs pass filters\s*:\s*(\d+)",
        "ZMWs_fail": r"^ZMWs fail filters\s*:\s*(\d+)",
        "Fail_lacking_full_passes": r"^Lacking full passes\s*:\s*(\d+)",
        "Fail_below_min_RQ": r"^CCS below minimum RQ\s*:\s*(\d+)",
        "Fail_coverage_drop": r"^Coverage drops \s*:\s*(\d+)",
        "Fail_draft_generation": r"^Draft generation error\s*:\s*(\d+)",
    }

    # Extract relevant sections to the metrics dictionary 
    with open(ccs_path) as fh:
        for line in fh:
            line = line.strip()
            for key, pattern in patterns.items():
                match = re.match(pattern, line)
                if match:
                    metrics[key] = int(match.group(1))

    # Check if all metrics are present 
    required = ["ZMWs_input","ZMWs_pass", "ZMWs_fail", "Fail_lacking_full_passes", "Fail_below_min_RQ", "Fail_coverage_drop", "Fail_draft_generation"]
    for key in required:
        if key not in metrics:
            raise ValueError(f"Missing CCS metric: {key}")
        
    return metrics


def parse_refine_report(refine_path):
    """
    Parse key QC metrics from the IsoSeq Refine report (.JSON)
    Input: Filepath to single SMRT cell Refine report
    Returns: A dictionary of relevant count metrics 
    """
    with open(refine_path) as fh:
        data = json.load(fh)

    metrics = {}

    # Parse relevant sections of the JSON file
    for attr in data.get("attributes", []):
        if attr["id"] == "num_reads_fl":
            metrics["FL_reads"] = int(attr["value"])
        elif attr["id"] == "num_reads_flnc":
            metrics["FLNC_reads"] = int(attr["value"])
        elif attr["id"] == "num_reads_flnc_polya":
            metrics["FLNC_polyA_reads"] = int(attr["value"])

    # Check if all requiredmetrics are present 
    required = ["FL_reads", "FLNC_reads", "FLNC_polyA_reads"]
    for key in required:
        if key not in metrics:
            raise ValueError(f"Missing refine metric: {key}")

    return metrics


## FUNCTION TO CALCULATE DERIVED QC METRICS
def pct_calc(numerator, denominator):
    """Percentage helper, returns float."""
    if denominator == 0:
        return 0.0
    return (numerator / denominator) * 100.0


## MAIN SCRIPT 
def main():
    parser = argparse.ArgumentParser(
        description="Extract Iso-Seq QC metrics from SMRT cell preprocessing reports"
    )
    parser.add_argument("--smrt-cell", required=True, help="SMRT cell ID")
    parser.add_argument("--ccs-report", required=True, help="Isoseq CCS report (.TXT)")
    parser.add_argument("--refine-report", required=True, help="Isoseq Refine report (.JSON)")

    args = parser.parse_args()

    # Run the above QC parser functions 
    ccs = parse_ccs_report(args.ccs_report)
    ref = parse_refine_report(args.refine_report)

    # Calculate derived metrics (percentages)
    ZMWs_pass_pct = pct_calc(ccs["ZMWs_pass"], ccs["ZMWs_input"])
    ZMWs_fail_pct = pct_calc(ccs["ZMWs_fail"], ccs["ZMWs_input"])

    Fail_full_passes_pct = pct_calc(ccs["Fail_lacking_full_passes"], ccs["ZMWs_fail"])
    Fail_below_RQ_pct = pct_calc(ccs["Fail_below_min_RQ"], ccs["ZMWs_fail"])
    Fail_coverage_pct = pct_calc(ccs["Fail_coverage_drop"], ccs["ZMWs_fail"])
    Fail_draft_gen_pct = pct_calc(ccs["Fail_draft_generation"], ccs["ZMWs_fail"])

    FLNC_rate_pct = pct_calc(ref["FLNC_reads"], ccs["ZMWs_input"])
    FLNC_polyA_pct_input = pct_calc(ref["FLNC_polyA_reads"], ccs["ZMWs_input"])

    # Output TSV row 
    out = [
        args.smrt_cell,
        ccs["ZMWs_input"],
        ccs["ZMWs_pass"],
        f"{ZMWs_pass_pct:.2f}",
        ccs["ZMWs_fail"],
        f"{ZMWs_fail_pct:.2f}",
        ccs["Fail_lacking_full_passes"],
        f"{Fail_full_passes_pct:.2f}",
        ccs["Fail_below_min_RQ"],
        f"{Fail_below_RQ_pct:.2f}",
        ccs["Fail_coverage_drop"],
        f"{Fail_coverage_pct:.2f}",
        ccs["Fail_draft_generation"],
        f"{Fail_draft_gen_pct:.2f}",
        ref["FL_reads"],
        ref["FLNC_reads"],
        f"{FLNC_rate_pct:.2f}",
        ref["FLNC_polyA_reads"],
        f"{FLNC_polyA_pct_input:.2f}",
    ]

    print("\t".join(map(str, out)))


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(f"ERROR: {e}\n")
        sys.exit(1)
