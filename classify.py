#!/opt/conda/bin/python3
import argparse
import os
import re
import subprocess
import sys
import time


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", type=str, default="krakenuniq", help="Change method (yet to implement for Bignorm+Krakenuniq or custom method)")
    parser.add_argument("--quick", action="store_true", help="Run krakenuniq with --quick flag")
    parser.add_argument("--file", type=str, help="Single ground truth file to evaluate")
    args = parser.parse_args()

    # PATHS
    DATA_DIR = "/mnt/data/"
    DATASET_NAME = "IMMSA"
    FASTA_DIR = os.path.join(DATA_DIR, DATASET_NAME)
    DB_DIR = os.path.join(DATA_DIR, "standard_db")
    GT_GENUS_DIR = os.path.join(FASTA_DIR, "truth_sets", "genus")
    OUTPUT_DIR = os.path.join(DATASET_NAME+"_evaluation", args.method + "_quick" if args.quick else "")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    if not os.path.isdir(FASTA_DIR):
        print(f"Error: {FASTA_DIR} is not a valid directory")
        return
    
    # If --file is given set as only file to make krakenuniq
    gt_files = sorted(os.listdir(GT_GENUS_DIR)) 
    if args.file and args.file in gt_files:
        gt_files = [args.file]

    for gt_filename in gt_files:
        # Skip temporary files
        if gt_filename.startswith("."):
            continue
        
        for fasta_filename in sorted(os.listdir(FASTA_DIR)):
            if fasta_filename.startswith(gt_filename.replace("_TRUTH.txt", "")):
                fasta_fp = os.path.join(FASTA_DIR, fasta_filename)
                if os.path.isfile(fasta_fp):
                    if "R2" in fasta_filename and os.path.isfile(fasta_fp.replace("R2", "R1")):
                        continue # avoid double evaluating the paired read files
                    if "R1" in fasta_filename and os.path.isfile(fasta_fp.replace("R1", "R2")):
                        dataset = gt_filename.replace("_TRUTH.txt", "")
                        report_fp = os.path.join(OUTPUT_DIR, dataset+"_REPORT.tsv")
                        if (os.path.isfile(report_fp)): # skip already classified read files
                            continue
                        cmd = [
                            "krakenuniq",
                            "--paired",
                            "--check-names",
                            "--db", DB_DIR,
                            "--threads", "8",
                            "--report", report_fp,
                            "--output", report_fp.replace("REPORT", "OUTPUT")
                        ]
                        if args.quick:
                            cmd.append("--quick")
                        cmd.extend([
                            fasta_fp,
                            fasta_fp.replace("R1", "R2")
                        ])
                        print(" ".join(cmd))
                        start_time = time.time()
                        subprocess.run(cmd)
                        runtime = time.time() - start_time
                    else:
                        dataset = gt_filename.replace("_TRUTH.txt", "")
                        report_fp = os.path.join(OUTPUT_DIR, re.sub(r"\.(fq|fast[aq])\.gz$", "", fasta_filename) + "_REPORT.tsv")
                        if (os.path.isfile(report_fp)): # skip already classified read files
                            continue
                        cmd = [
                            "krakenuniq",
                            "--db", DB_DIR,
                            "--threads", "8",
                            "--report", report_fp,
                            "--output", report_fp.replace("REPORT", "OUTPUT")
                        ]
                        if args.quick:
                            cmd.append("--quick")
                        cmd.append(fasta_fp)
                        print(" ".join(cmd))
                        start_time = time.time()
                        subprocess.run(cmd)
                        runtime = time.time() - start_time

     
if __name__ == '__main__':
    main()
    