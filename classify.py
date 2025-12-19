#!/opt/conda/bin/python3
import argparse
import os
import re
import subprocess
import sys
import time

# PATHS
DATA_DIR = "/mnt/data/"
DATASET_NAME = "IMMSA"
FASTA_DIR = os.path.join(DATA_DIR, DATASET_NAME)
DB_DIR = os.path.join(DATA_DIR, "standard_db")
GT_GENUS_DIR = os.path.join(FASTA_DIR, "truth_sets", "genus")


def load_runtimes(runtime_file):
    """Load existing runtimes from file."""
    runtimes = {}
    if not os.path.exists(runtime_file):
        return runtimes
    
    with open(runtime_file, 'r') as f:
        lines = f.readlines()
        if len(lines) <= 1:
            return runtimes
        
        for line in lines[1:]:  # Skip header
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                dataset = parts[0]
                runtime = float(parts[1])
                runtimes[dataset] = runtime
    
    return runtimes


def save_runtimes(runtime_file, runtimes):
    """Save runtimes to file."""
    with open(runtime_file, 'w') as f:
        f.write("dataset\truntime_s\n")
        for dataset in sorted(runtimes.keys()):
            f.write(f"{dataset}\t{runtimes[dataset]:.0f}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--method-suffix", type=str, default="krakenuniq", help="Change method (yet to implement for Bignorm+Krakenuniq or custom method)")
    parser.add_argument("--quick", action="store_true", help="Run krakenuniq with --quick flag")
    parser.add_argument("--preload-size", type=str, help="Run KrakenUniq with --preload-size <500M,1G,2G,...>")
    parser.add_argument("--datasets-to-classify", nargs='+', help="List of datasets (from ground truth naming convention) to classify")
    args = parser.parse_args()

    # Build output directory name based on arguments
    output_dir_name = "krakenuniq"
    if args.method_suffix:
        output_dir_name += f"_{args.method_suffix}"
    if args.quick:
        output_dir_name += "_quick"
    if args.preload_size:
        output_dir_name += f"_preload-size-{args.preload_size}"
    
    OUTPUT_DIR = os.path.join(DATASET_NAME+"_evaluation", output_dir_name)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Load existing runtimes
    runtime_file = os.path.join(OUTPUT_DIR, "runtimes.tsv")
    runtimes = load_runtimes(runtime_file)

    if not os.path.isdir(FASTA_DIR):
        print(f"Error: {FASTA_DIR} is not a valid directory")
        return
    
    gt_files = sorted(os.listdir(GT_GENUS_DIR)) 
    if args.datasets_to_classify:
        gt_files = [f + "_TRUTH.txt" for f in args.datasets_to_classify if f + "_TRUTH.txt" in gt_files]


    for gt_filename in gt_files:
        # Skip temporary files
        if gt_filename.startswith("."):
            continue
        
        dataset = gt_filename.replace("_TRUTH.txt", "")
        
        # Skip if already classified and runtime recorded
        if dataset in runtimes:
            print(f"Skipping {dataset} - already classified")
            continue
        
        for fasta_filename in sorted(os.listdir(FASTA_DIR)):
            # If method_suffix is specified, only process files with that suffix
            if args.method_suffix:
                if not fasta_filename.endswith(f".{args.method_suffix}.gz"):
                    continue
                 
            fasta_fp = os.path.join(FASTA_DIR, fasta_filename)
            if fasta_filename.startswith(dataset) and os.path.isfile(fasta_fp):
                if "R2" in fasta_filename and os.path.isfile(fasta_fp.replace("R2", "R1")):
                    continue # avoid double evaluating the paired read files
                if "R1" in fasta_filename and os.path.isfile(fasta_fp.replace("R1", "R2")):
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
                    if args.preload_size:
                        cmd.extend(["--preload-size", args.preload_size])
                    cmd.extend([
                        fasta_fp,
                        fasta_fp.replace("R1", "R2")
                    ])
                    print(" ".join(cmd))
                    start_time = time.time()
                    subprocess.run(cmd)
                    runtime = time.time() - start_time
                    
                    # Save runtime
                    runtimes[dataset] = runtime
                    save_runtimes(runtime_file, runtimes)
                    print(f"Runtime for {dataset}: {runtime:.0f}s")
                else:
                    report_fp = os.path.join(OUTPUT_DIR, dataset + "_REPORT.tsv")
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
                    if args.preload_size:
                        cmd.extend(["--preload-size", args.preload_size])
                    cmd.append(fasta_fp)
                    print(" ".join(cmd))
                    start_time = time.time()
                    subprocess.run(cmd)
                    runtime = time.time() - start_time
                    
                    # Save runtime
                    runtimes[dataset] = runtime
                    save_runtimes(runtime_file, runtimes)
                    print(f"Runtime for {dataset}: {runtime:.0f}s")
                break

     
if __name__ == '__main__':
    main()
    