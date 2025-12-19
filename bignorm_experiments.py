#!/opt/conda/bin/python3
"""
Bignorm Filtering Experiments
==============================

This script runs Bignorm normalization experiments with different parameter combinations
across various datasets (paired and unpaired reads, different file sizes).

Based on the Bignorm paper and Manual.txt, it:
1. Tests multiple parameter combinations for filtering options (-Q, -A, -B, -C, -N, -k)
2. Handles both paired-end and unpaired reads (fastq and fasta)
3. Parses Bignorm output statistics from stdout
4. Aggregates metrics into a single summary TSV file

References:
- Paper: "An Improved Filtering Algorithm for Big Read Datasets and its Application to Single-Cell Assembly"
- Manual: Bignorm/Manual.txt
- Decision function 6 (default): uses phred scores, requires fastq
- Decision function 3: simpler version without phred, works with fasta
"""

import os
import re
import subprocess
import time
import argparse
from pathlib import Path
from collections import defaultdict
import gzip

# PATHS
DATA_DIR = "/mnt/data"
DATASET_NAME = "IMMSA"
FASTA_DIR = os.path.join(DATA_DIR, DATASET_NAME)
BIGNORM_BIN = "/mgvd/Bignorm/Bignorm"
SUMMARY_FILE = None
LOGS_DIR = None


# Parameter combinations to test (based on paper and Manual.txt)
# Paper values: Q_0 (phred threshold), c_0 (rarity threshold A), B (contribution threshold), c_1 (abundance threshold C)
# From paper: bacteria use m=1024, t=10 (10 GB RAM); eukaryotes m=4096, t=10 (40 GB RAM)
PARAMETER_SETS = [
    # Default/baseline (from paper bacteria settings)
    {"k": 31, "Q": 20, "A": 3, "B": 5, "C": 20, "N": 10, "desc": "bignorm_default"},
    
    # Vary Q (minimum phred score threshold) - only for decision function 6 (fastq)
    {"k": 31, "Q": 10, "A": 3, "B": 5, "C": 20, "N": 10, "desc": "bignorm_low_phred"},
    {"k": 31, "Q": 30, "A": 3, "B": 5, "C": 20, "N": 10, "desc": "bignorm_high_phred"},
    
    # Vary A (rarity threshold c_0 in paper)
    {"k": 31, "Q": 20, "A": 2, "B": 5, "C": 20, "N": 10, "desc": "bignorm_low_rarity"},
    {"k": 31, "Q": 20, "A": 5, "B": 5, "C": 20, "N": 10, "desc": "bignorm_high_rarity"},
    
    # Vary B (contribution threshold)
    {"k": 31, "Q": 20, "A": 3, "B": 3, "C": 20, "N": 10, "desc": "bignorm_low_contribution"},
    {"k": 31, "Q": 20, "A": 3, "B": 10, "C": 20, "N": 10, "desc": "bignorm_high_contribution"},
    
    # Vary C (abundance threshold c_1 in paper)
    {"k": 31, "Q": 20, "A": 3, "B": 5, "C": 10, "N": 10, "desc": "bignorm_low_abundance"},
    {"k": 31, "Q": 20, "A": 3, "B": 5, "C": 50, "N": 10, "desc": "bignorm_high_abundance"},
    
    # Vary N (max number of N nucleotides before skipping read)
    {"k": 31, "Q": 20, "A": 3, "B": 5, "C": 20, "N": 5, "desc": "bignorm_strict_N"},
    {"k": 31, "Q": 20, "A": 3, "B": 5, "C": 20, "N": 20, "desc": "bignorm_relaxed_N"},
    
    # Conservative (high quality, strict filtering)
    {"k": 31, "Q": 30, "A": 5, "B": 10, "C": 50, "N": 5, "desc": "bignorm_conservative"},
    
    # Permissive (lower quality threshold, more reads kept)
    {"k": 31, "Q": 10, "A": 2, "B": 3, "C": 10, "N": 20, "desc": "bignorm_permissive"},
]

# For fasta files (decision function 3), we don't use Q parameter
PARAMETER_SETS_FASTA = [
    {"k": 31, "A": 3, "B": 5, "C": 20, "N": 10, "desc": "bignorm_default"},
    {"k": 31, "A": 2, "B": 5, "C": 20, "N": 10, "desc": "bignorm_low_rarity_fasta"},
    {"k": 31, "A": 5, "B": 5, "C": 20, "N": 10, "desc": "bignorm_high_rarity_fasta"},
    {"k": 31, "A": 3, "B": 3, "C": 20, "N": 10, "desc": "bignorm_ow_contribution_fasta"},
    {"k": 31, "A": 3, "B": 10, "C": 20, "N": 10, "desc": "bignorm_high_contribution_fasta"},
    {"k": 31, "A": 3, "B": 5, "C": 10, "N": 10, "desc": "bignorm_low_abundance_fasta"},
    {"k": 31, "A": 3, "B": 5, "C": 50, "N": 10, "desc": "bignorm_high_abundance_fasta"},
]


def get_file_size(filepath):
    """Get file size in bytes."""
    if os.path.exists(filepath):
        return os.path.getsize(filepath)
    return 0


def count_reads_in_file(filepath):
    """Count total reads in a fastq/fasta file (gzipped or not)."""
    try:
        if filepath.endswith('.gz'):
            opener = gzip.open
        else:
            opener = open
        
        count = 0
        with opener(filepath, 'rt') as f:
            # Determine format from first character
            first_char = f.read(1)
            f.seek(0)
            
            if first_char == '@':  # FASTQ
                for i, line in enumerate(f):
                    if i % 4 == 0:
                        count += 1
            elif first_char == '>':  # FASTA
                for line in f:
                    if line.startswith('>'):
                        count += 1
        return count
    except Exception as e:
        print(f"Warning: Could not count reads in {filepath}: {e}")
        return 0


def is_fastq(filepath):
    """Check if file is FASTQ format (vs FASTA)."""
    # Check by extension first
    if filepath.endswith('.fq.gz') or filepath.endswith('.fq') or filepath.endswith('.fastq.gz') or filepath.endswith('.fastq'):
        return True
    elif filepath.endswith('.fa.gz') or filepath.endswith('.fa') or filepath.endswith('.fasta.gz') or filepath.endswith('.fasta'):
        return False
    
    # Check by content (first character)
    try:
        if filepath.endswith('.gz'):
            with gzip.open(filepath, 'rt') as f:
                first_char = f.read(1)
        else:
            with open(filepath, 'r') as f:
                first_char = f.read(1)
        return first_char == '@'
    except:
        return True  # Default to fastq


def parse_bignorm_output(output_text):
    """
    Parse Bignorm output statistics from stdout.
    
    Example output format (from bignorm_output_example.txt):
        Reads processed:
            unpaired: 0
            paired  : 4162035
            ----------
            total   : 4162035
        
        Reads accepted:
            unpaired: 0 
            paired  : 2676627
            ----------
            total   : 2676627
        
        fp:     0.000 (base: 0.344)
        reads kept: 64.31 %
        time needed: 438488 s
    """
    stats = {
        'reads_processed': 0,
        'reads_accepted': 0,
        'reads_kept_pct': 0.0,
        'fp': 0.0,
        'time_s': 0,
    }
    
    # Parse reads processed (total)
    match = re.search(r'total\s*:\s*(\d+)', output_text)
    if match:
        stats['reads_processed'] = int(match.group(1))
    
    # Parse reads accepted (look for second occurrence after "Reads accepted:")
    accepted_section = re.search(r'Reads accepted:.*?total\s*:\s*(\d+)', output_text, re.DOTALL)
    if accepted_section:
        stats['reads_accepted'] = int(accepted_section.group(1))
    
    # Parse fp (false positive rate)
    fp_match = re.search(r'fp:\s*([\d\.]+)', output_text)
    if fp_match:
        stats['fp'] = float(fp_match.group(1))
    
    # Parse percentage of reads kept
    pct_match = re.search(r'reads kept:\s*([\d\.]+)\s*%', output_text)
    if pct_match:
        stats['reads_kept_pct'] = float(pct_match.group(1))
    
    # Parse time
    time_match = re.search(r'time needed:\s*(\d+)\s*s', output_text)
    if time_match:
        stats['time_s'] = int(time_match.group(1))
    
    return stats


def run_bignorm_experiment(dataset, params, paired_file=None):
    """
    Run a single Bignorm experiment with given parameters.
    
    Args:
        read_file: Path to the input read file (unpaired or first of pair)
        params: Dictionary of parameters (k, Q, A, B, C, N)
        paired_file: Path to paired read file if paired-end data
    
    Returns:
        Dictionary with experiment results and statistics
    """
    
    is_fq = is_fastq(read_file)
    
    # Determine decision function (6 for fastq, 3 for fasta)
    decision_func = 6 if is_fq else 3
    
    # Build Bignorm command
    cmd = [BIGNORM_BIN]
    
    # Input files
    if paired_file:
        cmd.extend(['-1', read_file, '-2', paired_file])
    else:
        cmd.extend(['-u', read_file])
    
    # Output compression
    cmd.append('-z')
    
    # Decision function
    cmd.extend(['-d', str(decision_func)])
    
    # Parameters
    cmd.extend(['-k', str(params['k'])])
    cmd.extend(['-A', str(params['A'])])
    cmd.extend(['-B', str(params['B'])])
    cmd.extend(['-C', str(params['C'])])
    cmd.extend(['-N', str(params['N'])])
    
    # Q parameter only for decision function 6 (fastq)
    if decision_func == 6 and 'Q' in params:
        cmd.extend(['-Q', str(params['Q'])])
    
    # Set memory and hash lines (use default auto-allocation)
    # For bacteria: -m 1024 -t 10 (10 GB)
    # Let Bignorm auto-detect for now
    
    # Log file for this experiment
    log_name = f"{base_name}_{params['desc']}_log.txt"
    log_path = os.path.join(LOGS_DIR, log_name)
    
    print(f"\n{'='*80}")
    print(f"Running experiment: {params['desc']}")
    print(f"  Input: {base_name}")
    if paired_file:
        print(f"  Paired: {os.path.basename(paired_file)}")
    print(f"  Parameters: k={params['k']}, A={params['A']}, B={params['B']}, C={params['C']}, N={params['N']}", end='')
    if 'Q' in params:
        print(f", Q={params['Q']}", end='')
    print(f"\n  Decision function: {decision_func}")
    print(f"  Command: {' '.join(cmd)}")
    print(f"  Log: {log_path}")
    
    # Run Bignorm and capture output
    start_time = time.time()
    try:
        result = subprocess.run(
            cmd,
            cwd=os.path.dirname(read_file),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        )
        elapsed_time = time.time() - start_time
        
        # Save log
        with open(log_path, 'w') as f:
            f.write(f"Command: {' '.join(cmd)}\n")
            f.write(f"Return code: {result.returncode}\n")
            f.write(f"Elapsed time: {elapsed_time:.2f} s\n\n")
            f.write(result.stdout)
        
        # Parse statistics from output
        stats = parse_bignorm_output(result.stdout)
        
        # Get output file sizes and rename from "_keep.gz" to experiment name
        # With -z flag, Bignorm always produces _keep.gz files
        if paired_file:
            # Bignorm creates _keep.gz files
            out1_keep = read_file + '_keep.gz'
            out2_keep = paired_file + '_keep.gz'
            
            # Rename to experiment name
            out1 = read_file + f'.{params["desc"]}.gz'
            out2 = paired_file + f'.{params["desc"]}.gz'
            
            if os.path.exists(out1_keep):
                os.rename(out1_keep, out1)
            else:
                print(f"{out1_keep} was not found and therefore not renamed")
            if os.path.exists(out2_keep):
                os.rename(out2_keep, out2)
            else:
                print(f"{out2_keep} was not found and therefore not renamed")
            
            output_size = get_file_size(out1) + get_file_size(out2)
            output_files = f"{os.path.basename(out1)}, {os.path.basename(out2)}"
        else:
            # Bignorm creates _keep.gz files
            out1_keep = read_file + '_keep.gz'
            
            # Rename to experiment name
            out1 = read_file + f'.{params["desc"]}.gz'
            
            if os.path.exists(out1_keep):
                os.rename(out1_keep, out1)
            else:
                print(f"{out1_keep} was not found and therefore not renamed")
            
            output_size = get_file_size(out1)
            output_files = os.path.basename(out1)
        
        input_size = get_file_size(read_file)
        if paired_file:
            input_size += get_file_size(paired_file)
        
        # Compile results
        experiment_result = {
            'read_file': base_name,
            'paired_file': os.path.basename(paired_file) if paired_file else 'N/A',
            'param_desc': params['desc'],
            'decision_func': decision_func,
            'k': params['k'],
            'Q': params.get('Q', 'N/A'),
            'A': params['A'],
            'B': params['B'],
            'C': params['C'],
            'N': params['N'],
            'input_size_bytes': input_size,
            'output_size_bytes': output_size,
            'output_files': output_files,
            'reads_processed': stats['reads_processed'],
            'reads_accepted': stats['reads_accepted'],
            'reads_kept_pct': stats['reads_kept_pct'],
            'fp': stats['fp'],
            'runtime_s': stats['time_s'] if stats['time_s'] > 0 else elapsed_time,
            'exit_code': result.returncode,
            'log_file': log_name,
        }
        
        print(f"  ✓ Completed in {elapsed_time:.1f}s")
        print(f"    Reads: {stats['reads_processed']:,} → {stats['reads_accepted']:,} ({stats['reads_kept_pct']:.2f}%)")
        print(f"    Size: {input_size / 1024**2:.1f} MB → {output_size / 1024**2:.1f} MB")
        print(f"    FP rate: {stats['fp']:.3f}")
        
        return experiment_result
        
    except subprocess.TimeoutExpired:
        elapsed_time = time.time() - start_time
        print(f"  ✗ TIMEOUT after {elapsed_time:.1f}s")
        return {
            'read_file': base_name,
            'paired_file': os.path.basename(paired_file) if paired_file else 'N/A',
            'param_desc': params['desc'],
            'decision_func': decision_func,
            'k': params['k'],
            'Q': params.get('Q', 'N/A'),
            'A': params['A'],
            'B': params['B'],
            'C': params['C'],
            'N': params['N'],
            'input_size_bytes': 0,
            'output_size_bytes': 0,
            'output_files': 'TIMEOUT',
            'reads_processed': 0,
            'reads_accepted': 0,
            'reads_kept_pct': 0.0,
            'fp': 0.0,
            'runtime_s': elapsed_time,
            'exit_code': -1,
            'log_file': log_name,
        }
    except Exception as e:
        print(f"  ✗ ERROR: {e}")
        return None


def write_summary_header(f):
    """Write TSV header for summary file."""
    headers = [
        'dataset',
        'param_set',
        'reads_processed',
        'reads_accepted',
        'reads_kept_pct',
        'fp',
        'runtime_s',
    ]
    f.write('\t'.join(headers) + '\n')


def get_existing_results(summary_file):
    """Load existing results from summary file."""
    results = {}
    if not os.path.exists(summary_file):
        return results
    
    with open(summary_file, 'r') as f:
        lines = f.readlines()
        if len(lines) <= 1:  # Only header or empty
            return results
        
        for line in lines[1:]:  # Skip header
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                dataset = parts[0]
                param_set = parts[1]
                key = (dataset, param_set)
                results[key] = line  # Store the full line
    
    return results


def write_summary_row(existing_results, result):
    """Add or update a result row in existing results dict."""
    if result is None or result['exit_code'] != 0:
        return
    
    # Extract dataset name from read_file (remove extensions and R1/R2)
    dataset = result['read_file']
    # Remove common extensions
    dataset = re.sub(r'\.(fq|fastq|fa|fasta)(\.gz)?$', '', dataset)
    # Remove R1/R2 markers
    dataset = re.sub(r'_R[12]$', '', dataset)
    
    param_set = result['param_desc']
    key = (dataset, param_set)
    
    # Format the row
    row = [
        dataset,
        param_set,
        str(result['reads_processed']),
        str(result['reads_accepted']),
        f"{result['reads_kept_pct']:.2f}",
        f"{result['fp']:.3f}",
        f"{result['runtime_s']:.0f}",
    ]
    
    existing_results[key] = '\t'.join(row) + '\n'


def main():
    global LOGS_DIR, SUMMARY_FILE
    parser = argparse.ArgumentParser(
        description='Run Bignorm normalization experiments with various parameter combinations.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
            Examples:
            # Run all experiments on all files in IMMSA_bignorm_experiments
            python3 bignorm_experiments.py

            # Run only on a specific file
            python3 bignorm_experiments.py --file sample_1.fq.gz

            # Run only baseline parameters
            python3 bignorm_experiments.py --baseline-only
        """
    )
    parser.add_argument('--file', type=str, help='Run only on a specific input file')
    parser.add_argument('--baseline-only', action='store_true', help='Run only baseline parameter set')
    parser.add_argument('--experiment-dir', type=str, help='Specify experiment directory (logs)')
    args = parser.parse_args()
    
    EXPERIMENT_DIR = args.experiment_dir
    SUMMARY_FILE = os.path.join(EXPERIMENT_DIR, "summary.tsv")
    LOGS_DIR = os.path.join(EXPERIMENT_DIR, "logs")
    # Create output directories
    os.makedirs(EXPERIMENT_DIR, exist_ok=True)
    os.makedirs(LOGS_DIR, exist_ok=True)

    # Check if Bignorm binary exists
    if not os.path.exists(BIGNORM_BIN):
        print(f"Error: Bignorm binary not found at {BIGNORM_BIN}")
        print("Please compile Bignorm first: cd Bignorm && make")
        return 1
    
    # Check if fasta directory exists
    if not os.path.isdir(FASTA_DIR):
        print(f"Error: Fasta directory not found: {FASTA_DIR}")
        print("Please ensure /mnt/data/IMMSA_bignorm_experiments exists with input files")
        return 1
    
    print(f"\n{'='*80}")
    print(f"Bignorm Normalization Experiments")
    print(f"{'='*80}")
    print(f"Input directory: {FASTA_DIR}")
    print(f"Output directory: {EXPERIMENT_DIR}")
    
    if args.baseline_only:
        print("\nRunning BASELINE ONLY experiments")
    else:
        print(f"\nRunning {len(PARAMETER_SETS)} parameter combinations per FASTQ file")
        print(f"Running {len(PARAMETER_SETS_FASTA)} parameter combinations per FASTA file")
    
    # Load existing results
    existing_results = get_existing_results(SUMMARY_FILE)
    print(f"Loaded {len(existing_results)} existing results")
    
    total_experiments = 0
    successful_experiments = 0
    
    # Process files directly
    for filename in sorted(os.listdir(FASTA_DIR)):
        filepath = os.path.join(FASTA_DIR, filename)
        
        # Skip non-sequence files
        if not (filename.endswith('.fq.gz') or filename.endswith('.fq') or 
                filename.endswith('.fastq.gz') or filename.endswith('.fastq') or
                filename.endswith('.fa.gz') or filename.endswith('.fa') or
                filename.endswith('.fasta.gz') or filename.endswith('.fasta')):
            continue
        
        # Skip if specific file requested and this isn't it
        if args.file and filepath != args.file:
            continue
        
        # Skip already-processed output files or R2 files
        if '_keep' in filename or 'R2' in filename:
            continue
        
        # Determine if paired
        paired_file = None
        if 'R1' in filename:
            r2_file = filepath.replace('R1', 'R2')
            if os.path.isfile(r2_file):
                paired_file = r2_file
        
        # Select parameter sets based on file format
        is_fq = is_fastq(filepath)
        param_sets = PARAMETER_SETS if is_fq else PARAMETER_SETS_FASTA
        
        if args.baseline_only:
            param_sets = [p for p in param_sets if p['desc'] == 'bignorm_default']
        
        # Run experiments for this file
        for params in param_sets:
            total_experiments += 1
            result = run_bignorm_experiment(filepath, params, paired_file)
            
            if result and result['exit_code'] == 0:
                successful_experiments += 1
                write_summary_row(existing_results, result)
    
    # Write all results to file (sorted by dataset and param_set)
    with open(SUMMARY_FILE, 'w') as f:
        write_summary_header(f)
        for key in sorted(existing_results.keys()):
            f.write(existing_results[key])
    
    print(f"\n{'='*80}")
    print(f"Experiments Complete")
    print(f"{'='*80}")
    print(f"Total experiments: {total_experiments}")
    print(f"Successful: {successful_experiments}")
    print(f"Failed: {total_experiments - successful_experiments}")
    print(f"\nResults saved to: {SUMMARY_FILE}")
    print(f"Logs saved to: {LOGS_DIR}")
    
    return 0


if __name__ == '__main__':
    exit(main())
