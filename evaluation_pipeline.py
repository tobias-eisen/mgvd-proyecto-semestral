#!/opt/conda/bin/python3
import argparse
import os
import re
import subprocess
import sys
import time
from collections import Counter
import numpy as np
from sklearn.metrics import precision_recall_curve, auc, precision_score, recall_score, f1_score
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt

# PATHS
DATA_DIR = "/mnt/data/"
DATASET_NAME = "IMMSA"
FASTA_DIR = os.path.join(DATA_DIR, DATASET_NAME)
DB_DIR = os.path.join(DATA_DIR, "standard_db")
GT_DIR = os.path.join(FASTA_DIR, "truth_sets")
GT_GENUS_DIR = os.path.join(GT_DIR, "genus")
GT_SPECIES_DIR = os.path.join(GT_DIR, "species")


def parse_krakenuniq_report(report_file):
    """Parse KrakenUniq REPORTFILE.tsv to get taxonomy information."""
    taxonomy = {}
    
    with open(report_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('%'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            reads = int(parts[1])
            taxReads = int(parts[2])
            taxid = int(parts[6])
            rank = parts[7]
            name = parts[8].strip()
            
            taxonomy[taxid] = {
                'rank': rank,
                'name': name,
                'reads': reads,
                'taxReads': taxReads
            }
    
    return taxonomy

def extract_taxa_by_rank(taxonomy, rank):
    """Extract all taxa at a specific rank from the taxonomy report."""
    rank_counts = Counter()
    
    for taxid, info in taxonomy.items():
        if info['rank'] == rank and info['taxReads'] > 0:
            rank_counts[taxid] = info['taxReads']
    
    return rank_counts

def parse_ground_truth(gt_file):
    """Parse ground truth file."""
    gt_taxa = {}
    gt_taxids = set()
    
    with open(gt_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            
            taxid = int(parts[0])
            fraction = float(parts[1])
            abundance = float(parts[2])
            rank = parts[3]
            name = parts[4]
            
            gt_taxa[taxid] = {
                'fraction': fraction,
                'abundance': abundance,
                'rank': rank,
                'name': name
            }
            gt_taxids.add(taxid)
    
    return gt_taxa, gt_taxids

def calc_metrics(pred_counts, gt_taxa, gt_taxids, total_reads):
    """Calculate precision, recall, F1 score, and AUPR using threshold-based approach."""
    pred_taxids = set(pred_counts.keys())
    
    if len(pred_taxids) == 0:
        return {
            'precision': 0.0,
            'recall': 0.0,
            'f1': 0.0,
            'aupr': 0.0,
            'true_positives': 0,
            'false_positives': 0,
            'false_negatives': 0
        }
    
    # Sort predictions by abundance (descending)
    sorted_preds = sorted(pred_counts.items(), key=lambda x: x[1], reverse=True)
    
    # Calculate metrics at each threshold (like in the manuscript)
    total_positives = len(gt_taxids)
    f1_scores = []
    recalls = []
    precisions = []
    
    for j in range(1, len(sorted_preds) + 1):
        # Consider top j predictions
        top_j_taxids = set([taxid for taxid, _ in sorted_preds[:j]])
        tp = len(top_j_taxids & gt_taxids)
        
        # Precision = TP / (TP + FP) = TP / j
        precision = tp / j
        # Recall = TP / total_positives
        recall = tp / total_positives if total_positives > 0 else 0
        
        precisions.append(precision)
        recalls.append(recall)
        
        # F1 score
        if precision + recall > 0:
            f1 = 2 * precision * recall / (precision + recall)
        else:
            f1 = 0.0
        f1_scores.append(f1)
    
    # Get maximum F1 score
    max_f1_idx = np.argmax(f1_scores)
    max_f1 = f1_scores[max_f1_idx]
    max_precision = precisions[max_f1_idx]
    max_recall = recalls[max_f1_idx]
    
    # Calculate final TP, FP, FN at max F1 threshold
    top_taxa = set([taxid for taxid, _ in sorted_preds[:max_f1_idx + 1]])
    tp = len(top_taxa & gt_taxids)
    fp = len(top_taxa - gt_taxids)
    fn = len(gt_taxids - top_taxa)
    
    # For AUPR, use all predictions with their scores
    all_taxids = pred_taxids | gt_taxids
    y_true = []
    scores = []
    
    for taxid in all_taxids:
        y_true.append(1 if taxid in gt_taxids else 0)
        score = pred_counts.get(taxid, 0) / total_reads if total_reads > 0 else 0
        scores.append(score)
    
    y_true = np.array(y_true)
    scores = np.array(scores)
    
    # Calculate AUPR if there are positive examples
    if np.sum(y_true) > 0:
        try:
            precisions_pr, recalls_pr, _ = precision_recall_curve(y_true, scores)
            aupr = auc(recalls_pr, precisions_pr)
        except:
            aupr = 0.0
    else:
        aupr = 0.0
    
    return {
        'precision': max_precision,
        'recall': max_recall,
        'f1': max_f1,
        'aupr': aupr,
        'true_positives': int(tp),
        'false_positives': int(fp),
        'false_negatives': int(fn)
    }

def eval_classification(dataset, eval_dir):
    """Evaluate KrakenUniq classification against ground truth and return metrics."""
    report_fp = os.path.join(eval_dir, dataset+"_REPORT.tsv")
    gt_genus_fp = os.path.join(GT_GENUS_DIR, dataset+"_TRUTH.txt")
    gt_species_fp = os.path.join(GT_SPECIES_DIR, dataset+"_TRUTH.txt")

    if not os.path.exists(report_fp):
        print(f"Warning: Report file not found: {report_fp}")
        return None
    
    # Parse KrakenUniq report
    taxonomy = parse_krakenuniq_report(report_fp)
    total_reads = taxonomy.get(1, {}).get('reads', 0)
    if total_reads == 0:
        total_reads = sum(info['reads'] for info in taxonomy.values() if info.get('reads', 0) > 0)
        
    # Process genus level
    genus_counts = extract_taxa_by_rank(taxonomy, 'genus')
    gt_genus, gt_genus_ids = parse_ground_truth(gt_genus_fp)
    genus_metrics = calc_metrics(genus_counts, gt_genus, gt_genus_ids, total_reads)
    
    # Process species level
    species_counts = extract_taxa_by_rank(taxonomy, 'species')
    gt_species, gt_species_ids = parse_ground_truth(gt_species_fp)
    species_metrics = calc_metrics(species_counts, gt_species, gt_species_ids, total_reads)

    return {
        'genus': genus_metrics,
        'species': species_metrics,
        'total_reads': total_reads
    }

def calc_summary_measures(metric, all_results):
    avg_genus = np.mean([r['genus'][metric] for r in all_results])
    std_genus = np.std([r['genus'][metric] for r in all_results])
    avg_species = np.mean([r['species'][metric] for r in all_results])
    std_species = np.std([r['species'][metric] for r in all_results])
    return avg_genus, std_genus, avg_species, std_species

def save_metrics(result, file_handle):
    """Write evaluation result to output file."""
    file_handle.write(f"{result['dataset']}\t")
    file_handle.write(f"{result['genus']['precision']:.4f}\t{result['genus']['recall']:.4f}\t")
    file_handle.write(f"{result['genus']['f1']:.4f}\t{result['genus']['aupr']:.4f}\t")
    file_handle.write(f"{result['species']['precision']:.4f}\t{result['species']['recall']:.4f}\t")
    file_handle.write(f"{result['species']['f1']:.4f}\t{result['species']['aupr']:.4f}\t")
    file_handle.write(f"{result['total_reads']}\t{result['runtime']}\n")
    file_handle.flush()
    
    # Print metrics
    print(f"  Genus   - P: {result['genus']['precision']:.3f}, R: {result['genus']['recall']:.3f}, F1: {result['genus']['f1']:.3f}, AUPR: {result['genus']['aupr']:.3f}")
    print(f"  Species - P: {result['species']['precision']:.3f}, R: {result['species']['recall']:.3f}, F1: {result['species']['f1']:.3f}, AUPR: {result['species']['aupr']:.3f}")
    print(f"  Runtime: {result['runtime']}s")

def get_already_evaluated(eval_path):
    if not os.path.exists(eval_path):
        return set()
    
    print(f"Found existing evaluation file: {eval_path}")
    already_evaluated = set()
    all_results = []
    with open(eval_path, 'r') as f:
        lines = f.readlines()
        if lines:
            for line in lines[1:]:  # Skip header
                if line.startswith("AVERAGE") or line.startswith("STD_DEV"):
                    continue
                parts = line.strip().split('\t')
                dataset_name = parts[0] + "_TRUTH.txt"
                already_evaluated.add(dataset_name)
                result = {
                    'dataset': parts[0],
                    'genus': {
                        'precision': float(parts[1]),
                        'recall': float(parts[2]),
                        'f1': float(parts[3]),
                        'aupr': float(parts[4])
                    },
                    'species': {
                        'precision': float(parts[5]),
                        'recall': float(parts[6]),
                        'f1': float(parts[7]),
                        'aupr': float(parts[8])
                    },
                    'total_reads': int(float(parts[9])),
                    'runtime': float(parts[10])
                }
                all_results.append(result)
    return already_evaluated, all_results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--method-suffix", type=str, required=False, help="Method suffix for compressed read files")
    parser.add_argument("--quick", action="store_true", help="Run krakenuniq with --quick flag")
    parser.add_argument("--files-to-evaluate", nargs='+', help="List of ground truth files to evaluate")
    args = parser.parse_args()
    
    # Build output directory name based on arguments
    output_dir_name = "krakenuniq"
    if args.method_suffix:
        output_dir_name += f"_{args.method_suffix}"
    if args.quick:
        output_dir_name += "_quick"
    
    OUTPUT_DIR = os.path.join(DATASET_NAME+"_evaluation", output_dir_name)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    EVAL_OUTPUT = os.path.join(OUTPUT_DIR, "evaluation_summary.tsv")
    PLOT_OUTPUT = os.path.join(OUTPUT_DIR, "evaluation_metrics.png")

    if not os.path.isdir(FASTA_DIR):
        print(f"Error: {FASTA_DIR} is not a valid directory")
        return
    
    # If --files-to-evaluate is given, only evaluate those files
    gt_files = sorted(os.listdir(GT_GENUS_DIR))
    if args.files_to_evaluate:
        gt_files = [f + "_TRUTH.txt" for f in args.files_to_evaluate if f + "_TRUTH.txt" in gt_files]

    # Check for existing evaluations
    already_evaluated, all_results = get_already_evaluated(EVAL_OUTPUT) if os.path.isfile(EVAL_OUTPUT) else set(), []

    # Only write header if no evaluations exist
    if len(already_evaluated) == 0:
        with open(EVAL_OUTPUT, 'w') as f:
            f.write("Dataset\tGenus_Precision\tGenus_Recall\tGenus_F1\tGenus_AUPR\t")
            f.write("Species_Precision\tSpecies_Recall\tSpecies_F1\tSpecies_AUPR\tTotal_Reads\tRuntime_sec\n")

    with open(EVAL_OUTPUT, 'a') as f:
        for gt_filename in gt_files:
            # Skip temporary files
            if gt_filename.startswith("."):
                continue
            
            # Skip if already evaluated
            if gt_filename in already_evaluated:
                print(f"Skipping already evaluated: {gt_filename}")
                continue
            
            for fasta_filename in sorted(os.listdir(FASTA_DIR)):
                # If method_suffix is specified, only process files with that suffix
                if args.method_suffix:
                    if not fasta_filename.endswith(f".{args.method_suffix}.gz"):
                        continue
                
                if fasta_filename.startswith(gt_filename.replace("_TRUTH.txt", "")):
                    fasta_fp = os.path.join(FASTA_DIR, fasta_filename)
                    if os.path.isfile(fasta_fp):
                        if "R2" in fasta_filename and os.path.isfile(fasta_fp.replace("R2", "R1")):
                            continue # avoid double evaluating the paired read files
                        if "R1" in fasta_filename and os.path.isfile(fasta_fp.replace("R1", "R2")):
                            dataset = gt_filename.replace("_TRUTH.txt", "")
                            report_fp = os.path.join(OUTPUT_DIR, dataset+"_REPORT.tsv")
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
                            report_fp = os.path.join(OUTPUT_DIR, dataset + "_REPORT.tsv")
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
                        
                        # Evaluate classification
                        print(f"Evaluating {gt_filename}...")
                        result = eval_classification(dataset, eval_dir=OUTPUT_DIR)
                        if result:
                            result['dataset'] = gt_filename.replace("_TRUTH.txt", "")
                            result['runtime'] = round(runtime)
                            all_results.append(result)
                            save_metrics(result, f)
    
if __name__ == '__main__':
    main()
    