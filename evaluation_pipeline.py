#!/opt/conda/bin/python3
import argparse
import os
import re
import subprocess
import sys
from collections import Counter
import numpy as np
from sklearn.metrics import precision_recall_curve, auc, precision_score, recall_score, f1_score
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt


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
        if info['rank'] == rank and info['reads'] > 0:
            rank_counts[taxid] = info['reads']
    
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
    """Calculate precision, recall, F1 score, and AUPR."""
    all_taxids = set(pred_counts.keys()) | gt_taxids
    
    y_true = []
    y_pred = []
    scores = []
    
    for taxid in all_taxids:
        y_true.append(1 if taxid in gt_taxids else 0)
        y_pred.append(1 if taxid in pred_counts else 0)
        score = pred_counts.get(taxid, 0) / total_reads if total_reads > 0 else 0
        scores.append(score)
    
    y_true = np.array(y_true)
    y_pred = np.array(y_pred)
    scores = np.array(scores)
    
    if len(y_true) == 0 or np.sum(y_true) == 0:
        return {
            'precision': 0.0,
            'recall': 0.0,
            'f1': 0.0,
            'aupr': 0.0,
            'true_positives': 0,
            'false_positives': 0,
            'false_negatives': 0
        }
    
    precision = precision_score(y_true, y_pred, zero_division=0)
    recall = recall_score(y_true, y_pred, zero_division=0)
    f1 = f1_score(y_true, y_pred, zero_division=0)
    
    if np.sum(y_true) > 0 and len(np.unique(scores)) > 1:
        precisions, recalls, _ = precision_recall_curve(y_true, scores)
        aupr = auc(recalls, precisions)
    else:
        aupr = 0.0
    
    tp = np.sum((y_true == 1) & (y_pred == 1))
    fp = np.sum((y_true == 0) & (y_pred == 1))
    fn = np.sum((y_true == 1) & (y_pred == 0))
    
    return {
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'aupr': aupr,
        'true_positives': int(tp),
        'false_positives': int(fp),
        'false_negatives': int(fn)
    }

def eval_classification(report_fp, gt_genus_fp, gt_species_fp):
    """Evaluate KrakenUniq classification against ground truth and return metrics."""
    if not os.path.exists(report_fp):
        print(f"Warning: Report file not found: {report_fp}")
        return None
    
    # Parse KrakenUniq report
    taxonomy = parse_krakenuniq_report(report_fp)
    total_reads = taxonomy.get(1, {}).get('reads', 0)
    if total_reads == 0:
        total_reads = sum(info['reads'] for info in taxonomy.values() if info.get('reads', 0) > 0)
    
    # Process species level
    species_counts = extract_taxa_by_rank(taxonomy, 'species')
    gt_species, gt_species_ids = parse_ground_truth(gt_species_fp)
    species_metrics = calc_metrics(species_counts, gt_species, gt_species_ids, total_reads)
    
    # Process genus level
    genus_counts = extract_taxa_by_rank(taxonomy, 'genus')
    gt_genus, gt_genus_ids = parse_ground_truth(gt_genus_fp)
    genus_metrics = calc_metrics(genus_counts, gt_genus, gt_genus_ids, total_reads)
    
    return {
        'species': species_metrics,
        'genus': genus_metrics,
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
    file_handle.write(f"{result['species']['precision']:.4f}\t{result['species']['recall']:.4f}\t")
    file_handle.write(f"{result['species']['f1']:.4f}\t{result['species']['aupr']:.4f}\t")
    file_handle.write(f"{result['genus']['precision']:.4f}\t{result['genus']['recall']:.4f}\t")
    file_handle.write(f"{result['genus']['f1']:.4f}\t{result['genus']['aupr']:.4f}\t")
    file_handle.write(f"{result['total_reads']}\n")
    file_handle.flush()
    
    # Print metrics
    print(f"  Genus   - P: {result['genus']['precision']:.3f}, R: {result['genus']['recall']:.3f}, F1: {result['genus']['f1']:.3f}, AUPR: {result['genus']['aupr']:.3f}")
    print(f"  Species - P: {result['species']['precision']:.3f}, R: {result['species']['recall']:.3f}, F1: {result['species']['f1']:.3f}, AUPR: {result['species']['aupr']:.3f}")

def plot_metrics(all_results, plot_output):
    """Create and save evaluation metric plots."""
    if not all_results:
        return
    
    # Collect plot data
    datasets = [r['dataset'] for r in all_results]
    species_precision = [r['species']['precision'] for r in all_results]
    species_recall = [r['species']['recall'] for r in all_results]
    species_f1 = [r['species']['f1'] for r in all_results]
    species_aupr = [r['species']['aupr'] for r in all_results]
    genus_precision = [r['genus']['precision'] for r in all_results]
    genus_recall = [r['genus']['recall'] for r in all_results]
    genus_f1 = [r['genus']['f1'] for r in all_results]
    genus_aupr = [r['genus']['aupr'] for r in all_results]
    
    # Create plots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('KrakenUniq Classification Metrics', fontsize=16, fontweight='bold')
    
    x_pos = np.arange(len(datasets))
    
    # Precision subplot
    axes[0, 0].bar(x_pos - 0.2, species_precision, 0.4, label='Species', alpha=0.8)
    axes[0, 0].bar(x_pos + 0.2, genus_precision, 0.4, label='Genus', alpha=0.8)
    axes[0, 0].set_ylabel('Precision', fontweight='bold')
    axes[0, 0].set_title('Precision')
    axes[0, 0].set_xticks(x_pos)
    axes[0, 0].set_xticklabels(datasets, rotation=45, ha='right', fontsize=8)
    axes[0, 0].legend()
    axes[0, 0].set_ylim(0, 1)
    axes[0, 0].grid(axis='y', alpha=0.3)
    
    # Recall subplot
    axes[0, 1].bar(x_pos - 0.2, species_recall, 0.4, label='Species', alpha=0.8)
    axes[0, 1].bar(x_pos + 0.2, genus_recall, 0.4, label='Genus', alpha=0.8)
    axes[0, 1].set_ylabel('Recall', fontweight='bold')
    axes[0, 1].set_title('Recall')
    axes[0, 1].set_xticks(x_pos)
    axes[0, 1].set_xticklabels(datasets, rotation=45, ha='right', fontsize=8)
    axes[0, 1].legend()
    axes[0, 1].set_ylim(0, 1)
    axes[0, 1].grid(axis='y', alpha=0.3)
    
    # F1 Score subplot
    axes[1, 0].bar(x_pos - 0.2, species_f1, 0.4, label='Species', alpha=0.8)
    axes[1, 0].bar(x_pos + 0.2, genus_f1, 0.4, label='Genus', alpha=0.8)
    axes[1, 0].set_ylabel('F1 Score', fontweight='bold')
    axes[1, 0].set_title('F1 Score')
    axes[1, 0].set_xticks(x_pos)
    axes[1, 0].set_xticklabels(datasets, rotation=45, ha='right', fontsize=8)
    axes[1, 0].legend()
    axes[1, 0].set_ylim(0, 1)
    axes[1, 0].grid(axis='y', alpha=0.3)
    
    # AUPR subplot
    axes[1, 1].bar(x_pos - 0.2, species_aupr, 0.4, label='Species', alpha=0.8)
    axes[1, 1].bar(x_pos + 0.2, genus_aupr, 0.4, label='Genus', alpha=0.8)
    axes[1, 1].set_ylabel('AUPR', fontweight='bold')
    axes[1, 1].set_title('Area Under Precision-Recall Curve')
    axes[1, 1].set_xticks(x_pos)
    axes[1, 1].set_xticklabels(datasets, rotation=45, ha='right', fontsize=8)
    axes[1, 1].legend()
    axes[1, 1].set_ylim(0, 1)
    axes[1, 1].grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(plot_output, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"\nPlot saved to: {plot_output}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true", help="Run krakenuniq with --quick flag")
    args = parser.parse_args()

    DB_DIR = "/mnt/data/standard_db"
    FASTA_DIR = "/mnt/data/IMMSA"
    GT_DIR = os.path.join(FASTA_DIR, "truth_sets")
    GT_GENUS_DIR = os.path.join(GT_DIR, "genus")
    GT_SPECIES_DIR = os.path.join(GT_DIR, "species")
    OUTPUT_DIR = os.path.join(FASTA_DIR, "krakenuniq")
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    EVAL_OUTPUT = os.path.join(OUTPUT_DIR, "evaluation_summary.tsv")
    PLOT_OUTPUT = os.path.join(OUTPUT_DIR, "evaluation_metrics.png")

    if not os.path.isdir(FASTA_DIR):
        print(f"Error: {FASTA_DIR} is not a valid directory")
        return
    
    # Store all evaluation results
    all_results = []
    
    with open(EVAL_OUTPUT, 'w') as f:
        # Write header
        f.write("Dataset\tSpecies_Precision\tSpecies_Recall\tSpecies_F1\tSpecies_AUPR\t")
        f.write("Genus_Precision\tGenus_Recall\tGenus_F1\tGenus_AUPR\tTotal_Reads\n")

        for gt_filename in os.listdir(GT_GENUS_DIR):
            if gt_filename.startswith("."):
                continue
            for fasta_filename in os.listdir(FASTA_DIR):
                if fasta_filename.startswith(gt_filename.replace("_TRUTH.txt", "")):
                    fasta_fp = os.path.join(FASTA_DIR, fasta_filename)
                    if os.path.isfile(fasta_fp):
                        if "R1" in fasta_filename and os.path.isfile(fasta_fp.replace("R1", "R2")):
                            report_fp = os.path.join(OUTPUT_DIR, re.sub(r"_*R1|\.(fq|fast[aq])\.gz$", "", fasta_filename) + "_REPORT.tsv")
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
                            subprocess.run(cmd)
                        else:
                            report_fp = os.path.join(OUTPUT_DIR, re.sub(r"\.(fq|fast[aq])\.gz$", "", fasta_filename) + "_REPORT.tsv")
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
                            print(" ".join)
                            subprocess.run(cmd)
                        
                        # Evaluate classification
                        gt_genus_fp = os.path.join(GT_GENUS_DIR, gt_filename)
                        gt_species_fp = os.path.join(GT_SPECIES_DIR, gt_filename)
                        
                        if os.path.exists(gt_genus_fp) and os.path.exists(gt_species_fp):
                            print(f"Evaluating {gt_filename}...")
                            result = eva_classification(report_fp, gt_genus_fp, gt_species_fp)
                            if result:
                                result['dataset'] = gt_filename.replace("_TRUTH.txt", "")
                                all_results.append(result)
                                save_metrics(result, f)
    
        # Calculate averages and standard deviations
        avg_genus_p, std_genus_p, avg_species_p, std_species_p = calc_summary_measures("precision", all_results)
        avg_genus_r, std_genus_r, avg_species_r, std_species_r = calc_summary_measures("recall", all_results)
        avg_genus_f1, std_genus_f1, avg_species_f1, std_species_f1 = calc_summary_measures("f1", all_results)
        avg_genus_aupr, std_genus_aupr, avg_species_aupr, std_species_aupr = calc_summary_measures("aupr", all_results)
        avg_total_reads = np.mean([r['total_reads'] for r in all_results])
        std_total_reads = np.mean([r['total_reads'] for r in all_results])
        
        # Write averages and standard deviations
        f.write(f"AVERAGE\t")
        f.write(f"{avg_genus_p:.4f}\t{avg_genus_r:.4f}\t{avg_genus_f1:.4f}\t{avg_genus_aupr:.4f}\t")
        f.write(f"{avg_species_p:.4f}\t{avg_species_r:.4f}\t{avg_species_f1:.4f}\t{avg_species_aupr:.4f}\t")
        f.write(f"{avg_total_reads:.0f}\n")
        f.write(f"STD_DEV\t")
        f.write(f"{std_genus_p:.4f}\t{std_genus_r:.4f}\t{std_genus_f1:.4f}\t{std_genus_aupr:.4f}\t")
        f.write(f"{std_species_p:.4f}\t{std_species_r:.4f}\t{std_species_f1:.4f}\t{std_species_aupr:.4f}\t")
        f.write(f"{std_total_reads:.0f}\n")
    
    # Plot evaluation metrics and print summary
    plot_metrics(all_results, PLOT_OUTPUT)
    print(f"\n=== Evaluation Summary ===")
    print(f"Evaluated {len(all_results)} datasets")
    print(f"Results written to: {EVAL_OUTPUT}")
    print(f"\nAverage Metrics (±StdDev):")
    print(f"  Genus   - P: {avg_genus_p:.3f}±{std_genus_p:.3f}, R: {avg_genus_r:.3f}±{std_genus_r:.3f}, F1: {avg_genus_f1:.3f}±{std_genus_f1:.3f}, AUPR: {avg_genus_aupr:.3f}±{std_genus_aupr:.3f}")
    print(f"  Species - P: {avg_species_p:.3f}±{std_species_p:.3f}, R: {avg_species_r:.3f}±{std_species_r:.3f}, F1: {avg_species_f1:.3f}±{std_species_f1:.3f}, AUPR: {avg_species_aupr:.3f}±{std_species_aupr:.3f}")
    

if __name__ == '__main__':
    main()
    