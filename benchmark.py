#!/opt/conda/bin/python3
import argparse
from collections import defaultdict, Counter
import numpy as np
import os
from sklearn.metrics import precision_recall_curve, auc, precision_score, recall_score, f1_score
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# PATHS - will be set in main()
GT_GENUS_DIR = None
GT_SPECIES_DIR = None


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

def convert_krakenuniq_to_ground_truth_format(tool_dir, dataset, tool_name):
    """Convert KrakenUniq report to ground truth format files in genus/ and species/ subdirs."""
    # Try both naming conventions for KrakenUniq reports
    report_fp = os.path.join(tool_dir, dataset + "_REPORT.tsv")
    if not os.path.exists(report_fp):
        for filename in os.listdir(tool_dir):
            if filename.startswith(dataset) and filename.endswith('.krakenuniq.report.tsv'):
                report_fp = os.path.join(tool_dir, filename)
                break
    
    if not os.path.exists(report_fp):
        print(f"  Warning: Report file not found for {dataset} in {tool_dir}")
        return False
    
    # Parse KrakenUniq report
    taxonomy = parse_krakenuniq_report(report_fp)
    total_reads = taxonomy.get(1, {}).get('reads', 0)
    if total_reads == 0:
        total_reads = sum(info.get('reads', 0) for info in taxonomy.values())
    
    if total_reads == 0:
        print(f"  Warning: No reads found in {report_fp}")
        return False
    
    # Create genus and species subdirectories
    genus_dir = os.path.join(tool_dir, 'genus')
    species_dir = os.path.join(tool_dir, 'species')
    os.makedirs(genus_dir, exist_ok=True)
    os.makedirs(species_dir, exist_ok=True)
    
    # Extract and write genus level
    genus_counts = extract_taxa_by_rank(taxonomy, 'genus')
    genus_file = os.path.join(genus_dir, f"{dataset}_{tool_name}.txt")
    with open(genus_file, 'w') as f:
        for taxid, reads in genus_counts.items():
            fraction = reads / total_reads if total_reads > 0 else 0
            abundance = reads
            rank = taxonomy[taxid]['rank']
            name = taxonomy[taxid]['name']
            f.write(f"{taxid}\t{fraction}\t{abundance}\t{rank}\t{name}\n")
    
    # Extract and write species level
    species_counts = extract_taxa_by_rank(taxonomy, 'species')
    species_file = os.path.join(species_dir, f"{dataset}_{tool_name}.txt")
    with open(species_file, 'w') as f:
        for taxid, reads in species_counts.items():
            fraction = reads / total_reads if total_reads > 0 else 0
            abundance = reads
            rank = taxonomy[taxid]['rank']
            name = taxonomy[taxid]['name']
            f.write(f"{taxid}\t{fraction}\t{abundance}\t{rank}\t{name}\n")
    
    return True


def plot_benchmark_results(all_tool_results, output_file):
    """Create a 4-panel barplot showing metrics with standard error bars."""
    tools = sorted(all_tool_results.keys())
    n_tools = len(tools)
    
    # Set up the figure with 2x2 subplots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Benchmark Results by Tool', fontsize=16, fontweight='bold')
    
    metrics = ['precision', 'recall', 'f1', 'aupr']
    metric_titles = ['Precision', 'Recall', 'F1 Score', 'AUPR']
    
    x = np.arange(n_tools)
    width = 0.35
    
    for idx, (metric, title) in enumerate(zip(metrics, metric_titles)):
        ax = axes[idx // 2, idx % 2]
        
        genus_means = [all_tool_results[tool][f'genus_{metric}_mean'] for tool in tools]
        genus_stds = [all_tool_results[tool][f'genus_{metric}_std'] for tool in tools]
        
        species_means = [all_tool_results[tool][f'species_{metric}_mean'] for tool in tools]
        species_stds = [all_tool_results[tool][f'species_{metric}_std'] for tool in tools]
        
        # Create bars
        bars1 = ax.bar(x - width/2, genus_means, width, label='Genus', 
                       alpha=0.8, color='steelblue')
        bars2 = ax.bar(x + width/2, species_means, width, label='Species',
                       alpha=0.8, color='coral')
        
        # Customize subplot
        ax.set_ylabel(title, fontsize=12, fontweight='bold')
        ax.set_title(f'{title} Comparison', fontsize=13, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(tools, rotation=45, ha='right', fontsize=10)
        
        # Color krakenuniq tool labels red
        for i, label in enumerate(ax.get_xticklabels()):
            if 'krakenuniq' in tools[i].lower():
                label.set_color('red')
        
        ax.legend(fontsize=10)
        ax.grid(axis='y', alpha=0.3, linestyle='--')
        ax.set_ylim(0, 1.0)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Benchmark plot saved to: {output_file}")


def eval_tool_classification(dataset, tool_dir, tool_name):
    """Evaluate a tool's classification against ground truth."""
    gt_genus_fp = os.path.join(GT_GENUS_DIR, dataset + "_TRUTH.txt")
    gt_species_fp = os.path.join(GT_SPECIES_DIR, dataset + "_TRUTH.txt")
    
    if not os.path.exists(gt_genus_fp) or not os.path.exists(gt_species_fp):
        print(f"  Warning: Ground truth not found for {dataset}")
        return None
    
    # Parse ground truth
    gt_genus, gt_genus_ids = parse_ground_truth(gt_genus_fp)
    gt_species, gt_species_ids = parse_ground_truth(gt_species_fp)
    
    # Find tool output files in genus/species subdirectories
    genus_file = os.path.join(tool_dir, 'genus', f"{dataset}_{tool_name}.txt")
    species_file = os.path.join(tool_dir, 'species', f"{dataset}_{tool_name}.txt")
    
    if not os.path.exists(genus_file) or not os.path.exists(species_file):
        print(f"  Warning: Tool output not found for {dataset}_{tool_name}")
        return None
    
    # Parse tool predictions
    pred_genus, pred_genus_ids = parse_ground_truth(genus_file)
    pred_species, pred_species_ids = parse_ground_truth(species_file)
    
    # Extract counts from predictions
    genus_counts = Counter()
    for taxid in pred_genus_ids:
        genus_counts[taxid] = int(pred_genus[taxid]['abundance'])
    
    species_counts = Counter()
    for taxid in pred_species_ids:
        species_counts[taxid] = int(pred_species[taxid]['abundance'])
    
    total_reads = sum(genus_counts.values())
    if total_reads == 0:
        total_reads = sum(species_counts.values())
    if total_reads == 0:
        total_reads = 1000000
    
    # Calculate metrics
    genus_metrics = calc_metrics(genus_counts, gt_genus, gt_genus_ids, total_reads)
    species_metrics = calc_metrics(species_counts, gt_species, gt_species_ids, total_reads)
    
    return {
        'dataset': dataset,
        'genus': genus_metrics,
        'species': species_metrics,
        'total_reads': total_reads
    }


def main():
    global GT_GENUS_DIR, GT_SPECIES_DIR
    
    parser = argparse.ArgumentParser(
        description='Benchmark classification tools against ground truth'
    )
    parser.add_argument(
        '--ground-truth-dir',
        required=True,
        help='Directory containing ground truth files (with genus/ and species/ subdirs)'
    )
    parser.add_argument(
        '--eval-dir',
        required=True,
        help='Directory containing tool outputs to benchmark'
    )
    parser.add_argument(
        '--files-to-evaluate',
        nargs='+',
        help='Optional list of dataset names to evaluate (without _TRUTH.txt suffix)'
    )
    
    args = parser.parse_args()
    
    GT_GENUS_DIR = os.path.join(args.ground_truth_dir, "genus")
    GT_SPECIES_DIR = os.path.join(args.ground_truth_dir, "species")
    
    if not os.path.isdir(GT_GENUS_DIR) or not os.path.isdir(GT_SPECIES_DIR):
        print(f"Error: Ground truth directory must have genus/ and species/ subdirectories")
        return 1
    
    if not os.path.isdir(args.eval_dir):
        print(f"Error: Evaluation directory not found: {args.eval_dir}")
        return 1
    
    # Create metrics directory
    metrics_dir = os.path.join(args.eval_dir, "metrics")
    os.makedirs(metrics_dir, exist_ok=True)
    
    # Get list of datasets to evaluate
    if args.files_to_evaluate:
        datasets = args.files_to_evaluate
    else:
        datasets = []
        for filename in os.listdir(GT_GENUS_DIR):
            if filename.endswith('_TRUTH.txt'):
                datasets.append(filename.replace('_TRUTH.txt', ''))
    
    if not datasets:
        print("Error: No datasets found to evaluate")
        return 1
    
    print(f"Evaluating {len(datasets)} datasets")
    
    # First pass: Convert KrakenUniq reports to ground truth format
    print("\nConverting KrakenUniq reports to ground truth format...")
    for dir_name in os.listdir(args.eval_dir):
        tool_dir = os.path.join(args.eval_dir, dir_name)
        if not os.path.isdir(tool_dir) or dir_name == 'metrics':
            continue
        
        is_krakenuniq = 'krakenuniq' in dir_name.lower()
        
        if is_krakenuniq:
            print(f"\nProcessing {dir_name}...")
            for dataset in datasets:
                convert_krakenuniq_to_ground_truth_format(tool_dir, dataset, dir_name)
    
    # Second pass: Evaluate all tools uniformly
    print("\n\nEvaluating all tools...")
    all_tool_results = {}
    
    for dir_name in os.listdir(args.eval_dir):
        tool_dir = os.path.join(args.eval_dir, dir_name)
        if not os.path.isdir(tool_dir) or dir_name == 'metrics':
            continue
        
        is_krakenuniq = 'krakenuniq' in dir_name.lower()
        
        # Check for genus/species subdirectories
        genus_dir = os.path.join(tool_dir, 'genus')
        species_dir = os.path.join(tool_dir, 'species')
        
        if not os.path.isdir(genus_dir) or not os.path.isdir(species_dir):
            print(f"\nWarning: Skipping {dir_name} - no genus/species subdirectories found")
            continue
        
        print(f"\nDiscovering tools in {dir_name}...")
        
        # Find all unique tool names from the genus directory files
        tool_names = set()
        for filename in os.listdir(genus_dir):
            if filename.endswith('.txt') and '_' in filename:
                # Extract tool name from filename (format: <dataset>_<tool>.txt)
                for dataset in datasets:
                    if filename.startswith(dataset + '_'):
                        tool_name = filename[len(dataset) + 1:].replace('.txt', '')
                        tool_names.add(tool_name)
                        break
        
        if not tool_names:
            print(f"  No tool outputs found in {dir_name}")
            continue
        
        print(f"  Found {len(tool_names)} tool(s): {', '.join(sorted(tool_names))}")
        
        # Evaluate each tool
        for tool_name in sorted(tool_names):
            print(f"\n  Evaluating {tool_name}...")
            
            tool_results = []
            
            for dataset in datasets:
                result = eval_tool_classification(dataset, tool_dir, tool_name)
                
                if result:
                    tool_results.append(result)
                    
                    # Save individual metrics only for KrakenUniq tools
                    if is_krakenuniq:
                        output_file = os.path.join(metrics_dir, f"{tool_name}_{dataset}_metrics.tsv")
                        with open(output_file, 'w') as f:
                            f.write("Metric\tGenus\tSpecies\n")
                            f.write(f"Precision\t{result['genus']['precision']:.4f}\t{result['species']['precision']:.4f}\n")
                            f.write(f"Recall\t{result['genus']['recall']:.4f}\t{result['species']['recall']:.4f}\n")
                            f.write(f"F1\t{result['genus']['f1']:.4f}\t{result['species']['f1']:.4f}\n")
                            f.write(f"AUPR\t{result['genus']['aupr']:.4f}\t{result['species']['aupr']:.4f}\n")
                            f.write(f"Total_Reads\t{result['total_reads']}\t{result['total_reads']}\n")
            
            if tool_results:
                all_tool_results[tool_name] = {
                    'n_datasets': len(tool_results),
                    'genus_precision_mean': np.mean([r['genus']['precision'] for r in tool_results]),
                    'genus_precision_std': np.std([r['genus']['precision'] for r in tool_results]),
                    'genus_recall_mean': np.mean([r['genus']['recall'] for r in tool_results]),
                    'genus_recall_std': np.std([r['genus']['recall'] for r in tool_results]),
                    'genus_f1_mean': np.mean([r['genus']['f1'] for r in tool_results]),
                    'genus_f1_std': np.std([r['genus']['f1'] for r in tool_results]),
                    'genus_aupr_mean': np.mean([r['genus']['aupr'] for r in tool_results]),
                    'genus_aupr_std': np.std([r['genus']['aupr'] for r in tool_results]),
                    'species_precision_mean': np.mean([r['species']['precision'] for r in tool_results]),
                    'species_precision_std': np.std([r['species']['precision'] for r in tool_results]),
                    'species_recall_mean': np.mean([r['species']['recall'] for r in tool_results]),
                    'species_recall_std': np.std([r['species']['recall'] for r in tool_results]),
                    'species_f1_mean': np.mean([r['species']['f1'] for r in tool_results]),
                    'species_f1_std': np.std([r['species']['f1'] for r in tool_results]),
                    'species_aupr_mean': np.mean([r['species']['aupr'] for r in tool_results]),
                    'species_aupr_std': np.std([r['species']['aupr'] for r in tool_results]),
                }
                
                print(f"    Evaluated {len(tool_results)} datasets")
                print(f"    Genus   - P: {all_tool_results[tool_name]['genus_precision_mean']:.3f}±{all_tool_results[tool_name]['genus_precision_std']:.3f}, "
                      f"R: {all_tool_results[tool_name]['genus_recall_mean']:.3f}±{all_tool_results[tool_name]['genus_recall_std']:.3f}, "
                      f"F1: {all_tool_results[tool_name]['genus_f1_mean']:.3f}±{all_tool_results[tool_name]['genus_f1_std']:.3f}, "
                      f"AUPR: {all_tool_results[tool_name]['genus_aupr_mean']:.3f}±{all_tool_results[tool_name]['genus_aupr_std']:.3f}")
                print(f"    Species - P: {all_tool_results[tool_name]['species_precision_mean']:.3f}±{all_tool_results[tool_name]['species_precision_std']:.3f}, "
                      f"R: {all_tool_results[tool_name]['species_recall_mean']:.3f}±{all_tool_results[tool_name]['species_recall_std']:.3f}, "
                      f"F1: {all_tool_results[tool_name]['species_f1_mean']:.3f}±{all_tool_results[tool_name]['species_f1_std']:.3f}, "
                      f"AUPR: {all_tool_results[tool_name]['species_aupr_mean']:.3f}±{all_tool_results[tool_name]['species_aupr_std']:.3f}")
            else:
                print(f"    No results for {tool_name}")
    
    # Write summary results to eval_dir
    output_file = os.path.join(args.eval_dir, "benchmark_summary.tsv")
    with open(output_file, 'w') as f:
        f.write("Tool\tN_Datasets\t")
        f.write("Genus_Precision_Mean\tGenus_Precision_Std\tGenus_Recall_Mean\tGenus_Recall_Std\t")
        f.write("Genus_F1_Mean\tGenus_F1_Std\tGenus_AUPR_Mean\tGenus_AUPR_Std\t")
        f.write("Species_Precision_Mean\tSpecies_Precision_Std\tSpecies_Recall_Mean\tSpecies_Recall_Std\t")
        f.write("Species_F1_Mean\tSpecies_F1_Std\tSpecies_AUPR_Mean\tSpecies_AUPR_Std\n")
        
        for tool_name, metrics in all_tool_results.items():
            f.write(f"{tool_name}\t{metrics['n_datasets']}\t")
            f.write(f"{metrics['genus_precision_mean']:.4f}\t{metrics['genus_precision_std']:.4f}\t")
            f.write(f"{metrics['genus_recall_mean']:.4f}\t{metrics['genus_recall_std']:.4f}\t")
            f.write(f"{metrics['genus_f1_mean']:.4f}\t{metrics['genus_f1_std']:.4f}\t")
            f.write(f"{metrics['genus_aupr_mean']:.4f}\t{metrics['genus_aupr_std']:.4f}\t")
            f.write(f"{metrics['species_precision_mean']:.4f}\t{metrics['species_precision_std']:.4f}\t")
            f.write(f"{metrics['species_recall_mean']:.4f}\t{metrics['species_recall_std']:.4f}\t")
            f.write(f"{metrics['species_f1_mean']:.4f}\t{metrics['species_f1_std']:.4f}\t")
            f.write(f"{metrics['species_aupr_mean']:.4f}\t{metrics['species_aupr_std']:.4f}\n")
    
    print(f"\nBenchmark summary written to: {output_file}")
    print(f"Individual metrics saved to: {metrics_dir}")
    print(f"Total tools benchmarked: {len(all_tool_results)}")
    
    # Create visualization
    plot_file = os.path.join(args.eval_dir, "benchmark_plot.png")
    plot_benchmark_results(all_tool_results, plot_file)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
