#!/opt/conda/bin/python3
import argparse
from collections import defaultdict, Counter
import numpy as np
import os
from sklearn.metrics import precision_recall_curve, auc, precision_score, recall_score, f1_score
import sys
import matplotlib
import re
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


def load_bignorm_runtimes(bignorm_experiments_dir):
    """Load Bignorm preprocessing runtimes from summary file."""
    bignorm_runtimes = {}
    
    if not bignorm_experiments_dir:
        return bignorm_runtimes
    
    summary_file = os.path.join(bignorm_experiments_dir, 'summary.tsv')
    
    if not os.path.exists(summary_file):
        print(f"Warning: Bignorm summary file not found: {summary_file}")
        return bignorm_runtimes
    
    with open(summary_file, 'r') as f:
        lines = f.readlines()
        if len(lines) <= 1:
            return bignorm_runtimes
        
        for line in lines[1:]:  # Skip header
            parts = line.strip().split('\t')
            dataset = parts[0]
            param_set = parts[1]
            runtime = int(parts[6])
            bignorm_runtimes[(dataset, param_set)] = runtime
    
    return bignorm_runtimes


def load_tool_runtimes(tool_dir, datasets):
    """Load runtimes from runtimes.tsv or log files in tool directory."""
    runtimes = {}
    runtime_file = os.path.join(tool_dir, 'runtimes.tsv')
    
    # First try to load from runtimes.tsv
    if os.path.exists(runtime_file):
        with open(runtime_file, 'r') as f:
            lines = f.readlines()
            if len(lines) > 1:
                for line in lines[1:]:  # Skip header
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        dataset = parts[0]
                        runtime = float(parts[1])
                        runtimes[dataset] = runtime
        return runtimes
    
    # If no runtimes.tsv, try to extract from log files
    log_dir = os.path.join(tool_dir, 'log')
    if not os.path.isdir(log_dir):
        return None  # No runtime data available
    
    for dataset in datasets:
        # Look for .log file and extract processing time
        log_file = os.path.join(log_dir, f"{dataset}.fastq.krakenuniq.log")
        if os.path.exists(log_file):
            with open(log_file, 'r') as f:
                for line in f:
                    # Look for "sequences processed in Xs" line
                    if 'sequences' in line and 'processed in' in line and 's (' in line:
                        # Format: "505962 sequences (123.11 Mbp) processed in 11.592s (2618.9 Kseq/m, 637.21 Mbp/m)."
                        match = re.search(r'processed in ([\d.]+)s', line)
                        if match:
                            runtime = float(match.group(1))
                            runtimes[dataset] = runtime
                            break
    
    return runtimes if runtimes else None


def plot_metrics_results(all_tool_results, output_file, datasets):
    """Create a 2-panel barplot showing F1 score and recall."""
    tools = sorted(all_tool_results.keys())
    n_tools = len(tools)
    
    # Set up the figure with 2x1 subplots (vertical alignment)
    fig, axes = plt.subplots(2, 1, figsize=(14, 12))
    fig.suptitle('Benchmark Metrics by Tool', fontsize=20, fontweight='bold')
    
    metrics = ['f1', 'recall']
    metric_titles = ['F1 Score', 'Recall']
    
    x = np.arange(n_tools)
    width = 0.35
    
    # Store handles and labels from first plot for shared legend
    legend_handles = None
    legend_labels = None
    
    for idx, (metric, title) in enumerate(zip(metrics, metric_titles)):
        ax = axes[idx]
        
        genus_means = [all_tool_results[tool][f'genus_{metric}_mean'] for tool in tools]
        genus_stds = [all_tool_results[tool][f'genus_{metric}_std'] for tool in tools]
        
        species_means = [all_tool_results[tool][f'species_{metric}_mean'] for tool in tools]
        species_stds = [all_tool_results[tool][f'species_{metric}_std'] for tool in tools]
        
        # Create bars
        bars1 = ax.bar(x - width/2, genus_means, width, label='Genus', 
                       alpha=0.8, color='steelblue')
        bars2 = ax.bar(x + width/2, species_means, width, label='Species',
                       alpha=0.8, color='coral')
        
        # Capture legend handles from first plot
        if idx == 0:
            legend_handles, legend_labels = ax.get_legend_handles_labels()
        
        # Customize subplot
        ax.set_ylabel(title, fontsize=16, fontweight='bold')
        ax.set_title(f'{title} Comparison', fontsize=18, fontweight='bold')
        ax.set_ylim(0, 1.0)
        
        ax.set_xticks(x)
        ax.set_xticklabels(tools, rotation=45, ha='right', fontsize=16)
        
        # Increase y-axis tick label size
        ax.tick_params(axis='y', labelsize=14)
        
        # Color krakenuniq tool labels red
        for i, label in enumerate(ax.get_xticklabels()):
            if 'krakenuniq' in tools[i].lower():
                label.set_color('red')
        
        ax.grid(axis='y', alpha=0.3, linestyle='--')
    
    # Add a single legend to the bottom
    if legend_handles and legend_labels:
        fig.legend(legend_handles, legend_labels, loc='lower center', fontsize=14, frameon=True, ncol=2, bbox_to_anchor=(0.5, -0.01))
    
    # Add datasets text box on the right
    datasets_text = "Evaluated Datasets:\n" + "\n".join([f"• {ds}" for ds in datasets])
    fig.text(1.02, 0.5, datasets_text, fontsize=16, va='center', ha='left',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
             transform=fig.transFigure)
    
    plt.tight_layout(rect=[0, 0.02, 1, 0.97])
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Metrics plot saved to: {output_file}")


def plot_runtime_results(all_tool_results, output_file, datasets):
    """Create a separate plot for runtime comparison."""
    # Filter to only include tools with runtime data
    tools_with_runtime = [tool for tool in sorted(all_tool_results.keys()) 
                          if 'runtime_mean' in all_tool_results[tool]]
    
    if not tools_with_runtime:
        print("No tools with runtime data to plot")
        return
    
    n_tools = len(tools_with_runtime)
    
    # Set up the figure
    fig = plt.figure(figsize=(14, 6))
    ax = plt.gca()
    
    x = np.arange(n_tools)
    width = 0.7
    
    # Separate base runtime and bignorm runtime
    base_runtime_means = [all_tool_results[tool]['runtime_mean'] for tool in tools_with_runtime]
    base_runtime_stds = [all_tool_results[tool]['runtime_std'] for tool in tools_with_runtime]
    bignorm_runtime_means = [all_tool_results[tool].get('bignorm_runtime_mean', 0) for tool in tools_with_runtime]
    bignorm_runtime_stds = [all_tool_results[tool].get('bignorm_runtime_std', 0) for tool in tools_with_runtime]
    
    # Calculate total runtimes for positioning labels
    total_runtimes = [base + bignorm for base, bignorm in zip(base_runtime_means, bignorm_runtime_means)]
    
    # Plot stacked bars
    bars_base = ax.bar(x, base_runtime_means, width, alpha=0.8, color='forestgreen', 
                       yerr=base_runtime_stds, label='KrakenUniq')
    bars_bignorm = ax.bar(x, bignorm_runtime_means, width, alpha=0.8, color='orange',
                          bottom=base_runtime_means, yerr=bignorm_runtime_stds, 
                          label='Bignorm')
    
    ax.set_ylabel('Runtime (seconds)', fontsize=16, fontweight='bold')
    ax.set_title('Runtime Comparison', fontsize=18, fontweight='bold')
    
    # Set y-axis limits to leave space for labels at the top
    max_runtime = max(total_runtimes)
    ax.set_ylim(top=max_runtime * 1.3)

    # Add time labels to the right of bars
    for i, runtime_sec in enumerate(total_runtimes):
        hours = int(runtime_sec // 3600)
        minutes = int((runtime_sec % 3600) // 60)
        seconds = int(runtime_sec % 60)
        
        # Build label based on what units are needed
        label_parts = []
        if hours > 0:
            label_parts.append(f"{hours}h")
        if hours > 0 or minutes > 0:
            label_parts.append(f"{minutes}m")
        label_parts.append(f"{seconds}s")
        
        label = ", ".join(label_parts)
        
        # Position label to the right of the bar
        height = runtime_sec + max_runtime * 0.1
        x_pos = i + width * 0.05
        ax.text(x_pos, height,
                label,
                ha='left', va='center', fontsize=16, rotation=0, fontweight='bold')
    
    
    ax.set_xticks(x)
    ax.set_xticklabels(tools_with_runtime, rotation=45, ha='right', fontsize=16)
    
    # Increase y-axis tick label size
    ax.tick_params(axis='y', labelsize=14)
    
    # Color krakenuniq tool labels red
    for i, label in enumerate(ax.get_xticklabels()):
        if 'krakenuniq' in tools_with_runtime[i].lower():
            label.set_color('red')
    
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    
    # Add legend
    ax.legend(loc='upper left', fontsize=14, frameon=True)
    
    # Add datasets text box on the right
    datasets_text = "Evaluated Datasets:\n" + "\n".join([f"• {ds}" for ds in datasets])
    fig.text(1.02, 0.5, datasets_text, fontsize=16, va='center', ha='left',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
             transform=fig.transFigure)
    
    plt.tight_layout(rect=[0, 0, 1, 1])
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Runtime plot saved to: {output_file}")


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
        '--bignorm-experiments-dir',
        help='Directory containing tool outputs to benchmark'
    )
    parser.add_argument(
        '--files-to-evaluate',
        nargs='+',
        help='Optional list of dataset names to evaluate (without _TRUTH.txt suffix)'
    )
    parser.add_argument(
        '--krakenuniq-only',
        action='store_true',
        help='Only evaluate and plot KrakenUniq variations (skip other tools)'
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
    
    # Load Bignorm preprocessing runtimes if directory is provided
    bignorm_runtimes = load_bignorm_runtimes(args.bignorm_experiments_dir)
    if args.bignorm_experiments_dir and bignorm_runtimes:
        print(f"Loaded {len(bignorm_runtimes)} Bignorm preprocessing runtimes from {args.bignorm_experiments_dir}")
    
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
        
        # Load runtimes for this tool
        tool_runtimes = load_tool_runtimes(tool_dir, datasets)
        
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
            # Skip special files from tool benchmark
            skip_keywords = ["average", "max", "median", "priority"]
            if any(keyword in filename for keyword in skip_keywords) or \
               filename.endswith('_TRUTH.txt') or filename.endswith('_COMMUNITY.txt'):
                continue
            
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
            # Skip non-KrakenUniq tools if --krakenuniq-only flag is set
            if args.krakenuniq_only and not is_krakenuniq:
                continue
            
            print(f"\n  Evaluating {tool_name}...")
            
            tool_results = []
            tool_runtimes_list = []
            tool_bignorm_runtimes_list = []
            
            for dataset in datasets:
                result = eval_tool_classification(dataset, tool_dir, tool_name)
                
                # Get runtime for this dataset
                if tool_runtimes is not None and dataset in tool_runtimes:
                    dataset_runtime = tool_runtimes[dataset]
                    bignorm_prep_time = 0
                    
                    # Check if this tool uses Bignorm preprocessing
                    match = re.search(r'bignorm_(\w+)', dir_name)
                    if match:
                        bignorm_param_set = match.group(0).replace("_quick", "")
                        print(bignorm_param_set)
                        # Get Bignorm preprocessing time if available
                        if bignorm_runtimes and (dataset, bignorm_param_set) in bignorm_runtimes:
                            bignorm_prep_time = bignorm_runtimes[(dataset, bignorm_param_set)]
                            print(f"    {dataset}: Bignorm preprocessing time: {bignorm_prep_time:.0f}s for {bignorm_param_set}")
                    
                    tool_runtimes_list.append(dataset_runtime)
                    tool_bignorm_runtimes_list.append(bignorm_prep_time)
                
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
            
            # Only include tools that have results for ALL datasets
            if tool_results and len(tool_results) == len(datasets):
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
                
                # Add runtime statistics if available
                if tool_runtimes_list:
                    all_tool_results[tool_name]['runtime_mean'] = np.mean(tool_runtimes_list)
                    all_tool_results[tool_name]['runtime_std'] = np.std(tool_runtimes_list)
                    all_tool_results[tool_name]['bignorm_runtime_mean'] = np.mean(tool_bignorm_runtimes_list)
                    all_tool_results[tool_name]['bignorm_runtime_std'] = np.std(tool_bignorm_runtimes_list)
                    total_mean = all_tool_results[tool_name]['runtime_mean'] + all_tool_results[tool_name]['bignorm_runtime_mean']
                    print(f"    Runtime - Mean: {total_mean:.0f}s (Base: {all_tool_results[tool_name]['runtime_mean']:.0f}s, Bignorm: {all_tool_results[tool_name]['bignorm_runtime_mean']:.0f}s)")
                else:
                    print(f"    No runtime data available for {tool_name}")
                
                print(f"    Evaluated {len(tool_results)} datasets")
                print(f"    Genus   - P: {all_tool_results[tool_name]['genus_precision_mean']:.3f}±{all_tool_results[tool_name]['genus_precision_std']:.3f}, "
                      f"R: {all_tool_results[tool_name]['genus_recall_mean']:.3f}±{all_tool_results[tool_name]['genus_recall_std']:.3f}, "
                      f"F1: {all_tool_results[tool_name]['genus_f1_mean']:.3f}±{all_tool_results[tool_name]['genus_f1_std']:.3f}, "
                      f"AUPR: {all_tool_results[tool_name]['genus_aupr_mean']:.3f}±{all_tool_results[tool_name]['genus_aupr_std']:.3f}")
                print(f"    Species - P: {all_tool_results[tool_name]['species_precision_mean']:.3f}±{all_tool_results[tool_name]['species_precision_std']:.3f}, "
                      f"R: {all_tool_results[tool_name]['species_recall_mean']:.3f}±{all_tool_results[tool_name]['species_recall_std']:.3f}, "
                      f"F1: {all_tool_results[tool_name]['species_f1_mean']:.3f}±{all_tool_results[tool_name]['species_f1_std']:.3f}, "
                      f"AUPR: {all_tool_results[tool_name]['species_aupr_mean']:.3f}±{all_tool_results[tool_name]['species_aupr_std']:.3f}")
            elif tool_results:
                print(f"    Skipping {tool_name} - only has results for {len(tool_results)}/{len(datasets)} datasets")
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
    
    # Create visualizations
    metrics_plot_file = os.path.join(args.eval_dir, "benchmark_metrics.png")
    plot_metrics_results(all_tool_results, metrics_plot_file, datasets)
    
    runtime_plot_file = os.path.join(args.eval_dir, "benchmark_runtime.png")
    plot_runtime_results(all_tool_results, runtime_plot_file, datasets)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
