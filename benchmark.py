#!/opt/conda/bin/python3
import sys
import argparse
from collections import defaultdict, Counter
import numpy as np
from sklearn.metrics import precision_recall_curve, auc, precision_score, recall_score, f1_score


def parse_krakenuniq_report(report_file):
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
            # Tax name has leading spaces indicating hierarchy - strip them
            name = parts[8].strip()
            
            taxonomy[taxid] = {
                'rank': rank,
                'name': name,
                'reads': reads,
                'taxReads': taxReads
            }
    
    return taxonomy

def extract_taxa_by_rank(taxonomy, rank):
    rank_counts = Counter()
    
    for taxid, info in taxonomy.items():
        if info['rank'] == rank and info['reads'] > 0:
            # Use 'reads' which includes reads assigned to this taxon and descendants
            rank_counts[taxid] = info['reads']
    
    return rank_counts

def format_like_ground_truth(rank_counts, taxonomy, total_reads, rank):
    """
    Format KrakenUniq results like ground truth files.
    
    Returns:
        list of tuples: (taxid, fraction, abundance, rank, name)
    """
    results = []
    
    for taxid, count in rank_counts.most_common():
        if taxid not in taxonomy:
            continue
        
        fraction = count / total_reads if total_reads > 0 else 0
        # Abundance is typically reads/total in ground truth
        abundance = fraction
        name = taxonomy[taxid]['name']
        
        results.append((taxid, fraction, abundance, rank, name))
    
    return results

def parse_ground_truth_format(gt_file):
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

def calculate_metrics(pred_counts, gt_taxa, gt_taxids, total_reads):
    # Get all unique taxids (union of predicted and ground truth)
    all_taxids = set(pred_counts.keys()) | gt_taxids
    
    # Build binary vectors and scores
    y_true = []
    y_pred = []
    scores = []
    
    for taxid in all_taxids:
        # Ground truth: 1 if taxid is in ground truth, 0 otherwise
        y_true.append(1 if taxid in gt_taxids else 0)
        
        # Prediction: 1 if taxid was predicted, 0 otherwise
        y_pred.append(1 if taxid in pred_counts else 0)
        
        # Score: fraction of reads assigned to this taxid
        score = pred_counts.get(taxid, 0) / total_reads if total_reads > 0 else 0
        scores.append(score)
    
    y_true = np.array(y_true)
    y_pred = np.array(y_pred)
    scores = np.array(scores)
    
    # Calculate metrics
    if len(y_true) == 0 or np.sum(y_true) == 0:
        return {
            'precision': 0.0,
            'recall': 0.0,
            'f1': 0.0,
            'aupr': 0.0,
            'true_positives': 0,
            'false_positives': 0,
            'false_negatives': 0,
            'support': 0
        }
    
    # Precision, Recall, F1
    precision = precision_score(y_true, y_pred, zero_division=0)
    recall = recall_score(y_true, y_pred, zero_division=0)
    f1 = f1_score(y_true, y_pred, zero_division=0)
    
    # AUPR (Area Under Precision-Recall Curve)
    if np.sum(y_true) > 0 and len(np.unique(scores)) > 1:
        precisions, recalls, _ = precision_recall_curve(y_true, scores)
        aupr = auc(recalls, precisions)
    else:
        aupr = 0.0
    
    # Confusion matrix components
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
        'false_negatives': int(fn),
        'support': int(np.sum(y_true))
    }


def evaluate_tool_directory(tool_dir, gt_genus_dir, gt_species_dir, parser_func, tool_name, eval_files=None):
    """
    Evaluate a tool's output directory against ground truth.
    
    Args:
        tool_dir: Directory containing tool output files
        gt_genus_dir: Directory containing genus-level ground truth
        gt_species_dir: Directory containing species-level ground truth
        parser_func: Function to parse tool output files
        tool_name: Name of the tool for display
        eval_files: Optional list of dataset names to evaluate (without _TRUTH.txt suffix)
    
    Returns:
        dict: Average metrics across all evaluated files
    """
    import os
    
    all_results = []
    
    # Iterate through ground truth files
    for gt_filename in os.listdir(gt_genus_dir):
        if gt_filename.startswith('.') or not gt_filename.endswith('_TRUTH.txt'):
            continue
        
        dataset_name = gt_filename.replace('_TRUTH.txt', '')
        
        # Filter by eval_files if provided
        if eval_files and dataset_name not in eval_files:
            continue
        gt_genus_fp = os.path.join(gt_genus_dir, gt_filename)
        gt_species_fp = os.path.join(gt_species_dir, gt_filename)
        
        if not os.path.exists(gt_species_fp):
            continue
        
        # Find corresponding tool output file
        tool_file = None
        for filename in os.listdir(tool_dir):
            if filename.startswith(dataset_name):
                potential_file = os.path.join(tool_dir, filename)
                if os.path.isfile(potential_file):
                    tool_file = potential_file
                    break
        
        if not tool_file:
            continue
        
        # Parse tool output
        try:
            taxonomy = parser_func(tool_file)
            
            # Calculate total reads
            total_reads = taxonomy.get(1, {}).get('reads', 0)
            if total_reads == 0:
                total_reads = sum(info.get('reads', 0) for info in taxonomy.values())
            
            # Evaluate genus level
            genus_counts = extract_taxa_by_rank(taxonomy, 'genus')
            gt_genus, gt_genus_ids = parse_ground_truth_format(gt_genus_fp)
            genus_metrics = calculate_metrics(genus_counts, gt_genus, gt_genus_ids, total_reads)
            
            # Evaluate species level
            species_counts = extract_taxa_by_rank(taxonomy, 'species')
            gt_species, gt_species_ids = parse_ground_truth_format(gt_species_fp)
            species_metrics = calculate_metrics(species_counts, gt_species, gt_species_ids, total_reads)
            
            all_results.append({
                'dataset': dataset_name,
                'genus': genus_metrics,
                'species': species_metrics
            })
            
        except Exception as e:
            print(f"Warning: Failed to process {tool_file}: {e}")
            continue
    
    if not all_results:
        return None
    
    # Calculate averages
    avg_metrics = {
        'tool': tool_name,
        'n_datasets': len(all_results),
        'genus_precision': np.mean([r['genus']['precision'] for r in all_results]),
        'genus_recall': np.mean([r['genus']['recall'] for r in all_results]),
        'genus_f1': np.mean([r['genus']['f1'] for r in all_results]),
        'genus_aupr': np.mean([r['genus']['aupr'] for r in all_results]),
        'species_precision': np.mean([r['species']['precision'] for r in all_results]),
        'species_recall': np.mean([r['species']['recall'] for r in all_results]),
        'species_f1': np.mean([r['species']['f1'] for r in all_results]),
        'species_aupr': np.mean([r['species']['aupr'] for r in all_results])
    }
    
    return avg_metrics


def parse_ground_truth_as_taxonomy(filepath):
    """
    Parse ground truth format file and convert to taxonomy dict.
    Returns taxonomy dict compatible with extract_taxa_by_rank.
    """
    taxonomy = {}
    total_reads = 0
    
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            
            try:
                taxid = int(parts[0])
                fraction = float(parts[1])
                abundance = float(parts[2])
                rank = parts[3]
                name = parts[4]
                
                # Convert fraction back to read count (approximate)
                # We'll use abundance * 1000000 as pseudo-read count
                reads = int(abundance * 1000000)
                total_reads += reads
                
                taxonomy[taxid] = {
                    'rank': rank,
                    'name': name,
                    'reads': reads,
                    'taxReads': reads
                }
            except (ValueError, IndexError):
                continue
    
    # Add root node for total reads
    if total_reads > 0:
        taxonomy[1] = {
            'rank': 'root',
            'name': 'root',
            'reads': total_reads,
            'taxReads': total_reads
        }
    
    return taxonomy


def main():
    import os
    
    parser = argparse.ArgumentParser(
        description='Benchmark classification tools against ground truth'
    )
    parser.add_argument(
        '--gt-genus-dir',
        required=True,
        help='Directory containing genus-level ground truth files'
    )
    parser.add_argument(
        '--gt-species-dir',
        required=True,
        help='Directory containing species-level ground truth files'
    )
    parser.add_argument(
        '--tool-dirs',
        nargs='+',
        required=True,
        help='Directories containing tool outputs to benchmark'
    )
    parser.add_argument(
        '--output',
        default='benchmark.tsv',
        help='Output TSV file for benchmark results'
    )
    parser.add_argument(
        '--eval-files',
        nargs='+',
        help='Optional list of dataset names to evaluate (without _TRUTH.txt suffix)'
    )
    
    args = parser.parse_args()
    
    benchmark_results = []
    
    for tool_dir in args.tool_dirs:
        if not os.path.isdir(tool_dir):
            print(f"Warning: {tool_dir} is not a directory, skipping")
            continue
        
        # Extract tool name from directory path
        tool_name = os.path.basename(tool_dir.rstrip('/'))
        
        print(f"\nEvaluating {tool_name}...")
        
        # Determine parser based on tool name
        if 'krakenuniq' in tool_name.lower():
            parser_func = parse_krakenuniq_report
        else:
            # Assume benchmark tool format (same as ground truth)
            parser_func = parse_ground_truth_as_taxonomy
        
        # Evaluate tool
        metrics = evaluate_tool_directory(
            tool_dir,
            args.gt_genus_dir,
            args.gt_species_dir,
            parser_func,
            tool_name,
            eval_files=args.eval_files
        )
        
        if metrics:
            benchmark_results.append(metrics)
            print(f"  Evaluated {metrics['n_datasets']} datasets")
            print(f"  Genus   - P: {metrics['genus_precision']:.3f}, R: {metrics['genus_recall']:.3f}, "
                  f"F1: {metrics['genus_f1']:.3f}, AUPR: {metrics['genus_aupr']:.3f}")
            print(f"  Species - P: {metrics['species_precision']:.3f}, R: {metrics['species_recall']:.3f}, "
                  f"F1: {metrics['species_f1']:.3f}, AUPR: {metrics['species_aupr']:.3f}")
        else:
            print(f"  No results found for {tool_name}")
    
    if not benchmark_results:
        print("\nNo benchmark results to write")
        return
    
    # Write results to TSV
    with open(args.output, 'w') as f:
        # Write header
        f.write("Tool\tN_Datasets\t")
        f.write("Genus_Precision\tGenus_Recall\tGenus_F1\tGenus_AUPR\t")
        f.write("Species_Precision\tSpecies_Recall\tSpecies_F1\tSpecies_AUPR\n")
        
        # Write results
        for result in benchmark_results:
            f.write(f"{result['tool']}\t{result['n_datasets']}\t")
            f.write(f"{result['genus_precision']:.4f}\t{result['genus_recall']:.4f}\t")
            f.write(f"{result['genus_f1']:.4f}\t{result['genus_aupr']:.4f}\t")
            f.write(f"{result['species_precision']:.4f}\t{result['species_recall']:.4f}\t")
            f.write(f"{result['species_f1']:.4f}\t{result['species_aupr']:.4f}\n")
    
    print(f"\nBenchmark results written to: {args.output}")
    print(f"Total tools benchmarked: {len(benchmark_results)}")


if __name__ == '__main__':
    main()
