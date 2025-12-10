#!/opt/conda/bin/python3
"""
Evaluate KrakenUniq classification results against ground truth.

This script:
1. Parses KrakenUniq report and classification output
2. Formats genus and species output
3. Calculates precision, recall, AUPR, and F1 score at species and genus levels
"""

import sys
import argparse
from collections import defaultdict, Counter
import numpy as np
from sklearn.metrics import precision_recall_curve, auc, precision_score, recall_score, f1_score


def parse_krakenuniq_report(report_file):
    """
    Parse KrakenUniq REPORTFILE.tsv to get taxonomy information.
    
    Returns:
        dict: taxid -> {rank, name, reads, taxReads}
    """
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
    """
    Extract all taxa at a specific rank from the taxonomy report.
    
    Args:
        taxonomy: dict of taxid -> taxonomy info from report
        rank: target rank (e.g., 'species', 'genus')
    
    Returns:
        Counter: taxid -> read count for taxa at the specified rank
    """
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


def parse_ground_truth(gt_file):
    """
    Parse ground truth file.
    
    Format: taxid  fraction  abundance  rank  name
    
    Returns:
        dict: taxid -> fraction, set of taxids
    """
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
    """
    Calculate precision, recall, F1 score, and AUPR.
    
    Args:
        pred_counts: Counter of predicted taxid -> read count
        gt_taxa: dict of ground truth taxid -> info
        gt_taxids: set of ground truth taxids
        total_reads: total number of reads
    
    Returns:
        dict with metrics
    """
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


def main():
    parser = argparse.ArgumentParser(
        description='Evaluate KrakenUniq classification against ground truth'
    )
    parser.add_argument(
        '--report',
        required=True,
        help='KrakenUniq REPORTFILE.tsv'
    )
    parser.add_argument(
        '--gt-species',
        default='gt_classification_species.tsv',
        help='Ground truth species classification'
    )
    parser.add_argument(
        '--gt-genus',
        default='gt_classification_genus.tsv',
        help='Ground truth genus classification'
    )
    parser.add_argument(
        '--output-species',
        default='krakenuniq_classification_species.tsv',
        help='Output file for formatted species results'
    )
    parser.add_argument(
        '--output-genus',
        default='krakenuniq_classification_genus.tsv',
        help='Output file for formatted genus results'
    )
    
    args = parser.parse_args()
    
    print("=== KrakenUniq Evaluation ===\n")
    
    # Parse KrakenUniq report
    print(f"Parsing KrakenUniq report: {args.report}")
    taxonomy = parse_krakenuniq_report(args.report)
    print(f"  Found {len(taxonomy)} taxonomic entries")
    
    # Calculate total reads from the report
    total_reads = taxonomy.get(1, {}).get('reads', 0)  # Root taxid=1 has all classified reads
    if total_reads == 0:
        # Fallback: sum all reads at root level
        total_reads = sum(info['reads'] for info in taxonomy.values() if info.get('reads', 0) > 0)
    print(f"  Total reads in report: {total_reads}")
    
    # Process species level
    print("\n--- Species Level ---")
    species_counts = extract_taxa_by_rank(taxonomy, 'species')
    print(f"Found {len(species_counts)} species in predictions")
    
    # Format species results
    species_formatted = format_like_ground_truth(
        species_counts, taxonomy, total_reads, 'species'
    )
    
    # Write formatted species output
    with open(args.output_species, 'w') as f:
        for taxid, fraction, abundance, rank, name in species_formatted:
            f.write(f"{taxid}\t{fraction:.5f}\t{abundance:.5f}\t{rank}\t{name}\n")
    print(f"Wrote formatted species results to: {args.output_species}")
    
    # Parse ground truth species
    gt_species, gt_species_ids = parse_ground_truth(args.gt_species)
    print(f"Ground truth contains {len(gt_species_ids)} species")
    
    # Calculate species metrics
    species_metrics = calculate_metrics(
        species_counts, gt_species, gt_species_ids, total_reads
    )
    
    print("\nSpecies Metrics:")
    print(f"  Precision: {species_metrics['precision']:.4f}")
    print(f"  Recall:    {species_metrics['recall']:.4f}")
    print(f"  F1 Score:  {species_metrics['f1']:.4f}")
    print(f"  AUPR:      {species_metrics['aupr']:.4f}")
    print(f"  True Positives:  {species_metrics['true_positives']}")
    print(f"  False Positives: {species_metrics['false_positives']}")
    print(f"  False Negatives: {species_metrics['false_negatives']}")
    
    # Process genus level
    print("\n--- Genus Level ---")
    genus_counts = extract_taxa_by_rank(taxonomy, 'genus')
    print(f"Found {len(genus_counts)} genera in predictions")
    
    # Format genus results
    genus_formatted = format_like_ground_truth(
        genus_counts, taxonomy, total_reads, 'genus'
    )
    
    # Write formatted genus output
    with open(args.output_genus, 'w') as f:
        for taxid, fraction, abundance, rank, name in genus_formatted:
            f.write(f"{taxid}\t{fraction:.5f}\t{abundance:.5f}\t{rank}\t{name}\n")
    print(f"Wrote formatted genus results to: {args.output_genus}")
    
    # Parse ground truth genus
    gt_genus, gt_genus_ids = parse_ground_truth(args.gt_genus)
    print(f"Ground truth contains {len(gt_genus_ids)} genera")
    
    # Calculate genus metrics
    genus_metrics = calculate_metrics(
        genus_counts, gt_genus, gt_genus_ids, total_reads
    )
    
    print("\nGenus Metrics:")
    print(f"  Precision: {genus_metrics['precision']:.4f}")
    print(f"  Recall:    {genus_metrics['recall']:.4f}")
    print(f"  F1 Score:  {genus_metrics['f1']:.4f}")
    print(f"  AUPR:      {genus_metrics['aupr']:.4f}")
    print(f"  True Positives:  {genus_metrics['true_positives']}")
    print(f"  False Positives: {genus_metrics['false_positives']}")
    print(f"  False Negatives: {genus_metrics['false_negatives']}")
    
    # Summary
    print("\n=== Summary ===")
    print(f"Total Reads: {total_reads}")
    print(f"\nSpecies - P: {species_metrics['precision']:.3f}, R: {species_metrics['recall']:.3f}, F1: {species_metrics['f1']:.3f}, AUPR: {species_metrics['aupr']:.3f}")
    print(f"Genus   - P: {genus_metrics['precision']:.3f}, R: {genus_metrics['recall']:.3f}, F1: {genus_metrics['f1']:.3f}, AUPR: {genus_metrics['aupr']:.3f}")


if __name__ == '__main__':
    main()
