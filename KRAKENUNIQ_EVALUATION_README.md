# KrakenUniq Output Format and Evaluation

## KrakenUniq Output Files

### 1. REPORTFILE.tsv (Hierarchical Summary)
This file provides a hierarchical summary of all taxa found, with taxonomy tree structure.

**Format:**
```
%  reads  taxReads  kmers  dup  cov  taxID  rank  taxName
```

**Columns:**
- `%`: Percentage of reads assigned to this taxon (including descendants)
- `reads`: Number of reads assigned to this taxon and descendants
- `taxReads`: Number of reads assigned directly to this taxon (not descendants)
- `kmers`: Number of unique k-mers mapped to this taxon
- `dup`: Duplication rate
- `cov`: Coverage estimate
- `taxID`: NCBI taxonomy ID
- `rank`: Taxonomic rank (species, genus, family, etc.)
- `taxName`: Taxonomic name (indentation shows hierarchy)

**Example:**
```
97.5  24375  58  2267212  1.02  6.727e-05  1  no rank  root
70.34 17585  12  1584496  1.04  0.0001067  1224  phylum  Proteobacteria
37.8  9451   5   876330   1.04  0.000121   1236  class   Gammaproteobacteria
```

### 2. READCLASSIFICATION.tsv (Per-Read Classification)
This file contains classification results for each individual read.

**Format:**
```
C/U  read_id  taxid  read_length  additional_info
```

**Columns:**
- Column 1: `C` (classified) or `U` (unclassified)
- Column 2: Read ID
- Column 3: Taxonomy ID assigned to the read
- Column 4: Read length in base pairs
- Column 5+: Additional information (k-mer mappings)

**Example:**
```
C  carma_52609  43989  284  Q:1
C  carma_52615  43989  292  Q:1
U  carma_99999  0      300  -
```

## Ground Truth Format

The ground truth files follow a simpler format:

```
taxid  fraction  abundance  rank  name
```

**Example (species):**
```
813   1.00000  0.04000  species  Chlamydia trachomatis
562   1.00000  0.04000  species  Escherichia coli
```

**Example (genus):**
```
810   1.00000  0.04000  genus  Chlamydia
561   1.00000  0.04000  genus  Escherichia
```

## Evaluation Script Usage

### Installation
First, install required Python packages:
```bash
pip install numpy scikit-learn
```

### Basic Usage
```bash
python3 evaluate_krakenuniq.py
```

This uses default file names:
- Input: `REPORTFILE.tsv`, `READCLASSIFICATION.tsv`
- Ground truth: `gt_classification_species.tsv`, `gt_classification_genus.tsv`
- Output: `krakenuniq_classification_species.tsv`, `krakenuniq_classification_genus.tsv`

### Custom File Paths
```bash
python3 evaluate_krakenuniq.py \
  --report REPORTFILE.tsv \
  --classification READCLASSIFICATION.tsv \
  --gt-species gt_classification_species.tsv \
  --gt-genus gt_classification_genus.tsv \
  --output-species krakenuniq_species_formatted.tsv \
  --output-genus krakenuniq_genus_formatted.tsv
```

### Output

The script will:
1. **Parse KrakenUniq outputs** and extract species/genus level classifications
2. **Format results** to match ground truth format (saved to output files)
3. **Calculate metrics** for both species and genus levels:
   - **Precision**: What fraction of predicted taxa are correct?
   - **Recall**: What fraction of true taxa were found?
   - **F1 Score**: Harmonic mean of precision and recall
   - **AUPR**: Area Under Precision-Recall Curve (accounts for confidence scores)

**Example output:**
```
=== KrakenUniq Evaluation ===

Parsing KrakenUniq report: REPORTFILE.tsv
  Found 568 taxonomic entries

Parsing KrakenUniq classifications: READCLASSIFICATION.tsv
  Total reads: 25000
  Classified reads: 24375 (97.50%)

--- Species Level ---
Found 45 species in predictions
Wrote formatted species results to: krakenuniq_classification_species.tsv
Ground truth contains 25 species

Species Metrics:
  Precision: 0.8889
  Recall:    0.8000
  F1 Score:  0.8421
  AUPR:      0.8523
  True Positives:  20
  False Positives: 5
  False Negatives: 5

--- Genus Level ---
Found 32 genera in predictions
Wrote formatted genus results to: krakenuniq_classification_genus.tsv
Ground truth contains 22 genera

Genus Metrics:
  Precision: 0.9062
  Recall:    0.8636
  F1 Score:  0.8844
  AUPR:      0.8956
  True Positives:  19
  False Positives: 3
  False Negatives: 3

=== Summary ===
Total Reads: 25000
Classified: 24375 (97.50%)

Species - P: 0.889, R: 0.800, F1: 0.842, AUPR: 0.852
Genus   - P: 0.906, R: 0.864, F1: 0.884, AUPR: 0.896
```

## Metrics Explained

### Precision
Proportion of predicted taxa that are actually present in the sample.
- **High precision**: Few false positives (taxa incorrectly predicted)
- Formula: TP / (TP + FP)

### Recall (Sensitivity)
Proportion of true taxa that were correctly identified.
- **High recall**: Few false negatives (taxa missed)
- Formula: TP / (TP + FN)

### F1 Score
Harmonic mean of precision and recall, balancing both metrics.
- Formula: 2 × (Precision × Recall) / (Precision + Recall)
- Range: 0.0 (worst) to 1.0 (best)

### AUPR (Area Under Precision-Recall Curve)
Evaluates performance across all confidence thresholds.
- Uses the fraction of reads as confidence score
- Particularly useful for imbalanced datasets
- Range: 0.0 (worst) to 1.0 (best)

## Understanding the Results

### What makes a good classifier?
- **High Precision & High Recall**: Best case - finds true taxa without false positives
- **High Precision, Low Recall**: Conservative - only reports taxa when very confident
- **Low Precision, High Recall**: Liberal - reports many taxa but includes false positives
- **High F1**: Good balance between precision and recall
- **High AUPR**: Good discrimination ability across confidence levels

### Common Issues
- **Many False Positives**: Classifier is too sensitive or database has related organisms
- **Many False Negatives**: Classifier is too conservative or missing database entries
- **Low AUPR with decent F1**: Confidence scores don't correlate well with correctness

## Tips for Interpretation

1. **Species level is harder** than genus level (more false positives/negatives expected)
2. **Check true positives**: Are the most abundant species correctly identified?
3. **Examine false positives**: Are they closely related to true positives?
4. **Look at classification rate**: High % classified is good, but check if accurate
5. **Compare AUPR vs F1**: AUPR accounts for abundance/confidence weighting
