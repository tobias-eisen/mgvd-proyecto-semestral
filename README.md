# Semestral project of the course Tópicos en Manejo de Grandes Volúmenes de Datos (MGVD) of 2025-2
### Goal: optimize efficiency in metagenomic classification for low-memory devices

To assess performance of state-of-the-art metagenomic classification tools, two streaming-based algorithms were analyzed: Bignorm (A. Wedemeyer et al., “An improved filtering algorithm for big read datasets and its application to single-cell assembly”, BMC Bioinformatics, 2017.) for normalization based on rare k-mers and KrakenUniq (F. P. Breitwieser, D. N. Baker, and S. L. Salzberg, “KrakenUniq: confident and fast metagenomics classification using unique k-mer counts,” Genome Biol, 2018.) for k-mer-based classification.

## Quick Setup

### 1. Download Data and Create Docker Volume

```bash
# Move to the location designated for the KrakenUniq std. database and benchmark dataset (Volume with at least 400Gb)
cd /path/to/metagenomics_data

# Download KrakenUniq std. database
wget -P ./standard_db \
     https://genome-idx.s3.amazonaws.com/kraken/uniq/krakendb-2022-06-16-STANDARD/database.kdb
wget -P ./standard_db \
	https://genome-idx.s3.amazonaws.com/kraken/uniq/krakendb-2022-06-16-STANDARD/kuniq_standard_minus_kdb.20220616.tgz
# unzip the .tar.gz index files
tar -xf kuniq_standard_minus_kdb.20220616.tgz 

# Download McIntyre IMMSA dataset
wget -r -nH --cut-dirs=2 -np \
     -P ./IMMSA \
     ftp://ftp.ncbi.nlm.nih.gov/nist-immsa/IMMSA/
cd IMMSA && unzip truth_sets.zip && cd .. # (Unzipping might already be done)

# Create docker volume (bind path to working directory / parent directory of IMMSA and standard_db)
docker volume create \
  --driver local \
  --opt type=none \
  --opt device=$(pwd) \
  --opt o=bind \
  metagenomics_data
```

### 2. Setup Environment

**Option A: VS Code Dev Container (Recommended)**
1. Open this folder in VS Code
2. Install the "Dev Containers" extension
3. Click "Reopen in Container" when prompted (or Ctrl+P -> "> Dev Containers: Rebuild and Reopen in Container")
4. Bignorm will auto-compile on container start

**Option B: Manual Docker Setup**
```bash
# Navigate to repository directory
cd /path/to/mgvd_proyecto-semestral

# Build Docker image from Dockerfile
docker build -t mgvd .devcontainer

# Run container
docker run -it \
	-v metagenomics_data:/mnt/data \
	-v $(pwd):/mgvd \
	--name mgvd mgvd /bin/bash

# Optional re-build of Bignorm (Bignorm executable is already in repository)
cd /mgvd/Bignorm && make && cd ..
```

## Usage

The used data directories are currently hard-coded for the mounted IMMSA volume (`/mnt/data/IMMSA`), the results are saved in `IMMSA_evaluation` in the root of the repository, where the KrakenUniq baselines and existing evaluations from the project are located.


### classify.py - Run KrakenUniq Classification

Classify metagenomic reads using KrakenUniq with various options (they can be combined).

```bash
# Basic classification (iterates with standard KrakenUniq over all IMMSA read files)
./classify.py 

# Quick mode (use first hit(s))
./classify.py --quick

# With memory preloading
./classify.py --preload-size 1G

# Process specific datasets only
./classify.py --datasets-to-classify JGI_SRR033549 ABRF_MGRG_1ng

# Use normalized data (Bignorm or custom), searches for <original_fasta_filename>.<method_suffix>.gz
./classify.py --method-suffix bignorm_default
```

**Outputs:**
- `IMMSA_evaluation/krakenuniq_*/` - Classification reports and outputs
- `runtimes.tsv` - Execution times per dataset

### bignorm_experiments.py - Bignorm Parameter Sweep

Test multiple Bignorm normalization parameters to find optimal filtering settings. For using the output together with KrakenUniq, look up the file in `/mnt/data/IMMSA` and check if the output was renamed from `.gz_keep.gz` to `.<param_set>.gz` where the default `param_set` is `bignorm_default`.

```bash
# Run all parameter combinations
./bignorm_experiments.py --experiment-dir IMMSA_evaluation/bignorm_experiments

# Run only default/baseline parameters
./bignorm_experiments.py \
  --experiment-dir IMMSA_evaluation/bignorm_experiments \
  --baseline-only

# Process single file
./bignorm_experiments.py \
  --experiment-dir IMMSA_evaluation/bignorm_experiments \
  --file /mnt/data/IMMSA//Volumes/SSDext/metagenomics_data/IMMSA/ABRF_MGRG_1ng_Repli_g_08142015_GTCCGC_L001_R1_001.fastq.gz
```

**Parameter sets tested:**
- Default (k=31, Q=20, A=3, B=5, C=20, N=10)
- Quality thresholds (Q: 10, 20, 30)
- Rarity thresholds (A: 2, 3, 5)
- Contribution thresholds (B: 3, 5, 10)
- Abundance thresholds (C: 10, 20, 50)
- N-tolerance (N: 5, 10, 20)

**Outputs:**
- `summary.tsv` - Aggregated statistics (reads kept, FP rate, runtime)
- `logs/*.txt` - Detailed logs per experiment
- `*.bignorm_*.gz` - Normalized read files

### benchmark.py - Performance Evaluation

Evaluate classification accuracy against ground truth at genus and species levels. To include runtime in the plotting from Bignorm-normalized runs, the `dataset` entry in `/mgvd/IMMSA_evaluation/bignorm_experiments/summary.tsv` might have to be shortened to the naming convention in `/mnt/data/IMMSA/truth_sets` (without `_TRUTH.txt`).

```bash
# Benchmark classification results (evaluates all tools in eval-dir)
./benchmark.py \
  --ground-truth-dir /mnt/data/IMMSA/truth_sets \
  --eval-dir IMMSA_evaluation

# Benchmark with Bignorm preprocessing times included (for combined runtime of Bignorm-normalized runs in benchmark_runtime.png)
./benchmark.py \
  --ground-truth-dir /mnt/data/IMMSA/truth_sets \
  --eval-dir IMMSA_evaluation \
  --bignorm-experiments-dir IMMSA_evaluation/bignorm_experiments

# Evaluate only specific datasets
./benchmark.py \
  --ground-truth-dir /mnt/data/IMMSA/truth_sets \
  --eval-dir IMMSA_evaluation \
  --files-to-evaluate JGI_SRR033549 ABRF_MGRG_1ng

# Only evaluate and plot KrakenUniq variations
./benchmark.py \
  --ground-truth-dir /mnt/data/IMMSA/truth_sets \
  --eval-dir IMMSA_evaluation \
  --krakenuniq-only
```

**Outputs:**
- `benchmark_summary.tsv` - Precision, recall, F1, AUPR per tool (mean ± std)
- `metrics/*_metrics.tsv` - Detailed per-dataset metrics for each tool
- `benchmark_metrics.png` - F1 and recall comparison plots
- `benchmark_runtime.png` - Runtime comparison with Bignorm preprocessing times

## Directory Structure

```
/mgvd/                          # Working directory
├── Bignorm/                    # Bignorm source code
├── classify.py                 # KrakenUniq classification script
├── bignorm_experiments.py      # Bignorm parameter sweep
├── benchmark.py                # Accuracy benchmarking
└── IMMSA_evaluation/           # Results directory
    ├── krakenuniq_*/           # Classification outputs
    ├── bignorm_experiments/    # Normalization experiments
    └── benchmark_summary.tsv   # Performance metrics

/mnt/data/                      # Mounted Docker volume
├── IMMSA/                      # Metagenomic datasets
│   └── truth_sets/             # Ground truth taxonomies
└── standard_db/                # KrakenUniq database
```

## Notes

- Scripts auto-detect paired-end reads (R1/R2 naming)
- FASTQ files use decision function 6 (quality-aware), FASTA uses function 3
- All scripts support incremental runs (skip completed experiments)