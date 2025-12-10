#!/bin/bash
# ==============================================================
# Run Bignorm + KrakenUniq on all datasets in a folder
# using background logging via bg_exe.sh
# ==============================================================
# Usage:
#   ./run_bignorm_krakenuniq.sh <dataset_dir> <krakenuniq_db>
# ==============================================================

set -e

# --- Check input arguments ---
if [ $# -ne 2 ]; then
  echo "Usage: $0 <dataset_dir> <krakenuniq_db>"
  exit 1
fi

dataset_dir="$1"
db_path="$2"

# --- Setup output directories ---
bignorm_out="bignorm_output"
krakenuniq_out="krakenuniq_output"
summary_file="run_summary.csv"

mkdir -p "$bignorm_out" "$krakenuniq_out"

# --- Initialize summary file ---
echo "Sample,Step,StartTime,EndTime,Duration(s),LogFile" > "$summary_file"

# --- Helper function to log time differences ---
log_duration() {
  local sample="$1"
  local step="$2"
  local start_time="$3"
  local end_time="$4"
  local log_file="$5"
  local duration=$((end_time - start_time))
  echo "$sample,$step,$(date -d @$start_time +'%F %T'),$(date -d @$end_time +'%F %T'),$duration,$log_file" >> "$summary_file"
}

# --- Main processing loop ---
for file in "$dataset_dir"/*.gz; do
  [ -f "$file" ] || continue

  sample=$(basename "$file")
  base="${sample%.*}"

  echo "=== Processing sample: $sample ==="

  # ------------------------------
  # Step 1: Run Bignorm
  # ------------------------------
  bignorm_start=$(date +%s)
  bignorm_output="$bignorm_out/${sample%.gz}.gz_keep.gz"
  ./bg_exe.sh ./Bignorm/Bignorm -u "$file" -z -k 31

  # Wait until Bignorm finishes (optional, if you want sequential execution)
  wait

  bignorm_end=$(date +%s)
  log_duration "$sample" "Bignorm" "$bignorm_start" "$bignorm_end"

  # ------------------------------
  # Step 2: Run KrakenUniq
  # ------------------------------
  kraken_start=$(date +%s)
  kraken_report="$krakenuniq_out/${base}_krakenuniq_report.txt"
  kraken_output="$krakenuniq_out/${base}_krakenuniq_output.fastq.gz"
  ./bg_exe.sh krakenuniq --db "$db_path" --preload-size 4G --threads 10 --report-file "$kraken_report" --output "$kraken_output" "$bignorm_output" 

  # Wait until KrakenUniq finishes (optional)
  wait

  kraken_end=$(date +%s)
  log_duration "$sample" "KrakenUniq" "$kraken_start" "$kraken_end" "$kraken_log"

done

echo "✅ All samples processed. Summary saved to: $summary_file"
