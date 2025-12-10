#!/bin/bash
# ==========================================================
# Run any command in the background and log its output
# Usage: ./run_bg.sh <your command and args>
# Example: ./run_bg.sh krakenuniq --db db_path --fastq input.fastq
# ==========================================================

set -e

# Ensure a command was provided
if [ $# -eq 0 ]; then
  echo "Usage: $0 <command-to-run>"
  exit 1
fi

# Get timestamp and log file setup
timestamp=$(date +"%Y%m%d_%H%M%S")
log_dir="process_logs"
mkdir -p "$log_dir"

# Extract a sensible base name for the log
first_word=$(basename "$1")
log_path="${log_dir}/${first_word}__${timestamp}.log"

# Run the command in background and redirect all output
# nohup prevents termination when the terminal closes
# "setsid" ensures it fully detaches from the shell session
nohup bash -c "$*" > "$log_path" 2>&1 &

# Capture the background PID
pid=$!

echo "Started '$*' in background (PID: $pid)"
echo "Logging to: $log_path"
