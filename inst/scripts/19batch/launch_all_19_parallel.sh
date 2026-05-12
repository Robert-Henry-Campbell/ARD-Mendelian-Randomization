#!/usr/bin/env bash
# Fan out all 19 batches in parallel. Assumes the shared cache has already
# been pre-warmed via inst/templates/setup_download.R (variant manifests,
# LD reference, Neale sumstats if any sex-stratified rows, Pan-UKB tabix
# indices). No serial pre-warm step here.
set -euo pipefail

: "${ARDMR_PKG_PATH:?ARDMR_PKG_PATH must be set (e.g. /mnt/sdg/robert/ardmr/ARD-Mendelian-Randomization)}"
: "${ARDMR_CACHE_DIR:?ARDMR_CACHE_DIR must be set (tens-of-GB scratch dir)}"
: "${OPENGWAS_JWT:?OPENGWAS_JWT must be set}"
export ARDMR_PKG_PATH ARDMR_CACHE_DIR OPENGWAS_JWT

cd "$ARDMR_PKG_PATH"

SCRIPT_DIR="$ARDMR_PKG_PATH/inst/scripts/19batch"
LOG_DIR="$ARDMR_CACHE_DIR/logs_19batch"
mkdir -p "$LOG_DIR"

ts() { date '+%Y-%m-%d %H:%M:%S'; }

ALL_BATCHES=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19)
declare -a pids=()
declare -a names=()

for nn in "${ALL_BATCHES[@]}"; do
    script="$SCRIPT_DIR/run_batch_${nn}.R"
    log="$LOG_DIR/batch_${nn}.log"
    nohup Rscript "$script" > "$log" 2>&1 &
    pid=$!
    pids+=("$pid")
    names+=("batch_${nn}")
    echo "[$(ts)] Launched batch_${nn} (pid $pid -> $log)"
done

echo "[$(ts)] Waiting for ${#pids[@]} parallel jobs..."
fail=0
for i in "${!pids[@]}"; do
    if ! wait "${pids[$i]}"; then
        echo "[$(ts)] FAIL: ${names[$i]} (pid ${pids[$i]}) -- see $LOG_DIR/${names[$i]}.log"
        fail=$((fail + 1))
    else
        echo "[$(ts)] OK:   ${names[$i]}"
    fi
done

echo "[$(ts)] Summary: $((${#pids[@]} - fail))/${#pids[@]} batches succeeded; $fail failed."
exit "$fail"
