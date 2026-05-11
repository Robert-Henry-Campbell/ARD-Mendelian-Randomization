#!/usr/bin/env bash
set -euo pipefail

: "${ARDMR_PKG_PATH:?ARDMR_PKG_PATH must be set (e.g. /mnt/sdg/robert/ardmr/ARD-Mendelian-Randomization)}"
: "${ARDMR_CACHE_DIR:?ARDMR_CACHE_DIR must be set (tens-of-GB scratch dir)}"
: "${OPENGWAS_JWT:?OPENGWAS_JWT must be set}"
export ARDMR_PKG_PATH ARDMR_CACHE_DIR OPENGWAS_JWT

SCRIPT_DIR="$ARDMR_PKG_PATH/inst/scripts/19batch"
LOG_DIR="$ARDMR_CACHE_DIR/logs_19batch"
mkdir -p "$LOG_DIR"

ts() { date '+%Y-%m-%d %H:%M:%S'; }

echo "[$(ts)] Pre-warming cache via batch 08 (serial)..."
Rscript "$SCRIPT_DIR/run_batch_08.R" \
    > "$LOG_DIR/batch_08.log" 2>&1
echo "[$(ts)] Pre-warm complete."

PARALLEL_BATCHES=(01 02 03 04 05 06 07 09 10 11 12 13 14 15 16 17 18 19)
declare -a pids=()
declare -a names=()

for nn in "${PARALLEL_BATCHES[@]}"; do
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

echo "[$(ts)] Summary: $((${#pids[@]} - fail))/${#pids[@]} parallel batches succeeded; $fail failed."
echo "[$(ts)] Pre-warm batch 08 log: $LOG_DIR/batch_08.log"
exit "$fail"
