#!/usr/bin/env bash
# Overnight Phase 4f production runner.
# Launched by launchd (com.egor.strc-phase4f.plist) at 00:00 local time.
#
# Safety guards:
#   - caffeinate -i: prevent system sleep (display can sleep, system can't)
#   - timeout 11h hard limit: kill the run before morning regardless of state
#   - log to a timestamped file in research/strc/logs/
#   - only write the final JSON via the python script's own success path
#   - taildrop the log + (if exists) the result JSON to iPhone on completion
#
# Manual run (for testing): bash ~/Brain/research/strc/models/run_phase4f_overnight.sh

set -u
TS=$(date -u +%Y%m%d-%H%M%S)
LOG_DIR=/Users/egorlyfar/Brain/research/strc/logs
LOG=$LOG_DIR/phase4f-overnight-${TS}.log
RESULTS_JSON=/Users/egorlyfar/Brain/research/strc/models/pharmacochaperone_phase4f_interface_rescue.json
SCRIPT=/Users/egorlyfar/Brain/research/strc/models/pharmacochaperone_phase4f_interface_rescue.py
PYTHON=/opt/miniconda3/envs/strc-mmgbsa/bin/python
TIMEOUT_SECONDS=36000   # 10 hours hard cap (production estimate ~8h at 50/200 iter)
NOTIFY_DEVICE=iphone

mkdir -p "$LOG_DIR"

echo "[overnight runner] start ts=$TS log=$LOG" | tee -a "$LOG"
echo "[overnight runner] timeout=${TIMEOUT_SECONDS}s python=$PYTHON" | tee -a "$LOG"

# Production iter counts. Per-CPU smoke timing extrapolation:
#   APO  50 iter  → ~22 min for all 3 chains
#   BOUND 200 iter → ~70 min per lead (complex+chainA; chain B cached)
#   Total: 22 + 5×70 = 372 min ≈ 6.2 hours for 5 leads.
# Hard timeout 10h gives generous margin for OS jitter + thermal throttling.
# Override via env if launching manually with different precision goals.
export PHASE4F_MINIMISE=1
export PHASE4F_PLATFORM=CPU
export PHASE4F_APO_ITER=${PHASE4F_APO_ITER:-50}
export PHASE4F_BOUND_ITER=${PHASE4F_BOUND_ITER:-200}
unset PHASE4F_SMOKE

START_EPOCH=$(date +%s)
status=running

# caffeinate -i: prevent idle sleep. -d would also block display sleep.
caffeinate -i /usr/bin/perl -e 'alarm shift; exec @ARGV' "$TIMEOUT_SECONDS" \
    "$PYTHON" "$SCRIPT" >> "$LOG" 2>&1
exit_code=$?
END_EPOCH=$(date +%s)
elapsed=$((END_EPOCH - START_EPOCH))

if [[ $exit_code -eq 0 ]]; then
    status=success
elif [[ $exit_code -eq 142 || $exit_code -eq 14 ]]; then
    status=timeout
else
    status="error_exit_${exit_code}"
fi

echo "[overnight runner] done status=$status elapsed=${elapsed}s exit=$exit_code" | tee -a "$LOG"

# Push log + result to iPhone via Tailscale Taildrop.
if /usr/local/bin/tailscale status >/dev/null 2>&1; then
    echo "[overnight runner] sending log to ${NOTIFY_DEVICE}" | tee -a "$LOG"
    /usr/local/bin/tailscale file cp "$LOG" "${NOTIFY_DEVICE}:" 2>>"$LOG" || true
    if [[ -f "$RESULTS_JSON" && $status == "success" ]]; then
        /usr/local/bin/tailscale file cp "$RESULTS_JSON" "${NOTIFY_DEVICE}:" 2>>"$LOG" || true
    fi
else
    echo "[overnight runner] tailscale not available; skipping push" | tee -a "$LOG"
fi

exit $exit_code
