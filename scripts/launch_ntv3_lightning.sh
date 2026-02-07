#!/usr/bin/env bash
#
# Launch the NTv3 TE fossil scoring pipeline on Lightning.ai
#
# Prerequisites:
#   pip install lightning-sdk
#
# Auth (pick one):
#   Option 1 (interactive, one-time): lightning login
#   Option 2 (env vars):
#     export LIGHTNING_USER_ID=<your-user-id>
#     export LIGHTNING_API_KEY=<your-api-key>
#     (find these at lightning.ai → Settings → Keys → Programmatic Login)
#
# Usage:
#   bash scripts/launch_ntv3_lightning.sh          # full run
#   bash scripts/launch_ntv3_lightning.sh --test    # 50 regions only
#
set -euo pipefail

# Activate the conda env that has lightning-sdk + torch + transformers
source activate te-contrastive 2>/dev/null || conda activate te-contrastive 2>/dev/null || true

STUDIO_NAME="ntv3-te-fossils"
TEAMSPACE="default"          # change to your teamspace name
MACHINE="A10G"               # A10G (24GB, Pro) or A100 (Teams plan)
REMOTE_DIR="/teamspace/studios/this_studio"

# ---------------------------------------------------------------------------
# Parse args
# ---------------------------------------------------------------------------
TEST_MODE=false
if [[ "${1:-}" == "--test" ]]; then
    TEST_MODE=true
    echo "=== TEST MODE: 50 regions only ==="
fi

# ---------------------------------------------------------------------------
# 1. Create + start studio
# ---------------------------------------------------------------------------
echo "=== Creating/starting studio: $STUDIO_NAME ==="
python3 -c "
from lightning_sdk import Studio, Machine
s = Studio('${STUDIO_NAME}', '${TEAMSPACE}', create_ok=True)
s.start(Machine.${MACHINE})
print('Studio running on', '${MACHINE}')
"

# ---------------------------------------------------------------------------
# 2. Get SSH target from studio
# ---------------------------------------------------------------------------
echo "=== Getting SSH config ==="
# lightning CLI sets up SSH config automatically; the host alias is the studio name
SSH_TARGET="${STUDIO_NAME}"

# Verify SSH works
ssh "${SSH_TARGET}" "echo 'SSH OK'" || {
    echo "SSH failed. Run 'lightning ssh-setup' or check Lightning.ai SSH docs."
    echo "Your SSH target may be s_<SESSION_ID>@ssh.lightning.ai"
    echo "Check the studio settings in the Lightning.ai web UI."
    exit 1
}

# ---------------------------------------------------------------------------
# 3. Clone repo + transfer data
# ---------------------------------------------------------------------------
echo "=== Setting up remote environment ==="
ssh "${SSH_TARGET}" "
    cd ${REMOTE_DIR}
    if [ ! -d repeat_finder ]; then
        git clone git@github.com:jacobboysen/repeat_finder.git
    else
        cd repeat_finder && git pull
    fi
    mkdir -p repeat_finder/results repeat_finder/data/references
"

echo "=== Transferring dataset (691 MB) ==="
rsync -avz --progress \
    results/te_fossil_lm_dataset/ \
    "${SSH_TARGET}:${REMOTE_DIR}/repeat_finder/results/te_fossil_lm_dataset/"

echo "=== Transferring genome FASTA (139 MB) ==="
rsync -avz --progress \
    data/references/dmel_genome.fasta \
    "${SSH_TARGET}:${REMOTE_DIR}/repeat_finder/data/references/dmel_genome.fasta"

# ---------------------------------------------------------------------------
# 4. Install dependencies + run tests
# ---------------------------------------------------------------------------
echo "=== Installing dependencies ==="
ssh "${SSH_TARGET}" "
    cd ${REMOTE_DIR}/repeat_finder
    pip install torch transformers pyfaidx numpy pandas
"

echo "=== Running tests ==="
ssh "${SSH_TARGET}" "
    cd ${REMOTE_DIR}/repeat_finder
    python scripts/test_ntv3_pipeline.py
"

# ---------------------------------------------------------------------------
# 5. Run the pipeline
# ---------------------------------------------------------------------------
echo "=== Running NTv3 scoring pipeline ==="
if $TEST_MODE; then
    ssh "${SSH_TARGET}" "
        cd ${REMOTE_DIR}/repeat_finder
        python scripts/score_te_fossils_ntv3.py -v --max-regions 50 --satmut --top-n 10
    "
else
    # Use nohup so it survives SSH disconnects
    ssh "${SSH_TARGET}" "
        cd ${REMOTE_DIR}/repeat_finder
        nohup python scripts/score_te_fossils_ntv3.py -v --satmut --top-n 100 \
            > ntv3_run.log 2>&1 &
        echo 'Pipeline running in background. PID:' \$!
        echo 'Monitor with: tail -f ntv3_run.log'
    "
    echo ""
    echo "Pipeline launched in background on ${MACHINE}."
    echo "Monitor:  ssh ${SSH_TARGET} 'tail -f ${REMOTE_DIR}/repeat_finder/ntv3_run.log'"
    echo "Results:  rsync -avz ${SSH_TARGET}:${REMOTE_DIR}/repeat_finder/results/ntv3_te_scoring/ results/ntv3_te_scoring/"
fi

# ---------------------------------------------------------------------------
# 6. (After completion) Pull results back
# ---------------------------------------------------------------------------
echo ""
echo "=== When done, pull results with: ==="
echo "rsync -avz ${SSH_TARGET}:${REMOTE_DIR}/repeat_finder/results/ntv3_te_scoring/ results/ntv3_te_scoring/"
echo ""
echo "=== Then stop the studio with: ==="
echo "python3 -c \"from lightning_sdk import Studio; Studio('${STUDIO_NAME}', '${TEAMSPACE}').stop()\""
