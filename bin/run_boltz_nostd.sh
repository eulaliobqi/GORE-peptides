#!/usr/bin/env bash
# Igual a run_boltz.sh, mas via boltz_nostd.py (cargas do SMILES preservadas). Uso: run_boltz_nostd.sh <dir_yamls> <out_dir> <seed>
set -uo pipefail
source ~/miniforge3/etc/profile.d/conda.sh
conda activate boltz2-env
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
IN=$1; OUT=$2; SEED=$3; shift 3
mkdir -p "$OUT"
python ~/GOREs-boltz/bin/boltz_nostd.py predict "$IN" --out_dir "$OUT" --model boltz2 \
  --diffusion_samples ${DS:-3} --recycling_steps 3 --sampling_steps 200 \
  --max_msa_seqs ${MSAMAX:-8192} --no_kernels --use_potentials \
  --output_format pdb --seed "$SEED" --num_workers 2 --preprocessing-threads ${PT:-1} "$@"
