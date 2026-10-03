#!/usr/bin/env bash
# Uso: run_boltz.sh <dir_yamls> <out_dir> <seed> [flags extras do boltz]
# GPU pode ser compartilhada com outro job do usuario (b27-corrected). Padrao: 5 amostras de difusao.
# Se faltar memoria: DS=1 MSAMAX=2048 (varias sementes no lugar de varias amostras).
set -uo pipefail
source ~/miniforge3/etc/profile.d/conda.sh
conda activate boltz2-env
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
IN=$1; OUT=$2; SEED=$3; shift 3
mkdir -p "$OUT"
boltz predict "$IN" --out_dir "$OUT" --model boltz2 \
  --diffusion_samples ${DS:-5} --recycling_steps 3 --sampling_steps 200 \
  --max_msa_seqs ${MSAMAX:-8192} --no_kernels --use_potentials \
  --output_format pdb --seed "$SEED" --num_workers 2 --preprocessing-threads ${PT:-1} "$@"
