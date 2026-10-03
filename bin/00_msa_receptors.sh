#!/usr/bin/env bash
# Etapa 0: dobra cada receptor sozinho com --use_msa_server UMA vez (MSA fica em cache para todo o lote).
# O dobramento de receptor isolado tambem alimenta o QC de receptor (etapa 1).
# Autorizado pelo usuario em 30/09 para o DN1937 (sequencia inedita, enviada uma unica vez).
set -euo pipefail
source ~/miniforge3/etc/profile.d/conda.sh
conda activate boltz2-env
ROOT=~/GOREs-boltz
mkdir -p $ROOT/inputs/receptors $ROOT/runs/receptors
python - <<'PY'
import os
root = os.path.expanduser('~/GOREs-boltz')
for l in open(f'{root}/data/private/receptors.tsv').read().splitlines()[1:]:
    name, seq = l.split('\t')
    open(f'{root}/inputs/receptors/{name}.yaml', 'w').write(
        f"version: 1\nsequences:\n- protein:\n    id: A\n    sequence: {seq}\n")
PY
boltz predict $ROOT/inputs/receptors --out_dir $ROOT/runs/receptors \
  --use_msa_server --diffusion_samples 5 --recycling_steps 3 --no_kernels --use_potentials \
  --output_format pdb --seed 101 2>&1 | tee $ROOT/logs/00_msa_receptors.log
