#!/usr/bin/env bash
# Conjuntos com ligante CARREGADO (boltz_nostd) em 8,2 e pH 10: pose + afinidade (EXPLORATORIO, fora da distribuicao de treino).
# Espera o lote 03 terminar (GPU e uma so). Ordem: semente 101 (4 conjuntos) -> semente 102 (C82, C10) -> P s102 (adiado do lote 03).
# Blocos de 8 YAMLs; se um bloco falhar, refaz cada YAML sozinho (igual a 03). Retomavel: conjuntos .done sao pulados.
cd ~/GOREs-boltz || exit 1
while [ ! -f logs/03.all.done ]; do sleep 120; done

run_set() {  # run_set <dir_inputs> <nome_run> <seed> <DS> <runner>
  local IN=$1 NAME=$2 SEED=$3 DS=$4 RUN=$5
  local OUT=runs/${NAME}_s${SEED} CH=inputs/chunks/${NAME} TAG=${NAME}_s${SEED}
  if [ -f logs/$TAG.done ]; then echo "$TAG ja concluido, pulando" >> logs/05.status; return; fi
  rm -rf $CH; mkdir -p $CH $OUT
  local i=0 n=0
  for f in $(ls $IN/*.yaml); do
    d=$CH/c$(printf %02d $i); mkdir -p $d; ln -sf $PWD/$f $d/; n=$((n+1))
    [ $((n % 8)) -eq 0 ] && i=$((i+1))
  done
  date +%s > logs/$TAG.t0
  for d in $CH/c*; do
    c=$(basename $d)
    PT=1 DS=$DS $RUN $d $OUT $SEED > logs/${TAG}_$c.log 2>&1
    rc=$?
    ok=$(ls $OUT/boltz_results_$c/predictions 2>/dev/null | wc -l); tot=$(ls $d | wc -l)
    echo "$TAG $c rc=$rc feitos=$ok de $tot" >> logs/05.status
    if [ $ok -lt $tot ]; then
      for y in $d/*.yaml; do
        nm=$(basename $y .yaml)
        [ -d $OUT/boltz_results_$c/predictions/$nm ] && continue
        s=inputs/single/$TAG/$nm; mkdir -p $s; ln -sf $(readlink -f $y) $s/
        PT=1 DS=$DS $RUN $s $OUT $SEED > logs/${TAG}_single_$nm.log 2>&1
        echo "$TAG single $nm rc=$?" >> logs/05.status
      done
    fi
  done
  date +%s > logs/$TAG.t1
  echo "$TAG FIM $(date '+%F %T')" >> logs/05.status
  touch logs/$TAG.done
}

NS=bin/run_boltz_nostd.sh
echo "INICIO $(date '+%F %T')" >> logs/05.status
for S in C82 C10 C10Kn C10Ym; do run_set inputs/$S $S 101 3 $NS; done
for S in C82 C10; do run_set inputs/$S $S 102 3 $NS; done
rm -f logs/P_s102.done
run_set inputs/P P 102 5 bin/run_boltz.sh
echo "TUDO FIM $(date '+%F %T')" >> logs/05.status
touch logs/05.all.done
