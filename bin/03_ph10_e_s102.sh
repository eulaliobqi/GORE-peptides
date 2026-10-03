#!/usr/bin/env bash
# Fase 1b: pH 10 (braço S, só ionização do ligante) + semente 102 (ruído de base). Rodar dentro de `screen -S gore-ph10`.
# Ordem por prioridade: pH 10 s101 -> variantes Kn/Ym s101 -> S 8,2 s102 (piso de ruido) -> pH 10 s102 -> P 8,2 s102.
# Blocos de 8 YAMLs; se um bloco falhar (segfault no pre-processamento), refaz cada YAML sozinho (igual a 02_S_chunks.sh).
# Flags identicas as da fase 1: braço S DS=3, braço P DS=5, PT=1, sem --override (retomavel: predicoes prontas sao puladas).
cd ~/GOREs-boltz || exit 1
mkdir -p logs

run_set() {  # run_set <dir_inputs> <nome_run> <seed> <DS>
  local IN=$1 NAME=$2 SEED=$3 DS=$4
  local OUT=runs/${NAME}_s${SEED} CH=inputs/chunks/${NAME}
  local TAG=${NAME}_s${SEED}
  if [ -f logs/$TAG.done ]; then echo "$TAG ja concluido, pulando" >> logs/03.status; return; fi
  rm -rf $CH; mkdir -p $CH $OUT
  local i=0 n=0
  for f in $(ls $IN/*.yaml); do
    d=$CH/c$(printf %02d $i); mkdir -p $d; ln -sf $PWD/$f $d/; n=$((n+1))
    [ $((n % 8)) -eq 0 ] && i=$((i+1))
  done
  date +%s > logs/$TAG.t0
  for d in $CH/c*; do
    c=$(basename $d)
    PT=1 DS=$DS bin/run_boltz.sh $d $OUT $SEED > logs/${TAG}_$c.log 2>&1
    rc=$?
    ok=$(ls $OUT/boltz_results_$c/predictions 2>/dev/null | wc -l); tot=$(ls $d | wc -l)
    echo "$TAG $c rc=$rc feitos=$ok de $tot" >> logs/03.status
    if [ $ok -lt $tot ]; then
      for y in $d/*.yaml; do
        nm=$(basename $y .yaml)
        [ -d $OUT/boltz_results_$c/predictions/$nm ] && continue
        s=inputs/single/$TAG/$nm; mkdir -p $s; ln -sf $(readlink -f $y) $s/
        PT=1 DS=$DS bin/run_boltz.sh $s $OUT $SEED > logs/${TAG}_single_$nm.log 2>&1
        echo "$TAG single $nm rc=$?" >> logs/03.status
      done
    fi
  done
  date +%s > logs/$TAG.t1
  echo "$TAG FIM $(date '+%F %T')" >> logs/03.status
  touch logs/$TAG.done
}

echo "INICIO $(date '+%F %T')" >> logs/03.status
run_set inputs/S_pH10   S_pH10   101 3
run_set inputs/S_pH10Kn S_pH10Kn 101 3
run_set inputs/S_pH10Ym S_pH10Ym 101 3
run_set inputs/S        S        102 3
run_set inputs/S_pH10   S_pH10   102 3
run_set inputs/P        P        102 5
echo "TUDO FIM $(date '+%F %T')" >> logs/03.status
touch logs/03.all.done
