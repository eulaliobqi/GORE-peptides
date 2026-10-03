"""Conjuntos com ligante CARREGADO (rodados com boltz_nostd): 3 receptores nativos x 17 ligantes (13 alvos + N1-N4).
  C82   = SMILES zwitteriônico de 8,2 (os mesmos YAMLs de inputs/S)        -> 51
  C10   = pH 10 majoritário (N-term neutro; Lys+; Tyr neutra)               -> 51
  C10Kn = pH 10 com Lys épsilon neutra (ligantes com Lys; lista lida de data/private)  -> 15
  C10Ym = pH 10 com Tyr fenolato (ligantes com Tyr; lista lida de data/private) -> 9
Fora (de propósito): D189G, Pig/X1-X3 e benzamidina (controles; o pedido é ranking e mapa nos 3 receptores nativos)."""
import os, glob, shutil
ROOT = os.path.expanduser('~/GOREs-boltz')
RECS = ('Ag', 'Sf', 'Bov'); LIGS = {f'L{i:02d}' for i in range(1, 12)} | {'G1', 'G2', 'N1', 'N2', 'N3', 'N4'}
SETS = {'C82': 'S', 'C10': 'S_pH10', 'C10Kn': 'S_pH10Kn', 'C10Ym': 'S_pH10Ym'}
for new, src in SETS.items():
    d = f'{ROOT}/inputs/{new}'; shutil.rmtree(d, ignore_errors=True); os.makedirs(d)
    n = 0
    for f in sorted(glob.glob(f'{ROOT}/inputs/{src}/*.yaml')):
        rec, lig = os.path.basename(f)[:-5].split('__')
        if rec in RECS and lig in LIGS:
            shutil.copy(f, d); n += 1
    print(new, n)
