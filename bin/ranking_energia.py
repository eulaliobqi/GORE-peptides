"""Ranking dos peptídeos por energia de ligação predita pelo Boltz-2 (braço S). Roda em qualquer máquina com results/affinity.csv.

Energia = ΔG = -1,364 * pIC50 kcal/mol (298 K), pIC50 = 6 - log10(IC50 em uM) (conversão do próprio Boltz-2; IC50 tratado
como proxy de Ki). É energia livre PREDITA por rede neural, não energia calculada fisicamente; mais negativo = mais forte.
O ligante entra NEUTRO no modelo (o Boltz-2 neutraliza o SMILES no braço de afinidade; ver Erratum em RESULTADOS_vivo.md).
Ranking só entre os 13 alvos (L01-L11, G1, G2); negativos (N1-N4) e benzamidina ficam como referência, fora do ranking.
Incerteza: média das sementes presentes; "meia-amplitude" = metade da diferença entre sementes (só com >= 2 sementes).
Saída: results/ranking_energia.csv e results/ranking_energia.md (só códigos de ligante).
"""
import os
import numpy as np, pandas as pd

RES = os.path.join(os.path.dirname(__file__), '..', 'results')
if not os.path.exists(os.path.join(RES, 'affinity.csv')):
    RES = os.path.expanduser('~/GOREs-boltz/results')
KT = 1.364
TARGETS = [f'L{i:02d}' for i in range(1, 12)] + ['G1', 'G2']
REF = ['N1', 'N2', 'N3', 'N4', 'BEN']
RECS = ['Ag', 'Sf', 'Bov']

a = pd.read_csv(os.path.join(RES, 'affinity.csv'))
a['dG'] = -KT * (6 - a.log10_ic50_uM)
g = a.groupby(['receptor', 'ligand']).dG
t = pd.DataFrame({'dG': g.mean(), 'n_sementes': g.count(), 'meia_amplitude': g.apply(lambda x: (x.max() - x.min()) / 2 if len(x) > 1 else np.nan)}).reset_index()

out = []
for r in RECS + ['Ag_D189G', 'Sf_D189G']:
    d = t[(t.receptor == r) & t.ligand.isin(TARGETS)].copy()
    d['posto'] = d.dG.rank(method='min').astype(int)
    out.append(d)
R = pd.concat(out); R.to_csv(os.path.join(RES, 'ranking_energia.csv'), index=False)

W = R[R.receptor.isin(RECS)].pivot(index='ligand', columns='receptor', values='dG').reindex(TARGETS)
P = R[R.receptor.isin(RECS)].pivot(index='ligand', columns='receptor', values='posto').reindex(TARGETS)
W['posto_medio'] = P[RECS].mean(axis=1); W = W.sort_values('posto_medio')
L = ['**ΔG predito pelo Boltz-2 (kcal/mol; mais negativo = mais forte) e posto médio nos 3 receptores nativos**', '',
     '| Ligante | ΔG Ag | ΔG Sf | ΔG Bov | Posto Ag | Posto Sf | Posto Bov | Posto médio |', '|---|---|---|---|---|---|---|---|']
for l, row in W.iterrows():
    L.append(f'| {l} | ' + ' | '.join(f'{row[r]:.2f}' for r in RECS) + ' | ' + ' | '.join(str(int(P.loc[l, r])) for r in RECS) + f' | {row.posto_medio:.1f} |')
ns = int(t.n_sementes.max())
L += ['', f'Sementes por complexo: até {ns}. ' + ('Meia-amplitude mediana entre sementes: '
      f'{t.meia_amplitude.median():.2f} kcal/mol; diferenças menores que isso entre dois ligantes não separam postos.' if ns > 1 else
      'Com 1 semente não há estimativa de incerteza; o ruído observado entre execuções idênticas em Ag/Bov (DP do Δ ≈ 0,14 log10 ≈ 0,19 kcal/mol) é o piso a considerar.'), '',
     '**Referências (fora do ranking):** ' + '; '.join(
         f'{x} ' + '/'.join(f'{t[(t.receptor == r) & (t.ligand == x)].dG.mean():.2f}' for r in RECS) for x in REF) + ' (Ag/Sf/Bov, kcal/mol).']
open(os.path.join(RES, 'ranking_energia.md'), 'w', encoding='utf-8').write('\n'.join(L) + '\n')
print('\n'.join(L))
