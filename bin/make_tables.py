"""Tabelas em markdown (só códigos) para o documento comparativo. Saída: results/tabelas.md"""
import os
import numpy as np, pandas as pd
from scipy import stats

RES = os.path.join(os.path.dirname(__file__), '..', 'results')
S = pd.read_csv(f'{RES}/complex_summary.csv'); A = pd.read_csv(f'{RES}/affinity_por_complexo.csv')
cl = pd.read_csv(os.path.join(RES, '..', 'data', 'cluspro.tsv'), sep='\t').set_index('code')
T = ['L%02d' % i for i in range(1, 12)] + ['G1', 'G2']
NEG = ['N1', 'N2', 'N3', 'N4']
P = S[S.arm == 'P']; Ss = S[S.arm == 'S']
out = []
def md(df, fmt=None):
    cols = list(df.columns); s = '| ' + ' | '.join(cols) + ' |\n|' + '|'.join(['---'] * len(cols)) + '|\n'
    for _, r in df.iterrows():
        s += '| ' + ' | '.join(('' if (isinstance(v, float) and np.isnan(v)) else (f'{v:.2f}' if isinstance(v, float) else str(v))) for v in r.values) + ' |\n'
    return s

# Tabela A: por ligante, 3 receptores nativos
rows = []
for l in T + NEG:
    r = {'Ligante': l}
    for rec in ['Ag', 'Sf', 'Bov']:
        p = P[(P.receptor == rec) & (P.ligand == l)]; a = A[(A.receptor == rec) & (A.ligand == l)]
        r[f'ipSAE {rec}'] = float(p.ipsae.iloc[0]) if len(p) else np.nan
    for rec in ['Ag', 'Sf', 'Bov']:
        a = A[(A.receptor == rec) & (A.ligand == l)]
        r[f'log10 IC50 {rec}'] = float(a.y.iloc[0]) if len(a) else np.nan
    for rec in ['Ag', 'Sf']:
        a = A[(A.receptor == rec) & (A.ligand == l)]
        r[f'p(lig) {rec}'] = float(a.pb.iloc[0]) if len(a) else np.nan
    r['ClusPro LE'] = float(cl.lowest_energy[l]) if l in cl.index else np.nan
    rows.append(r)
ta = pd.DataFrame(rows)
out.append('### Tabela A. Por ligante (receptores nativos)\n\n' + md(ta))

# Tabela B: resumo por receptor
rb = []
for rec in ['Ag', 'Sf', 'Bov']:
    p = P[P.receptor == rec].set_index('ligand'); a = A[A.receptor == rec].set_index('ligand')
    t = p.loc[T]; n = p.loc[NEG]
    rb.append({'Receptor': rec, 'ipTM alvos': t.iptm.mean(), 'ipTM negativos': n.iptm.mean(), 'ipSAE alvos': t.ipsae.mean(), 'ipSAE negativos': n.ipsae.mean(),
               'log10 IC50 alvos': a.loc[T].y.mean(), 'log10 IC50 negativos': a.loc[NEG].y.mean(), 'p(ligante) alvos': a.loc[T].pb.mean(), 'p(ligante) negativos': a.loc[NEG].pb.mean(),
               'benzamidina log10 IC50': float(a.loc['BEN', 'y']) if 'BEN' in a.index else np.nan})
out.append('### Tabela B. Resumo por receptor (média dos 13 alvos e dos 4 negativos)\n\n' + md(pd.DataFrame(rb)))

# Tabela C: magnitude absoluta vs Ki experimental (G2 e L01)
rc = []
ki = {'G2': 100.0, 'L01': 4000.0}  # uM; Ki experimental informado pelo usuario (0,10 mM e 4,00 mM)
for l, k in ki.items():
    for rec in ['Ag', 'Sf', 'Bov']:
        y = float(A[(A.receptor == rec) & (A.ligand == l)].y.iloc[0]); ic = 10 ** y
        rc.append({'Ligante': l, 'Receptor': rec, 'IC50 previsto (µM)': ic, 'Ki experimental (µM)': k, 'previsto mais forte que o experimental (vezes)': k / ic})
out.append('### Tabela C. Magnitude prevista × Ki experimental conhecido\n\n' + md(pd.DataFrame(rc)))

# Tabela D: concordância entre braços P e S (estrutura: ocupação do S1 pelo P1)
rd = []
# ligantes com Lys/Arg: lidos do arquivo privado (sequências não são versionadas)
kr = [l.split('	')[0] for l in open(os.path.join(RES, '..', 'data', 'private', 'ligands.tsv')).read().splitlines()[1:] if any(a in l.split('	')[1] for a in 'KR') and l.split('	')[0][0] in 'LG']
for rec in ['Ag', 'Sf', 'Bov', 'Ag_D189G', 'Sf_D189G']:
    p = P[(P.receptor == rec) & P.ligand.isin(kr)]; s = Ss[(Ss.receptor == rec) & Ss.ligand.isin(kr)]
    rd.append({'Receptor': rec, 'P1 no S1, braço P': p.s1_pocket_frac.mean(), 'P1 no S1, braço S (SMILES)': s.s1_pocket_frac.mean(),
               'ipTM braço P': p.iptm.mean(), 'ipTM braço S': s.iptm.mean()})
out.append('### Tabela D. Braço P × braço S (10 peptídeos com Lys/Arg)\n\n' + md(pd.DataFrame(rd)))
open(f'{RES}/tabelas.md', 'w', encoding='utf-8').write('\n\n'.join(out))
print('\n\n'.join(out[1:]))
