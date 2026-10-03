"""Mapa de resíduos de interação peptídeo-tripsina a < 5,0 Å. Roda no servidor.

Definição (única): resíduo do receptor está "em contato" numa pose se QUALQUER átomo pesado seu fica a < 5,0 Å de QUALQUER
átomo pesado do ligante. Numeração = tripsina (1BTY), via alinhamento (gcommon.number_map). O corte de QC de 4,5 Å
(qc_thresholds.yaml) NÃO é alterado; este é um mapa descritivo separado.
Consenso por complexo: freq = fração das poses (que passam no QC) com o resíduo a < 5,0 Å; "contato de consenso" = freq >= 0,5.
Saída só com número/nome do resíduo do RECEPTOR e posição do peptídeo (P1..Pn do N-terminal, só no braço P); nunca a sequência.

Modos:
  python bin/contatos_5A.py               braço P (peptídeo como cadeia), poses do QC validado (poses.csv, pose_qc); todas as sementes
  python bin/contatos_5A.py --set C10     conjunto com ligante carregado (runs/C10_s*): pose do braço S; QC leve (sem choque < 2,2 Å,
                                          tríade His57-Asp102 <= 3,3 Å e His57-Ser195 <= 3,8 Å, os mesmos limites de qc_thresholds.yaml)
Saídas: results/contatos_5A[_<SET>].csv (longo) e results/contatos_5A[_<SET>]_<rec>.csv (matriz resíduo x ligante, freq >= 0,5 em algum).
"""
import os, re, glob, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from gcommon import *

CUT = 5.0
RES = f'{ROOT}/results'
SET = sys.argv[2] if len(sys.argv) > 2 and sys.argv[1] == '--set' else None
SUF = f'_{SET}' if SET else ''
nmaps = {}


def numbered(m, rec):
    A = receptor_residues(m, 'A')
    if rec not in nmaps:
        nmaps[rec] = number_map(chain_seq(m['A']))
    nm = nmaps[rec]
    return A, nm, {nm[i]: r for i, r in enumerate(A) if i in nm}


def light_qc(A_by_num, Bc, Ac_all):
    D = np.linalg.norm(Bc[:, None, :] - Ac_all[None, :, :], axis=2)
    if (D < 2.2).any():
        return False
    try:
        h, d, s = A_by_num[57], A_by_num[102], A_by_num[195]
        hd = min(dist(h['ND1'], d[a]) for a in ('OD1', 'OD2')); hs = dist(h['NE2'], s['OG'])
    except KeyError:
        return False
    return hd <= 3.3 and hs <= 3.8


def pose_contacts(pdb, rec, arm, check_qc):
    m = load(pdb)
    A, nm, by_num = numbered(m, rec)
    if arm == 'P':
        B = [x for x in m['B'] if x.id[0] == ' ']
        Bat = [(a.coord, i + 1) for i, x in enumerate(B) for a in x if a.element != 'H']
    else:
        Bat = [(a.coord, 0) for x in m['B'] for a in x if a.element != 'H']
    Bc = np.array([c for c, _ in Bat]); Bp = np.array([p for _, p in Bat])
    if check_qc:
        Ac_all = np.array([a.coord for r in A for a in r if a.element != 'H'])
        if not light_qc(by_num, Bc, Ac_all):
            return None
    d = {}
    for i, res in enumerate(A):
        if i not in nm:
            continue
        Ac = np.array([a.coord for a in res if a.element != 'H'])
        D = np.linalg.norm(Bc[:, None, :] - Ac[None, :, :], axis=2)
        dm = float(D.min())
        if dm < CUT:
            d[nm[i]] = (res.get_resname(), dm, int(Bp[np.unravel_index(D.argmin(), D.shape)[0]]))
    return d


# lista de poses a processar: (rec, lig, seed, pdb, arm, check_qc)
jobs = []
if SET is None:
    P = pd.read_csv(f'{RES}/poses.csv'); P = P[(P.arm == 'P') & (P.pose_qc == True)]  # noqa: E712
    for r in P.itertuples():
        fs = glob.glob(f'{ROOT}/runs/P_s{int(r.seed)}/boltz_results_*/predictions/{r.receptor}__{r.ligand}/{r.receptor}__{r.ligand}_model_{int(r.model)}.pdb')
        if fs:
            jobs.append((r.receptor, r.ligand, int(r.seed), fs[0], 'P', False))
else:
    for f in sorted(glob.glob(f'{ROOT}/runs/{SET}_s*/boltz_results_*/predictions/*/*_model_*.pdb')):
        nm = os.path.basename(os.path.dirname(f)); rec, lig = nm.split('__')
        jobs.append((rec, lig, int(re.search(rf'/{SET}_s(\d+)/', f).group(1)), f, 'S', True))

per = {}
n_total = n_fail = 0
for rec, lig, seed, pdb, arm, chk in jobs:
    n_total += 1
    d = pose_contacts(pdb, rec, arm, chk)
    if d is None:
        n_fail += 1; continue
    per.setdefault((rec, lig), []).append(d)

rows = []
for (rec, lig), poses in per.items():
    n = len(poses)
    for resnum in sorted({k for d in poses for k in d}):
        hit = [d[resnum] for d in poses if resnum in d]
        pos = [h[2] for h in hit]
        rows.append(dict(receptor=rec, ligand=lig, resnum=resnum, resname=hit[0][0], n_poses=n, n_contact=len(hit),
                         freq=len(hit) / n, dmin_mediana_A=float(np.median([h[1] for h in hit])),
                         pos_pep_mais_proxima=(int(pd.Series(pos).mode().iloc[0]) if SET is None else np.nan)))
df = pd.DataFrame(rows); df.to_csv(f'{RES}/contatos_5A{SUF}.csv', index=False)
for rec, g in df.groupby('receptor'):
    keep = g[g.freq >= 0.5].resnum.unique()
    M = g[g.resnum.isin(keep)].pivot_table(index='resnum', columns='ligand', values='freq').fillna(0.0)
    names = g.drop_duplicates('resnum').set_index('resnum').resname
    M.insert(0, 'resname', [names.get(i, '') for i in M.index])
    M.to_csv(f'{RES}/contatos_5A{SUF}_{rec}.csv')
print(f'modo={SET or "P (QC validado)"} | poses lidas={n_total} reprovadas no QC={n_fail} | complexos com mapa={len(per)} | linhas={len(df)}')
cons = df[df.freq >= 0.5].groupby(['receptor', 'ligand']).size().groupby('receptor').describe()[['count', 'mean', 'min', 'max']]
print('resíduos de consenso (freq>=0,5) por complexo, por receptor:\n', cons.round(1).to_string())
