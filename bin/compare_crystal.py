"""Controles cristalográficos (pré-corte do Boltz-2): X1 (PRN, 5XW8), X2 (PRY, 5XW9), X3 (TRE, 5XWJ) em tripsina suína.
Compara cada pose predita com o cristal: LRMSD do peptídeo (átomos pesados pareados por nome, cap acetil excluído),
fnat (contatos resíduo-resíduo a <= 4,5 A) e ocupação do S1 (P1 vs Asp189) no cristal e na predição."""
import os, re, sys, glob
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from gcommon import *

PAIRS = {'X1': ('5XW8', 'C'), 'X2': ('5XW9', 'C'), 'X3': ('5XWJ', 'C')}
CUT = 4.5


def res_dict(chains):
    """Cristal: numeração do PDB não é a de tripsina (cadeias fragmentadas); mapeia por alinhamento ao 1BTY."""
    res = [r for ch in chains for r in ch if r.id[0] == ' ']
    seq = ''.join(THREE.get(r.get_resname(), 'X') for r in res)
    nm = number_map(seq)
    return {nm[i]: r for i, r in enumerate(res) if i in nm}


def contacts(pep_res, rec_by_num):
    out = set()
    for i, pr in enumerate(pep_res):
        for n, rr in rec_by_num.items():
            if min(np.linalg.norm(a.coord - b.coord) for a in pr if a.element != 'H' for b in rr if b.element != 'H') <= CUT:
                out.add((i, n))
    return out


def p1_dist(pep_res, rec_by_num, seq):
    ks = [i for i, c in enumerate(seq) if c in 'KR']
    if not ks or 189 not in rec_by_num:
        return np.nan
    r = pep_res[ks[-1]]; d = rec_by_num[189]
    names = ['NZ'] if r.get_resname() == 'LYS' else ['NE', 'CZ', 'NH1', 'NH2']
    return float(min(np.linalg.norm(r[n].coord - d[o].coord) for n in names if n in r for o in ('OD1', 'OD2')))


rows = []
crystal = {}
for x, (pdbid, pch) in PAIRS.items():
    m = load(f'{ROOT}/refs/{pdbid}.pdb')
    rec = res_dict([m['A']] + ([m['B']] if 'B' in [c.id for c in m] and pdbid != '5XWJ' else []))
    pep = [r for r in m[pch] if r.id[0] == ' ' and r.get_resname() != 'ACE']
    crystal[x] = (rec, pep)
    assert rec.get(195) is not None and rec[195].get_resname() == 'SER', f'{pdbid}: numeracao 195 nao e Ser'

for pdb in sorted(glob.glob(f'{ROOT}/runs/P_s*/boltz_results_P/predictions/Pig__X*/Pig__X*_model_*.pdb')):
    x = re.search(r'Pig__(X\d)', pdb).group(1); seed = re.search(r'P_s(\d+)', pdb).group(1); k = int(re.search(r'_model_(\d+)', pdb).group(1))
    rec_c, pep_c = crystal[x]
    m = load(pdb); A = receptor_residues(m, 'A'); B = [r for r in m['B'] if r.id[0] == ' ']
    nmap = number_map(chain_seq(m['A'])); rn = {nmap[i]: r for i, r in enumerate(A) if i in nmap}
    common = sorted(set(rn) & set(rec_c))
    P = np.array([rn[n]['CA'].coord for n in common]); Q = np.array([rec_c[n]['CA'].coord for n in common])
    R, t = kabsch(P, Q)
    ca_rmsd = rmsd(P @ R.T + t, Q)
    # LRMSD do peptídeo: átomos pesados com mesmo nome resíduo a resíduo
    pa, qa = [], []
    for rp, rc in zip(B, pep_c):
        for a in rp:
            if a.element != 'H' and a.get_id() in rc:
                pa.append(a.coord); qa.append(rc[a.get_id()].coord)
    lrmsd = rmsd(np.array(pa) @ R.T + t, np.array(qa))
    # fnat
    for r in B:
        for a in r:
            a.coord = a.coord @ R.T + t
    for r in A:
        for a in r:
            a.coord = a.coord @ R.T + t
    nat = contacts(pep_c, {n: rec_c[n] for n in rec_c}); pre = contacts(B, rn)
    fnat = len(nat & pre) / len(nat) if nat else np.nan
    seq = ''.join(THREE[r.get_resname()] for r in pep_c)
    conf = json.load(open(pdb.replace(f'Pig__{x}_model_{k}.pdb', f'confidence_Pig__{x}_model_{k}.json')))
    rows.append(dict(control=x, seed=seed, model=k, ca_rmsd_A=ca_rmsd, lrmsd_A=lrmsd, fnat=fnat,
                     p1_asp189_cristal_A=p1_dist(pep_c, rec_c, seq), p1_asp189_pred_A=p1_dist(B, rn, seq),
                     iptm=conf['iptm'], confidence=conf['confidence_score']))
df = pd.DataFrame(rows)
os.makedirs(f'{ROOT}/results', exist_ok=True)
df.to_csv(f'{ROOT}/results/crystal_controls.csv', index=False)
print(df.groupby('control').agg(n=('model', 'count'), lrmsd_mean=('lrmsd_A', 'mean'), lrmsd_min=('lrmsd_A', 'min'),
                                fnat_mean=('fnat', 'mean'), fnat_max=('fnat', 'max'), p1_pred=('p1_asp189_pred_A', 'mean'),
                                p1_cristal=('p1_asp189_cristal_A', 'mean'), iptm=('iptm', 'mean')).round(3))
