"""QC dos receptores dobrados pelo Boltz-2 (antes de qualquer docking): confiança, tríade, dissulfetos,
concordância com o cristal (1BTY), estereoquímica (MolProbity ramalyze). Limiares em qc_thresholds.yaml (pré-registrados)."""
import os, re, sys, glob, subprocess, json
import numpy as np, pandas as pd, yaml
sys.path.insert(0, os.path.dirname(__file__))
from gcommon import *

TM = os.path.expanduser('~/miniforge3/envs/structure/bin/TMalign')
RAMA = os.path.expanduser('~/miniforge3/envs/molprobity/bin/molprobity.ramalyze')
thr = yaml.safe_load(open(f'{ROOT}/qc_thresholds.yaml'))['receptor']
refm, rseq, rnum, rres = ref_1bty()
ref_by_num = {r.id[1]: r for r in rres}


def by_num(residues, nmap):
    return {nmap[i]: r for i, r in enumerate(residues) if i in nmap}


def ramalyze(path):
    out = subprocess.run([RAMA, path], capture_output=True, text=True).stdout
    fav = re.search(r'([\d.]+)% favored', out); outl = re.search(r'([\d.]+)% outliers', out)
    return (float(fav.group(1)) if fav else np.nan, float(outl.group(1)) if outl else np.nan)


def tmscore(model, ref):
    out = subprocess.run([TM, model, ref], capture_output=True, text=True).stdout
    m = re.findall(r'TM-score=\s*([\d.]+)', out)
    r = re.search(r'RMSD=\s*([\d.]+)', out)
    return (float(m[1]) if len(m) > 1 else np.nan, float(r.group(1)) if r else np.nan)  # normalizado pelo comprimento da referência


rows = []
for pdb in sorted(glob.glob(f'{ROOT}/runs/recqc_s*/boltz_results_receptors_msa/predictions/*/*_model_*.pdb')):
    seed = re.search(r'recqc_s(\d+)', pdb).group(1)
    name = os.path.basename(os.path.dirname(pdb)); k = int(re.search(r'_model_(\d+)', pdb).group(1))
    m = load(pdb); res = receptor_residues(m, 'A'); seq = chain_seq(m['A'])
    nmap = number_map(seq); rn = by_num(res, nmap)
    plddt = np.array([r['CA'].bfactor for r in res])
    ser = rn.get(195)
    pocket = [r for r in res if ser is not None and min(np.linalg.norm(a.coord - ser['OG'].coord) for a in r) <= 10.0]
    ppl = float(np.mean([r['CA'].bfactor for r in pocket])) if pocket else np.nan
    ident = {n: (rn[n].get_resname() if n in rn else None) for n in (57, 102, 189, 195)}
    d_hd = min(dist(rn[57]['ND1'], rn[102][a]) for a in ('OD1', 'OD2')) if 57 in rn and 102 in rn else np.nan
    d_hs = dist(rn[57]['NE2'], rn[195]['OG']) if 57 in rn and 195 in rn else np.nan
    cys = [r for r in res if r.get_resname() == 'CYS']
    paired = set(); ss = 0
    for i, a in enumerate(cys):
        for b in cys[i + 1:]:
            if 1.8 <= dist(a['SG'], b['SG']) <= 2.3:
                paired.update([a.id[1], b.id[1]]); ss += 1
    # alinhamento com 1BTY por numeração de tripsina
    common = sorted(set(rn) & set(ref_by_num))
    P = np.array([rn[n]['CA'].coord for n in common]); Q = np.array([ref_by_num[n]['CA'].coord for n in common])
    R, t = kabsch(P, Q)
    s1 = [n for n in S1_SET if n in common]
    s1rmsd = rmsd(np.array([rn[n]['CA'].coord for n in s1]) @ R.T + t, np.array([ref_by_num[n]['CA'].coord for n in s1]))
    tm, tmr = tmscore(pdb, f'{ROOT}/refs/1BTY.pdb')
    fav, outl = ramalyze(pdb) if k == 0 else (np.nan, np.nan)
    rows.append(dict(receptor=name, seed=seed, model=k, n_res=len(res), mean_plddt=plddt.mean(), pocket_plddt=ppl,
                     res57=ident[57], res102=ident[102], res189=ident[189], res195=ident[195],
                     his_asp_A=d_hd, his_ser_A=d_hs, n_cys=len(cys), n_ss=ss, all_cys_paired=(ss == len(cys) // 2),
                     s1_ca_rmsd_A=s1rmsd, tm_vs_1BTY=tm, rmsd_tmalign=tmr, rama_fav=fav, rama_out=outl))
df = pd.DataFrame(rows)
df['mean_plddt'] = df['mean_plddt']  # escala 0-100 (B-factor)
def ok(r):
    c = [r.mean_plddt >= thr['mean_plddt_min'], r.pocket_plddt >= thr['pocket_plddt_min'],
         r.tm_vs_1BTY >= thr['tm_score_vs_1BTY_min'], r.s1_ca_rmsd_A <= thr['s1_ca_rmsd_max_A'],
         r.his_asp_A <= thr['his57_ND1_asp102_OD_max_A'], r.his_ser_A <= thr['his57_NE2_ser195_OG_max_A'],
         bool(r.all_cys_paired)]
    if not np.isnan(r.rama_fav):
        c += [r.rama_fav >= thr['rama_favored_min_pct'], r.rama_out <= thr['rama_outliers_max_pct']]
    return all(c)
df['qc_pass'] = df.apply(ok, axis=1)
os.makedirs(f'{ROOT}/results', exist_ok=True)
df.to_csv(f'{ROOT}/results/qc_receptors.csv', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
print(df.groupby('receptor').agg(mean_plddt=('mean_plddt', 'mean'), pocket=('pocket_plddt', 'mean'), tm=('tm_vs_1BTY', 'mean'),
                                 s1rmsd=('s1_ca_rmsd_A', 'mean'), hisasp=('his_asp_A', 'mean'), hisser=('his_ser_A', 'mean'),
                                 ss=('n_ss', 'mean'), qc_pass=('qc_pass', 'mean'), n=('model', 'count')).round(3))
