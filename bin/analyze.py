"""Análise dos complexos Boltz-2 (braço P = peptídeo como cadeia; braço S = ligante SMILES + afinidade).
Por pose: confiança, ipSAE/LIS, contatos (S1, catalíticos, S'2), ocupação do S1 pelo P1, geometria da cadeia,
QC da pose (clashes, quiralidade, ômega, tríade do receptor). Por complexo: agregação entre poses/sementes e clusterização.
Saídas: results/poses.csv, results/complex_summary.csv, results/affinity.csv
"""
import os, re, sys, glob, json
import numpy as np, pandas as pd, yaml
sys.path.insert(0, os.path.dirname(__file__))
from gcommon import *
from Bio.PDB.vectors import calc_dihedral

thr = yaml.safe_load(open(f'{ROOT}/qc_thresholds.yaml'))
CUT = thr['contact']['cutoff_A']; S1OCC = thr['contact']['s1_occupied_A']; CLU = thr['cluster']['pose_rmsd_A']
rt = thr['receptor']
lig = {}
for l in open(f'{ROOT}/data/private/ligands.tsv').read().splitlines()[1:]:
    c = l.split('\t'); lig[c[0]] = c[1]

try:
    from rdkit import Chem
    from pep_smiles import peptide_smiles
except Exception:
    Chem = None


def basic_atoms_smiles(smi):
    """índices (na ordem do SMILES) dos N/C básicos do P1 (Lys NZ; Arg N,N,N,C guanidínio)."""
    m = Chem.MolFromSmiles(smi)
    idx = []
    for patt, pos in (('[NH3+][CH2][CH2][CH2][CH2]', [0]), ('[NX3;H1]([CH2])C(=[NH2+])[NH2]', [0, 2, 3, 1])):
        q = Chem.MolFromSmarts(patt)
        for mt in m.GetSubstructMatches(q):
            idx += [mt[i] for i in pos] if patt.startswith('[NH3+]') else list(mt[:4])
    return idx, [a.GetSymbol() for a in m.GetAtoms()]


def chir_sign(res):
    try:
        n, ca, c, cb = res['N'].coord, res['CA'].coord, res['C'].coord, res['CB'].coord
    except KeyError:
        return None
    return float(np.sign(np.dot(np.cross(n - ca, c - ca), cb - ca)))


def analyze_pose(pdb, name, arm, seed, k, recseq_len):
    rec_name, lig_code = name.split('__')
    m = load(pdb)
    A = receptor_residues(m, 'A'); seqA = chain_seq(m['A']); nmap = number_map(seqA); rn = by_num = {nmap[i]: r for i, r in enumerate(A) if i in nmap}
    if arm == 'P':
        B = [r for r in m['B'] if r.id[0] == ' ']
        Batoms = [a for r in B for a in r if a.element != 'H']
    else:
        Batoms = [a for r in m['B'] for a in r if a.element != 'H']
        B = None
    Aatoms = [a for r in A for a in r if a.element != 'H']
    Acoords = np.array([a.coord for a in Aatoms]); Bcoords = np.array([a.coord for a in Batoms])
    D = np.linalg.norm(Bcoords[:, None, :] - Acoords[None, :, :], axis=2)
    out = dict(receptor=rec_name, ligand=lig_code, arm=arm, seed=seed, model=k)
    # contatos por resíduo (numeração de tripsina)
    rid_of = {}
    for j, a in enumerate(Aatoms):
        rid_of[j] = a.get_parent()
    contact_res = set()
    for i, j in zip(*np.where(D <= CUT)):
        contact_res.add(rid_of[j])
    inv = {r: n for n, r in rn.items()}
    cn = sorted(inv[r] for r in contact_res if r in inv)
    out.update(n_contact_res=len(cn), contacts_S1=sum(n in S1_SET for n in cn), contacts_cat=sum(n in CATAL for n in cn),
               contacts_oxy=sum(n in OXY for n in cn), contact_S2p=int(151 in cn), contact_H57=int(57 in cn),
               contact_S195=int(195 in cn), contact_D189=int(189 in cn), n_contacts_atoms=int((D <= CUT).sum()))
    out['contact_list'] = ';'.join(map(str, cn))
    out['min_dist_A'] = float(D.min())
    out['severe_clashes'] = int((D < 2.2).sum())
    # ocupação do S1 pelo P1 (Lys/Arg). Grupo básico do P1: Lys NZ; Arg NE/CZ/NH1/NH2 (braço P) ou átomos do SMILES (braço S).
    #  - p1_asp189_A: distância mínima ao carboxilato do Asp189 (só receptor nativo, que tem Asp189)
    #  - p1_ca189_A: distância mínima ao Cα da posição 189 (vale também no mutante D189G, sem cadeia lateral)
    d189 = rn.get(189)
    grp = []
    if arm == 'P':
        seq = lig[lig_code] if lig_code in lig else None
        if seq:
            ks = [i for i, c in enumerate(seq) if c in 'KR']
            if ks:
                r = B[ks[-1]]
                names = ['NZ'] if r.get_resname() == 'LYS' else ['NE', 'CZ', 'NH1', 'NH2']
                grp = [r[n].coord for n in names if n in r]
                out['p1_res'] = f"{seq[ks[-1]]}{ks[-1] + 1}"
    elif Chem is not None and lig_code in lig:
        smi = peptide_smiles(lig[lig_code])
        idx, elems = basic_atoms_smiles(smi)
        if idx and [a.element for a in Batoms] == elems:
            grp = [Batoms[i].coord for i in idx]
    p1_d = np.nan
    out['p1_ca189_A'] = float(min(np.linalg.norm(g - d189['CA'].coord) for g in grp)) if (grp and d189 is not None) else np.nan
    if grp and d189 is not None and d189.get_resname() == 'ASP':
        oxy = [d189['OD1'].coord, d189['OD2'].coord]
        p1_d = float(min(np.linalg.norm(g - o) for g in grp for o in oxy))
    # qualquer cadeia lateral do peptídeo no S1 (para peptídeos sem Lys/Arg): menor distância cadeia lateral -> Cα189
    if arm == 'P' and d189 is not None:
        best = (np.inf, None)
        for ri, r in enumerate(B):
            for at in r:
                if at.get_id() not in ('N', 'CA', 'C', 'O', 'OXT') and at.element != 'H':
                    dd = float(np.linalg.norm(at.coord - d189['CA'].coord))
                    if dd < best[0]:
                        best = (dd, ri + 1)
        out['sc_ca189_min_A'] = best[0]; out['sc_ca189_res'] = best[1]
    out['p1_asp189_A'] = p1_d
    out['s1_occupied'] = (p1_d <= S1OCC) if not np.isnan(p1_d) else np.nan
    # geometria da cadeia (braço P)
    if arm == 'P':
        omegas = []
        for a, b in zip(B, B[1:]):
            omegas.append(np.degrees(calc_dihedral(a['CA'].get_vector(), a['C'].get_vector(), b['N'].get_vector(), b['CA'].get_vector())))
        out['omega_min_abs'] = float(min(abs(o) for o in omegas)) if omegas else np.nan
        out['omega_trans'] = all(abs(o) >= thr['pose']['omega_trans_min_abs_deg'] for o in omegas) if omegas else True
        rs = [chir_sign(r) for r in A if r.get_resname() != 'GLY']; rs = [x for x in rs if x is not None]
        maj = np.sign(np.sum(rs)) if rs else 1.0
        ps = [chir_sign(r) for r in B if r.get_resname() != 'GLY']; ps = [x for x in ps if x is not None]
        out['chirality_ok'] = all(x == maj for x in ps)
        out['pep_plddt'] = float(np.mean([r['CA'].bfactor for r in B]))
    else:
        out['pep_plddt'] = float(np.mean([a.bfactor for a in Batoms]))
    # QC do receptor na pose
    ser = rn.get(195); his = rn.get(57); asp = rn.get(102)
    hd = min(dist(his['ND1'], asp[a]) for a in ('OD1', 'OD2')) if his is not None and asp is not None else np.nan
    hs = dist(his['NE2'], ser['OG']) if his is not None and ser is not None else np.nan
    out.update(his_asp_A=hd, his_ser_A=hs, rec_plddt=float(np.mean([r['CA'].bfactor for r in A])))
    out['triad_ok'] = bool(hd <= rt['his57_ND1_asp102_OD_max_A'] and hs <= rt['his57_NE2_ser195_OG_max_A'])
    # scissile: C do carbonila do P1 -> OG Ser195 (braço P; P1 no C-terminal usa C do carboxilato)
    if arm == 'P' and out.get('p1_res') and ser is not None:
        i = int(out['p1_res'][1:]) - 1
        out['p1CO_serOG_A'] = dist(B[i]['C'], ser['OG'])
    # coordenadas para clusterização: CA do receptor e átomos pesados do ligante
    caA = ca_coords(A)
    return out, caA, Bcoords


def kabsch_apply(P, Q, X):
    R, t = kabsch(P, Q)
    return X @ R.T + t


def load_conf(dirn, name, k):
    c = json.load(open(f'{dirn}/confidence_{name}_model_{k}.json'))
    pae = np.load(f'{dirn}/pae_{name}_model_{k}.npz')['pae']
    return c, pae


rows = []; coords = {}
runs = sorted(glob.glob(f'{ROOT}/runs/[PS]_s*/boltz_results_*/predictions/*'))
for d in runs:
    m = re.search(r'/([PS])_s(\d+)/', d); arm, seed = m.group(1), m.group(2)
    name = os.path.basename(d)
    for pdb in sorted(glob.glob(f'{d}/{name}_model_*.pdb')):
        k = int(re.search(r'_model_(\d+)\.pdb', pdb).group(1))
        try:
            row, caA, Bc = analyze_pose(pdb, name, arm, seed, k, None)
            c, pae = load_conf(d, name, k)
        except Exception as e:
            print('ERRO', pdb, repr(e)); continue
        nA = len(caA); nB = pae.shape[0] - nA
        ips, lis = ipsae_lis(pae, nA, nB, 10.0)
        pc = c.get('pair_chains_iptm', {})
        row.update(confidence=c['confidence_score'], iptm=c['iptm'], ptm=c['ptm'], lig_iptm=c.get('ligand_iptm'),
                   prot_iptm=c.get('protein_iptm'), complex_plddt=c['complex_plddt'], complex_iplddt=c['complex_iplddt'],
                   pair_iptm_AB=(pc.get('0', {}).get('1') if isinstance(pc.get('0'), dict) else None),
                   ipsae=ips, lis=lis)
        rows.append(row); coords[(name, arm, seed, k)] = (caA, Bc)
poses = pd.DataFrame(rows)
poses['pose_qc'] = (poses.severe_clashes <= thr['pose']['severe_clash_pairs_max']) & poses.triad_ok & \
                   (poses.rec_plddt >= rt['mean_plddt_min']) & \
                   (poses.get('chirality_ok', True) != False) & (poses.get('omega_trans', True) != False)
os.makedirs(f'{ROOT}/results', exist_ok=True)


# clusterização por complexo/braço (RMSD do ligante após alinhar Calfa do receptor à 1ª pose)
def cluster_stats(g):
    keys = [(r.receptor + '__' + r.ligand, r.arm, r.seed, r.model) for r in g.itertuples()]
    ref_ca, _ = coords[keys[0]]
    X = []
    for kk in keys:
        ca, Bc = coords[kk]
        X.append(kabsch_apply(ca, ref_ca, Bc))
    n = len(X); Dm = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            Dm[i, j] = Dm[j, i] = rmsd(X[i], X[j])
    un = set(range(n)); clusters = []
    # cluster guloso: centro = pose com mais vizinhos <= CLU
    while un:
        u = sorted(un)
        cnt = {i: sum(1 for j in u if Dm[i, j] <= CLU) for i in u}
        c = max(cnt, key=cnt.get)
        mem = [j for j in u if Dm[c, j] <= CLU]
        clusters.append(mem); un -= set(mem)
    clusters.sort(key=len, reverse=True)
    return dict(n_poses=n, n_clusters=len(clusters), top_cluster_frac=len(clusters[0]) / n,
                mean_pairwise_rmsd=float(Dm[np.triu_indices(n, 1)].mean()) if n > 1 else np.nan), clusters, g.index.tolist()


# limiar geométrico do S1 por Cα189, calibrado POR BRAÇO só nas poses do receptor nativo (Ag, Sf, Bov) com ponte salina P1-Asp189 (percentil 95)
CA_CUTS = {}
poses['s1_pocket'] = np.nan
for arm_ in ('P', 'S'):
    wt = poses[(poses.arm == arm_) & poses.receptor.isin(['Ag', 'Sf', 'Bov']) & (poses.s1_occupied == True) & poses.pose_qc]
    CA_CUTS[arm_] = float(np.percentile(wt.p1_ca189_A.dropna(), 95)) if len(wt) else np.nan
    m_ = (poses.arm == arm_) & poses.p1_ca189_A.notna()
    poses.loc[m_, 's1_pocket'] = (poses.loc[m_, 'p1_ca189_A'] <= CA_CUTS[arm_]).astype(float)
print('limiares p1-CA189 (percentil 95 por braço, poses nativas com ponte salina):', {k: round(v, 2) for k, v in CA_CUTS.items()}, 'A')
summ = []
for (rec, ligc, arm), g in poses.groupby(['receptor', 'ligand', 'arm']):
    good = g[g.pose_qc]
    cs, clusters, idx = cluster_stats(g.reset_index(drop=True))
    top = g.loc[g.confidence.idxmax()]
    def mean_(col, df=good):
        return float(df[col].mean()) if len(df) else np.nan
    s = dict(receptor=rec, ligand=ligc, arm=arm, n_qc_pass=len(good), frac_qc_pass=len(good) / len(g),
             iptm=mean_('iptm'), iptm_sd=float(good.iptm.std()) if len(good) > 1 else np.nan, ipsae=mean_('ipsae'), lis=mean_('lis'),
             confidence=mean_('confidence'), pep_plddt=mean_('pep_plddt'), complex_iplddt=mean_('complex_iplddt'),
             s1_frac=float(good.s1_occupied.astype(float).mean()) if len(good) and good.s1_occupied.notna().any() else np.nan,
             s1_pocket_frac=float(good.s1_pocket.astype(float).mean()) if len(good) and good.s1_pocket.notna().any() else np.nan,
             p1_ca189_A=mean_('p1_ca189_A'), sc_ca189_min_A=mean_('sc_ca189_min_A') if 'sc_ca189_min_A' in good else np.nan,
             p1_asp189_A=mean_('p1_asp189_A'), contacts_S1=mean_('contacts_S1'), contact_S2p_frac=mean_('contact_S2p'),
             contact_H57_frac=mean_('contact_H57'), p1CO_serOG_A=mean_('p1CO_serOG_A') if 'p1CO_serOG_A' in good else np.nan,
             best_model=f"s{top.seed}_m{top.model}", best_confidence=float(top.confidence), **cs)
    summ.append(s)
summ = pd.DataFrame(summ)

# afinidade (braço S): 1 arquivo por complexo/semente
arows = []
for f in glob.glob(f'{ROOT}/runs/S_s*/boltz_results_*/predictions/*/affinity_*.json'):
    seed = re.search(r'/S_s(\d+)/', f).group(1)
    name = os.path.basename(os.path.dirname(f)); rec, ligc = name.split('__')
    a = json.load(open(f))
    arows.append(dict(receptor=rec, ligand=ligc, seed=seed, log10_ic50_uM=a.get('affinity_pred_value'),
                      p_binder=a.get('affinity_probability_binary'),
                      log10_ic50_uM_1=a.get('affinity_pred_value1'), log10_ic50_uM_2=a.get('affinity_pred_value2'),
                      p_binder_1=a.get('affinity_probability_binary1'), p_binder_2=a.get('affinity_probability_binary2')))
aff = pd.DataFrame(arows)
if len(aff):
    aff['pIC50'] = 6 - aff['log10_ic50_uM']  # valor do Boltz-2 = log10(IC50 em uM)
    aff['pIC50_kcal_doc'] = aff['pIC50'] * 1.364  # conversao documentada pelo Boltz: (6 - y) * 1.364
aff.to_csv(f'{ROOT}/results/affinity.csv', index=False)
poses.to_csv(f'{ROOT}/results/poses.csv', index=False)
summ.to_csv(f'{ROOT}/results/complex_summary.csv', index=False)
print(len(poses), 'poses;', len(summ), 'complexos-braco;', len(aff), 'afinidades')
print(poses.groupby('arm').pose_qc.mean())
