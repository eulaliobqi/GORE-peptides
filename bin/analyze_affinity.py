"""Testes da afinidade do Boltz-2 (braço S). Lê results/affinity.csv; escreve results/affinity_tests.txt.
Testes: (1) alvos x negativos; (2) ordem G2 > L01 (Ki 0,10 mM vs 4,00 mM, único controle de ordem); (3) efeito do D189G;
(4) concordância entre os 2 modelos do ensemble; (5) afinidade x ipSAE; (6) afinidade x ClusPro; (7) benzamidina e controles cristalográficos."""
import os, sys, io
import numpy as np, pandas as pd
from scipy import stats

RES = os.path.join(os.path.dirname(__file__), '..', 'results')
aff = pd.read_csv(f'{RES}/affinity.csv')
summ = pd.read_csv(f'{RES}/complex_summary.csv'); P = summ[summ.arm == 'P']
T = ['L%02d' % i for i in range(1, 12)] + ['G1', 'G2']
NEG = ['N1', 'N2', 'N3', 'N4']
out = io.StringIO()
def w(*a):
    print(*a); print(*a, file=out)

A = aff.groupby(['receptor', 'ligand']).agg(y=('log10_ic50_uM', 'mean'), pb=('p_binder', 'mean'),
                                              y1=('log10_ic50_uM_1', 'mean'), y2=('log10_ic50_uM_2', 'mean')).reset_index()
A.to_csv(f'{RES}/affinity_por_complexo.csv', index=False)
w('n complexos com afinidade:', len(A))
def boot(x, n=3000):
    rs = np.random.RandomState(1)
    return np.percentile([rs.choice(x, len(x)).mean() for _ in range(n)], [2.5, 97.5])
def auc(a, b):  # prob. de a > b
    return stats.mannwhitneyu(a, b, alternative='greater').statistic / (len(a) * len(b))

w('\n(1) ALVOS x NEGATIVOS (y = log10 IC50 em uM; menor = mais forte; pb = probabilidade de ligante)')
for r in ['Ag', 'Sf', 'Bov']:
    d = A[A.receptor == r].set_index('ligand')
    t = d.loc[[x for x in T if x in d.index]]; n = d.loc[[x for x in NEG if x in d.index]]
    if len(t) and len(n):
        w(f'  {r}: y alvos {t.y.mean():.2f}+-{t.y.std():.2f} | negativos {n.y.mean():.2f} | AUC(y menor em alvos) {auc(-t.y, -n.y):.2f} | pb alvos {t.pb.mean():.2f} negativos {n.pb.mean():.2f} | AUC(pb) {auc(t.pb, n.pb):.2f}')
    if 'BEN' in d.index:
        w(f'     benzamidina: y={d.loc["BEN","y"]:.2f} pb={d.loc["BEN","pb"]:.2f}')

w('\n(2) ORDEM G2 > L01 (Ki experimental 0,10 mM vs 4,00 mM; previsto correto se y(G2) < y(L01))')
for r in ['Ag', 'Sf', 'Bov']:
    d = A[A.receptor == r].set_index('ligand')
    if {'G2', 'L01'} <= set(d.index):
        w(f'  {r}: y(G2)={d.loc["G2","y"]:.2f} y(L01)={d.loc["L01","y"]:.2f} dif={d.loc["G2","y"]-d.loc["L01","y"]:+.2f} -> {"CORRETO" if d.loc["G2","y"] < d.loc["L01","y"] else "INVERTIDO"} | pb G2={d.loc["G2","pb"]:.2f} L01={d.loc["L01","pb"]:.2f}')
    if 'G1' in d.index:
        w(f'     G1: y={d.loc["G1","y"]:.2f} (Ki experimental de G1 ainda nao extraido)')

w('\n(3) EFEITO DO D189G (13 alvos, pareado nativo -> mutante; esperado: y sobe = afinidade cai)')
for r in ['Ag', 'Sf']:
    a = A[A.receptor == r].set_index('ligand'); b = A[A.receptor == r + '_D189G'].set_index('ligand')
    idx = [x for x in T if x in a.index and x in b.index]
    if len(idx) >= 5:
        dl = (b.loc[idx, 'y'] - a.loc[idx, 'y']).values
        wil = stats.wilcoxon(dl)
        w(f'  {r}: n={len(idx)} dy medio {dl.mean():+.2f} (IC95 boot {boot(dl).round(2)}) Wilcoxon p={wil.pvalue:.4f}; pb {a.loc[idx,"pb"].mean():.2f} -> {b.loc[idx,"pb"].mean():.2f}')

w('\n(4) CONCORDANCIA ENTRE OS 2 MODELOS DO ENSEMBLE (y1 vs y2, todos os complexos)')
ok = A.dropna(subset=['y1', 'y2'])
if len(ok) > 3:
    rho, p = stats.spearmanr(ok.y1, ok.y2); w(f'  Spearman {rho:.2f} (p={p:.1e}, n={len(ok)}); media |y1-y2| = {(ok.y1-ok.y2).abs().mean():.2f}')

w('\n(5) AFINIDADE x ipSAE (alvos, por receptor)')
for r in ['Ag', 'Sf', 'Bov']:
    d = A[(A.receptor == r) & A.ligand.isin(T)].set_index('ligand'); p_ = P[P.receptor == r].set_index('ligand')
    idx = [x for x in d.index if x in p_.index]
    if len(idx) >= 5:
        rho, p = stats.spearmanr(d.loc[idx, 'y'], p_.loc[idx, 'ipsae'])
        w(f'  {r}: Spearman(y, ipSAE) = {rho:+.2f} (p={p:.3f}, n={len(idx)}); negativo = afinidade prevista mais forte com ipSAE maior')

w('\n(6) AFINIDADE x ClusPro (energia mais negativa = melhor; positivo = concorda)')
c = pd.read_csv(os.path.join(RES, '..', 'data', 'cluspro.tsv'), sep='\t').set_index('code')
for r in ['Ag', 'Sf', 'Bov']:
    d = A[A.receptor == r].set_index('ligand'); idx = [x for x in c.index if x in d.index]
    if len(idx) >= 5:
        for cm in ['score_center', 'lowest_energy']:
            rho, p = stats.spearmanr(-d.loc[idx, 'y'], -c.loc[idx, cm])
            w(f'  {r} {cm}: rho={rho:+.2f} (p={p:.3f}, n={len(idx)})')

w('\n(7) RANKING DE AFINIDADE DOS ALVOS (y; menor = mais forte)')
for r in ['Ag', 'Sf', 'Bov']:
    d = A[(A.receptor == r) & A.ligand.isin(T)].sort_values('y')
    w(f'  {r}: ' + ' '.join(f'{l}:{v:.2f}' for l, v in zip(d.ligand, d.y)))
d = A[A.receptor == 'Pig']
if len(d):
    w('  Pig (cristalograficos + benzamidina): ' + ' '.join(f'{l}:{v:.2f}/pb{p:.2f}' for l, v, p in zip(d.ligand, d.y, d.pb)))

w('\n(8) VIES DE TAMANHO: afinidade prevista x comprimento do peptideo (so o comprimento entra no texto; sequencias ficam em data/private)')
lig = pd.read_csv(os.path.join(RES, '..', 'data', 'private', 'ligands.tsv'), sep='\t').set_index('code')
L = lig.seq.str.len()
cl = pd.read_csv(os.path.join(RES, '..', 'data', 'cluspro.tsv'), sep='\t').set_index('code')
for r in ['Ag', 'Sf', 'Bov']:
    d = A[(A.receptor == r) & A.ligand.isin(T)].set_index('ligand')
    if len(d) >= 6:
        rho, p = stats.spearmanr(d.y, L[d.index]); w(f'  {r}: Spearman(y, comprimento) = {rho:+.2f} (p={p:.3f}, n={len(d)}); negativo = peptideo maior tem afinidade prevista mais forte')
        pp = P[P.receptor == r].set_index('ligand').loc[[x for x in d.index if x in P[P.receptor == r].ligand.values]]
        rho2, p2 = stats.spearmanr(pp.ipsae, L[pp.index]); w(f'     Spearman(ipSAE, comprimento) = {rho2:+.2f} (p={p2:.3f})')
        # correlacao parcial por postos: remove o efeito do comprimento nos dois lados
        from scipy.stats import rankdata
        def resid(x, z):
            rx, rz = rankdata(x), rankdata(z); b = np.polyfit(rz, rx, 1); return rx - np.polyval(b, rz)
        idx = [x for x in cl.index if x in d.index]
        for cm in ['lowest_energy']:
            a1 = resid(-d.loc[idx, 'y'].values, L[idx].values); a2 = resid(-cl.loc[idx, cm].values, L[idx].values)
            rr, pr = stats.pearsonr(a1, a2)
            w(f'     ClusPro x afinidade, bruto rho={stats.spearmanr(-d.loc[idx,"y"], -cl.loc[idx,cm])[0]:+.2f}; parcial controlando comprimento r={rr:+.2f} (p={pr:.3f}, n={len(idx)})')
cl_rho, cl_p = stats.spearmanr(-cl.lowest_energy, L[cl.index]); w(f'  ClusPro (lowest energy) x comprimento: rho={cl_rho:+.2f} (p={cl_p:.3f}, n={len(cl)}); positivo = peptideo maior tem energia mais favoravel')
open(f'{RES}/affinity_tests.txt', 'w', encoding='utf-8').write(out.getvalue())
