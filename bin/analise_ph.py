"""Ranking por energia e mapa de resíduos <5 Å em DOIS estados de ionização do ligante (8,2 e 10), com ligante CARREGADO.
Roda no servidor (lê runs/ e results/). Uso: python bin/analise_ph.py [--a C82] [--b C10] [--neutral S]

EXPLORATÓRIO: o Boltz-2 padrão neutraliza o SMILES (standardize); aqui os conjuntos C* foram rodados com boltz_nostd (cargas
preservadas), entrada FORA da distribuição de treino do módulo de afinidade. A coluna 'neutro' (conjunto S) é o padrão do Boltz
e a referência validada. Só o ligante muda de estado de ionização; o receptor é o mesmo nos dois pH (o Boltz não titula).
Energia: ΔG = -1,364 * (6 - log10 IC50[uM]) kcal/mol (conversão do Boltz-2; predição de rede neural, não energia calculada).
Saídas: results/ranking_energia_ph.{csv,md}, results/contatos_5A_ph.{csv,md}, results/ph_ruido.txt
"""
import os, re, sys, glob, json, io
import numpy as np, pandas as pd
from scipy import stats

ROOT = os.path.expanduser('~/GOREs-boltz'); RES = f'{ROOT}/results'
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
A, B, N = arg('--a', 'C82'), arg('--b', 'C10'), arg('--neutral', 'S')
KT = 1.364
TARGETS = [f'L{i:02d}' for i in range(1, 12)] + ['G1', 'G2']; NEG = ['N1', 'N2', 'N3', 'N4']
RECS = ['Ag', 'Sf', 'Bov']
log = io.StringIO()
def P(*a): print(*a, file=log)


def load_aff(setname):
    rows = []
    for f in glob.glob(f'{ROOT}/runs/{setname}_s*/boltz_results_*/predictions/*/affinity_*.json'):
        seed = int(re.search(rf'/{setname}_s(\d+)/', f).group(1))
        rec, lig = os.path.basename(os.path.dirname(f)).split('__')
        y = json.load(open(f)).get('affinity_pred_value')
        rows.append(dict(receptor=rec, ligand=lig, seed=seed, dG=-KT * (6 - y)))
    return pd.DataFrame(rows).drop_duplicates(['receptor', 'ligand', 'seed'])


aff = {s: load_aff(s) for s in (N, A, B)}
for s, d in aff.items():
    P(f'{s}: {len(d)} afinidades; sementes {sorted(d.seed.unique()) if len(d) else []}')
mean = {s: d.groupby(['receptor', 'ligand']).dG.mean() for s, d in aff.items() if len(d)}

# ---- 1. tabela de energia e postos por pH ----
rk = {}
tab = pd.DataFrame(index=pd.MultiIndex.from_product([RECS, TARGETS], names=['receptor', 'ligand']))
for s in (N, A, B):
    if s in mean:
        tab[f'dG_{s}'] = mean[s].reindex(tab.index)
        tab[f'posto_{s}'] = tab.groupby('receptor')[f'dG_{s}'].rank(method='min')
if f'dG_{A}' in tab and f'dG_{B}' in tab:
    tab[f'delta_{B}-{A}'] = tab[f'dG_{B}'] - tab[f'dG_{A}']
tab.reset_index().to_csv(f'{RES}/ranking_energia_ph.csv', index=False)

md = [f'**ΔG predito (kcal/mol; mais negativo = mais forte). neutro = {N} (padrão Boltz, validado); {A} = pH 8,2 carregado; {B} = pH 10 carregado (exploratório)**', '']
for r in RECS:
    t = tab.loc[r].copy()
    cols = list(dict.fromkeys(c for c in (f'dG_{N}', f'dG_{A}', f'dG_{B}', f'delta_{B}-{A}', f'posto_{N}', f'posto_{A}', f'posto_{B}') if c in t))
    t = t[cols].sort_values(f'posto_{B}' if f'posto_{B}' in t else cols[0])
    md += [f'*{r}*', '', '| Ligante | ' + ' | '.join(cols) + ' |', '|' + '---|' * (len(cols) + 1)]
    for l, row in t.iterrows():
        md.append(f'| {l} | ' + ' | '.join('' if pd.isna(row[c]) else (f'{row[c]:.0f}' if c.startswith('posto') else f'{row[c]:+.2f}' if c.startswith('delta') else f'{row[c]:.2f}') for c in cols) + ' |')
    md.append('')
for s in (N, A, B):
    if f'posto_{s}' in tab:
        pm = tab[f'posto_{s}'].groupby('ligand').mean().reindex(TARGETS).sort_values()
        md.append(f'Posto médio nos 3 receptores ({s}): ' + ', '.join(f'{l} {v:.1f}' for l, v in pm.items()))
open(f'{RES}/ranking_energia_ph.md', 'w', encoding='utf-8').write('\n'.join(md) + '\n')

# ---- 2. o efeito de pH supera o ruído? ----
P('\n== Δ(pH) vs ruído entre sementes (ΔG, kcal/mol) ==')
for r in RECS:
    if f'delta_{B}-{A}' in tab:
        d = tab.loc[r][f'delta_{B}-{A}'].dropna()
        w = None
        try:
            w = stats.wilcoxon(d).pvalue if len(d) >= 6 else None
        except ValueError:
            pass
        P(f'{r}: Δ(pH {B}-{A}) n={len(d)} média={d.mean():+.2f} DP={d.std():.2f} |Δ| mediano={d.abs().median():.2f}' + (f' Wilcoxon p={w:.3f}' if w else ''))
        rho = stats.spearmanr(tab.loc[r][f'dG_{A}'], tab.loc[r][f'dG_{B}'], nan_policy='omit')
        P(f'    preservação da ordem entre {A} e {B}: rho={rho.statistic:.2f}')
for s in (A, B):
    d = aff.get(s)
    if d is not None and {101, 102} <= set(d.seed):
        w = d.pivot_table(index=['receptor', 'ligand'], columns='seed', values='dG').dropna()
        w = w[w.index.get_level_values(1).isin(TARGETS)]
        dd = (w[102] - w[101]).abs()
        P(f'ruído {s} (semente 102 vs 101, alvos): n={len(dd)} |Δ| mediano={dd.median():.2f} DP(Δ)={(w[102] - w[101]).std():.2f}; '
          f'rho entre sementes={stats.spearmanr(w[101], w[102]).statistic:.2f}')
    else:
        P(f'ruído {s}: AUSENTE (precisa das sementes 101 e 102): não concluir efeito de pH sem isso')
for v, lab in (('C10Kn', 'Lys neutra'), ('C10Ym', 'Tyr fenolato')):
    d = load_aff(v)
    if len(d) and B in mean:
        m = d.groupby(['receptor', 'ligand']).dG.mean(); c = sorted(set(m.index) & set(mean[B].index))
        x = (m.reindex(c) - mean[B].reindex(c))
        P(f'variante {v} ({lab}) - {B}: n={len(x)} média={x.mean():+.2f} DP={x.std():.2f} max|Δ|={x.abs().max():.2f}')
P('\n== controles de leitura ==')
for s in (N, A, B):
    if s in mean:
        for r in RECS:
            g2, l1 = mean[s].get((r, 'G2'), np.nan), mean[s].get((r, 'L01'), np.nan)
            a = np.nanmean([mean[s].get((r, l), np.nan) for l in TARGETS]); n = np.nanmean([mean[s].get((r, l), np.nan) for l in NEG])
            P(f'{s:5s} {r:3s} G2={g2:+.2f} L01={l1:+.2f} ordem G2>L01 {"correta" if g2 < l1 else "INVERTIDA"} | alvos={a:+.2f} negativos={n:+.2f} separação={n - a:+.2f}')

# ---- 3. mapa de resíduos <5 Å: o que muda entre os dois pH ----
ca, cb = f'{RES}/contatos_5A_{A}.csv', f'{RES}/contatos_5A_{B}.csv'
cm = ['| Receptor | Ligante | Jaccard consenso | ganhos em pH 10 | perdidos em pH 10 |', '|---|---|---|---|---|']
out = []
if os.path.exists(ca) and os.path.exists(cb):
    ta, tb = pd.read_csv(ca), pd.read_csv(cb)
    for (r, l), ga in ta.groupby(['receptor', 'ligand']):
        gb = tb[(tb.receptor == r) & (tb.ligand == l)]
        if gb.empty:
            continue
        sa, sb = set(ga[ga.freq >= .5].resnum), set(gb[gb.freq >= .5].resnum)
        j = len(sa & sb) / max(len(sa | sb), 1)
        out.append(dict(receptor=r, ligante=l, jaccard=j, ganhos=';'.join(map(str, sorted(sb - sa))), perdidos=';'.join(map(str, sorted(sa - sb)))))
        cm.append(f'| {r} | {l} | {j:.2f} | {";".join(map(str, sorted(sb - sa)))} | {";".join(map(str, sorted(sa - sb)))} |')
    D = pd.DataFrame(out); D.to_csv(f'{RES}/contatos_5A_ph.csv', index=False)
    P('\n== Mapa <5 Å: consenso (freq>=0,5) em ' + A + ' vs ' + B + ' ==')
    for r, g in D.groupby('receptor'):
        chg = pd.Series([x for s in pd.concat([g.ganhos, g.perdidos]) for x in s.split(';') if x]).value_counts().head(6)
        P(f'{r}: Jaccard médio={g.jaccard.mean():.2f} (mín {g.jaccard.min():.2f}); resíduos que mais mudam (nº de ligantes): ' + (', '.join(f'{k}({v})' for k, v in chg.items()) or 'nenhum'))
    P('Sem baseline de ruído por semente para o mapa: mudanças de 1-2 resíduos com freq perto de 0,5 não são interpretáveis.')
    open(f'{RES}/contatos_5A_ph.md', 'w', encoding='utf-8').write('\n'.join(cm) + '\n')
else:
    P(f'\nmapa <5 Å: faltam {os.path.basename(ca)} e/ou {os.path.basename(cb)} (rodar contatos_5A.py --set {A} / --set {B})')

open(f'{RES}/ph_ruido.txt', 'w', encoding='utf-8').write(log.getvalue()); print(log.getvalue())
