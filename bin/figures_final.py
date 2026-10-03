"""Figuras finais: ranking por energia, mapa de resíduos <5 Å, dois estados de ionização. Só códigos de ligante."""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RES = os.path.join(os.path.dirname(__file__), '..', 'results')
FIG = os.path.join(RES, 'figuras')
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'figure.dpi': 130})
COL = {'Ag': '#1b9e77', 'Sf': '#d95f02', 'Bov': '#7570b3'}
LIG = ['L01', 'L02', 'L03', 'L04', 'L05', 'L06', 'L07', 'L08', 'L09', 'L10', 'L11', 'G1', 'G2']
NEG = ['N1', 'N2', 'N3', 'N4']
S1 = {189, 190, 191, 192, 213, 214, 215, 216, 217, 219, 220, 226}
CAT = {57, 102, 195}
EXPL = 'exploratório: entrada carregada fora da distribuição de treino do módulo de afinidade; sem validação'

rk = pd.read_csv(f'{RES}/ranking_energia.csv')
ph = pd.read_csv(f'{RES}/ranking_energia_ph.csv')
ct = pd.read_csv(f'{RES}/contatos_5A.csv')


def fig8():
    """ΔG por ligante nos 3 receptores nativos, ordenado pelo posto médio."""
    p = rk.pivot(index='ligand', columns='receptor', values='dG')
    h = rk.pivot(index='ligand', columns='receptor', values='meia_amplitude')
    p = p.reindex([l for l in LIG if l in p.index])
    order = p.assign(m=rk[rk.ligand.isin(LIG)].pivot(index='ligand', columns='receptor', values='posto').mean(axis=1)).sort_values('m').index
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.8), gridspec_kw={'width_ratios': [1.5, 1]})
    x = np.arange(len(order)); w = 0.27
    for k, r in enumerate(['Ag', 'Sf', 'Bov']):
        ax[0].bar(x + (k - 1) * w, p.loc[order, r], w, yerr=h.loc[order, r], color=COL[r], label=r, capsize=1.5)
    ax[0].set_xticks(x); ax[0].set_xticklabels(order, rotation=45)
    ax[0].set_ylabel('ΔG previsto (kcal/mol)'); ax[0].invert_yaxis(); ax[0].legend(frameon=False, ncol=3)
    ax[0].set_title('Ranking por energia (neutro, 2 sementes; barras = meia-amplitude)', fontsize=9)
    M = rk.pivot(index='ligand', columns='receptor', values='posto').reindex(order)[['Ag', 'Sf', 'Bov']]
    im = ax[1].imshow(M.values, cmap='viridis_r', aspect='auto')
    ax[1].set_xticks(range(3)); ax[1].set_xticklabels(M.columns); ax[1].set_yticks(range(len(order))); ax[1].set_yticklabels(order)
    for i in range(M.shape[0]):
        for j in range(3): ax[1].text(j, i, int(M.values[i, j]), ha='center', va='center', color='white' if M.values[i, j] > 7 else 'black', fontsize=8)
    ax[1].set_title('Posto por receptor (1 = mais forte)', fontsize=9)
    fig.tight_layout(); fig.savefig(f'{FIG}/fig8_ranking_energia.png'); plt.close(fig)


def fig9():
    """Mapa de contatos <5 Å (braço P, receptores nativos): frequência de poses por resíduo×ligante."""
    fig, axs = plt.subplots(1, 3, figsize=(13, 5.5))
    for a, r in zip(axs, ['Ag', 'Sf', 'Bov']):
        d = ct[(ct.receptor == r)]
        M = d.pivot_table(index='resnum', columns='ligand', values='freq', aggfunc='max').reindex(columns=LIG + NEG)
        M = M[(M.max(axis=1) >= 0.5)]
        inv = int((M[LIG].min(axis=1) >= 0.9).sum())
        M = M[(M[LIG].max(axis=1) - M[LIG].min(axis=1) >= 0.3)]
        names = d.drop_duplicates('resnum').set_index('resnum').resname
        im = a.imshow(M.values, aspect='auto', cmap='magma_r', vmin=0, vmax=1)
        a.set_xticks(range(M.shape[1])); a.set_xticklabels(M.columns, rotation=90, fontsize=7)
        a.set_yticks(range(M.shape[0]))
        a.set_yticklabels([f'{n}{names[n].title()}' + ('*' if n in S1 else '') + ('†' if n in CAT else '') for n in M.index], fontsize=7)
        a.axvline(12.5, color='k', lw=0.8)
        a.set_title(f'{r}: {M.shape[0]} resíduos variáveis ({inv} constantes omitidos)', fontsize=9)
    fig.colorbar(im, ax=axs, fraction=0.015, label='fração de poses com contato')
    fig.text(0.01, 0.005, '* bolsão S1; † tríade catalítica (57/102/195). Numeração de tripsina (1BTY). Só resíduos com contato ≥0,5 e amplitude ≥0,3 entre os 13 ligantes; constantes (≥0,9 em todos) omitidos. Braço P (peptídeo como cadeia), 2 sementes × 2 modelos.', fontsize=7)
    fig.savefig(f'{FIG}/fig9_mapa_5A.png', bbox_inches='tight'); plt.close(fig)


def fig10():
    """Dois estados de ionização: ΔG C82 × C10 e ordem dos ligantes."""
    d = pd.read_csv(f'{RES}/ranking_energia_ph.csv') if 'dG_C82' in pd.read_csv(f'{RES}/ranking_energia_ph.csv').columns else None
    if d is None: return
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.7), sharex=True, sharey=True)
    for a, r in zip(ax, ['Ag', 'Sf', 'Bov']):
        g = d[d.receptor == r] if 'receptor' in d.columns else None
        a.scatter(g.dG_C82, g.dG_C10, c=COL[r], s=22)
        for _, t in g.iterrows(): a.annotate(t.ligand, (t.dG_C82, t.dG_C10), fontsize=6, xytext=(2, 2), textcoords='offset points')
        lo, hi = -9.5, -5; a.plot([lo, hi], [lo, hi], 'k:', lw=0.8)
        rho = g[['dG_C82', 'dG_C10']].corr(method='spearman').iloc[0, 1]
        a.set_title(f'{r}  ρ={rho:.2f}', fontsize=9); a.set_xlabel('ΔG pH 8,2 carregado (kcal/mol)')
    ax[0].set_ylabel('ΔG pH 10 carregado (kcal/mol)')
    fig.suptitle('Dois estados de ionização do ligante (' + EXPL + ')', fontsize=8, y=1.0)
    fig.tight_layout(); fig.savefig(f'{FIG}/fig10_dois_pH.png', bbox_inches='tight'); plt.close(fig)


def fig11():
    """Neutro (validado) × carregado 8,2: o deslocamento por ligante frente ao ruído entre sementes."""
    d = pd.read_csv(f'{RES}/ranking_energia_ph.csv')
    if 'dG_S' not in d.columns: return
    fig, ax = plt.subplots(figsize=(8, 3.6))
    x = np.arange(len(LIG)); w = 0.27
    for k, r in enumerate(['Ag', 'Sf', 'Bov']):
        g = d[d.receptor == r].set_index('ligand').reindex(LIG)
        ax.bar(x + (k - 1) * w, g.dG_C82 - g.dG_S, w, color=COL[r], label=r)
    ax.axhspan(-0.12, 0.12, color='grey', alpha=0.25, label='ruído entre sementes (±1 DP)')
    ax.set_xticks(x); ax.set_xticklabels(LIG, rotation=45); ax.set_ylabel('ΔG(carregado 8,2) − ΔG(neutro) (kcal/mol)')
    ax.axhline(0, color='k', lw=0.6); ax.legend(frameon=False, ncol=4, fontsize=7)
    ax.set_title('O que muda só por declarar a carga (' + 'exploratório' + ')', fontsize=9)
    fig.tight_layout(); fig.savefig(f'{FIG}/fig11_neutro_vs_carregado.png'); plt.close(fig)


for f in (fig8, fig9, fig10, fig11): f()
print('ok')
