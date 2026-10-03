"""Figuras do lote GORE x Boltz-2 (roda no notebook, lê results/*.csv copiados do servidor). Só códigos de ligante."""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

RES = os.path.join(os.path.dirname(__file__), '..', 'results')
FIG = os.path.join(RES, 'figuras'); os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'figure.dpi': 130})
ORDER = ['L01', 'L02', 'L03', 'L04', 'L05', 'L06', 'L07', 'L08', 'L09', 'L10', 'L11', 'G1', 'G2', 'N1', 'N2', 'N3', 'N4']
RECS = ['Ag', 'Sf', 'Bov', 'Ag_D189G', 'Sf_D189G']
COL = {'Ag': '#1b9e77', 'Sf': '#d95f02', 'Bov': '#7570b3'}


def heat(ax, M, rows, cols, title, cmap='viridis', fmt='{:.2f}', vmin=None, vmax=None):
    im = ax.imshow(M, aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=35, ha='right')
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            if not np.isnan(v):
                ax.text(j, i, fmt.format(v), ha='center', va='center', fontsize=6.5, color='white' if (vmax is not None and v < (vmin + vmax) / 2 or vmax is None) and cmap == 'viridis' else 'black')
    ax.set_title(title, fontsize=9)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.02)


def pivot(df, col, arm, recs=RECS, ligs=ORDER):
    d = df[df.arm == arm].pivot(index='ligand', columns='receptor', values=col)
    return d.reindex(index=ligs, columns=recs).values.astype(float)


def fig_qc_receptors():
    f = os.path.join(RES, 'qc_receptors.csv')
    if not os.path.exists(f): return
    d = pd.read_csv(f)
    fig, ax = plt.subplots(1, 4, figsize=(11, 2.9))
    for a, (col, lab, line) in zip(ax, [('mean_plddt', 'pLDDT médio', 90), ('tm_vs_1BTY', 'TM-score vs 1BTY', 0.9),
                                        ('s1_ca_rmsd_A', 'RMSD Cα do S1 (Å)', 1.0), ('his_asp_A', 'His57 Nδ1···Asp102 Oδ (Å)', 3.3)]):
        for i, r in enumerate(['Ag', 'Sf', 'Bov', 'Pig']):
            v = d[d.receptor == r][col]
            a.scatter(np.full(len(v), i) + np.random.RandomState(0).uniform(-.12, .12, len(v)), v, s=14, alpha=.7,
                      color=[('#2ca02c' if p else '#d62728') for p in d[d.receptor == r].qc_pass])
        a.axhline(line, ls='--', c='grey', lw=.8); a.set_xticks(range(4)); a.set_xticklabels(['Ag', 'Sf', 'Bov', 'Pig']); a.set_title(lab, fontsize=9)
    fig.suptitle('QC dos receptores (verde = passa em todos os critérios; linha = limiar pré-registrado)', fontsize=9)
    fig.tight_layout(); fig.savefig(f'{FIG}/fig1_qc_receptores.png'); plt.close(fig)


def fig_heat_interaction(summ):
    fig, ax = plt.subplots(1, 3, figsize=(12, 5.2))
    heat(ax[0], pivot(summ, 'iptm', 'P'), ORDER, RECS, 'ipTM (braço P)', vmin=0.2, vmax=1.0)
    heat(ax[1], pivot(summ, 'ipsae', 'P'), ORDER, RECS, 'ipSAE (PAE<10 Å)', vmin=0.0, vmax=0.8)
    heat(ax[2], pivot(summ, 's1_pocket_frac', 'P'), ORDER, RECS, 'Fração de poses com P1 no bolsão S1 (≤5,5 Å do Cα189)', vmin=0, vmax=1)
    fig.tight_layout(); fig.savefig(f'{FIG}/fig2_interacao_heatmaps.png'); plt.close(fig)


def fig_affinity(summ, aff):
    if aff.empty: return
    a = aff.groupby(['receptor', 'ligand']).agg(y=('log10_ic50_uM', 'mean'), pb=('p_binder', 'mean')).reset_index()
    fig, ax = plt.subplots(1, 2, figsize=(11, 5))
    for k, (col, ttl, cm, vm) in enumerate([('y', 'log10(IC50, µM) previsto (menor = mais forte)', 'viridis_r', (-1, 4)), ('pb', 'Probabilidade de ligante', 'viridis', (0, 1))]):
        M = a.pivot(index='ligand', columns='receptor', values=col).reindex(index=ORDER + ['BEN'], columns=RECS).values.astype(float)
        heat(ax[k], M, ORDER + ['BEN'], RECS, ttl, cmap=cm, vmin=vm[0], vmax=vm[1])
    fig.tight_layout(); fig.savefig(f'{FIG}/fig3_afinidade_heatmaps.png'); plt.close(fig)


def fig_ag_vs_sf(summ, aff):
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.6))
    P = summ[summ.arm == 'P']
    for a_, col, lab in [(ax[0], 'iptm', 'ipTM'), (ax[1], 'ipsae', 'ipSAE')]:
        x = P[P.receptor == 'Ag'].set_index('ligand')[col]; y = P[P.receptor == 'Sf'].set_index('ligand')[col]
        idx = [l for l in ORDER if l in x.index and l in y.index]
        a_.scatter(x[idx], y[idx], c=['#d62728' if l.startswith('N') else '#1f77b4' for l in idx], s=20)
        for l in idx: a_.annotate(l, (x[l], y[l]), fontsize=6, xytext=(2, 2), textcoords='offset points')
        lo, hi = min(x[idx].min(), y[idx].min()) - .03, max(x[idx].max(), y[idx].max()) + .03
        a_.plot([lo, hi], [lo, hi], 'k--', lw=.6); a_.set_xlabel(f'Ag {lab}'); a_.set_ylabel(f'Sf {lab}')
        it = [l for l in idx if not l.startswith('N')]
        if len(it) > 3:
            r_, p_ = stats.spearmanr(x[it], y[it]); a_.set_title(f'{lab}: ρ={r_:.2f} (13 alvos; vermelho = negativos, fora do ρ)', fontsize=8)
    if not aff.empty:
        A = aff.groupby(['receptor', 'ligand']).log10_ic50_uM.mean().unstack(0)
        idx = [l for l in ORDER if l in A.index and not pd.isna(A.loc[l, 'Ag']) and not pd.isna(A.loc[l, 'Sf'])]
        ax[2].scatter(A.loc[idx, 'Ag'], A.loc[idx, 'Sf'], c=['#d62728' if l.startswith('N') else '#1f77b4' for l in idx], s=20)
        for l in idx: ax[2].annotate(l, (A.loc[l, 'Ag'], A.loc[l, 'Sf']), fontsize=6, xytext=(2, 2), textcoords='offset points')
        lo, hi = A.loc[idx, ['Ag', 'Sf']].values.min() - .1, A.loc[idx, ['Ag', 'Sf']].values.max() + .1
        ax[2].plot([lo, hi], [lo, hi], 'k--', lw=.6); ax[2].set_xlabel('Ag log10 IC50 (µM)'); ax[2].set_ylabel('Sf log10 IC50 (µM)')
        it = [l for l in idx if not l.startswith('N')]
        if len(it) > 3:
            r_, p_ = stats.spearmanr(A.loc[it, 'Ag'], A.loc[it, 'Sf']); ax[2].set_title(f'Afinidade: ρ={r_:.2f} (13 alvos)', fontsize=8)
    fig.tight_layout(); fig.savefig(f'{FIG}/fig4_Ag_vs_Sf.png'); plt.close(fig)


def fig_controls(summ, crystal):
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.6))
    P = summ[(summ.arm == 'P') & summ.s1_pocket_frac.notna()]
    lig_kr = [l for l in ORDER if l in P.ligand.unique()]
    w = .2
    for k, r in enumerate(['Ag', 'Ag_D189G', 'Sf', 'Sf_D189G']):
        d = P[P.receptor == r].set_index('ligand').reindex(lig_kr).s1_pocket_frac
        ax[0].bar(np.arange(len(lig_kr)) + (k - 1.5) * w, d.values, w, label=r)
    ax[0].set_xticks(range(len(lig_kr))); ax[0].set_xticklabels(lig_kr, rotation=45); ax[0].set_ylabel('fração de poses com P1 no bolsão S1')
    ax[0].legend(fontsize=6, loc='lower left', ncol=2); ax[0].set_title('Ocupação do S1 e controle D189G do receptor', fontsize=9)
    if crystal is not None and len(crystal):
        for i, c in enumerate(['X1', 'X2', 'X3']):
            d = crystal[crystal.control == c]
            ax[1].scatter(np.full(len(d), i), d.lrmsd_A, c='#1f77b4', s=14)
        ax[1].axhline(2.0, ls='--', c='grey', lw=.8); ax[1].set_xticks(range(3)); ax[1].set_xticklabels(['X1', 'X2', 'X3'])
        ax[1].set_ylabel('LRMSD do peptídeo vs cristal (Å)'); ax[1].set_title('Controles cristalográficos (pré-corte)', fontsize=9)
    fig.tight_layout(); fig.savefig(f'{FIG}/fig5_controles.png'); plt.close(fig)


def fig_cluspro(summ, aff):
    f = os.path.join(RES, '..', 'data', 'cluspro.tsv')
    c = pd.read_csv(f, sep='\t').set_index('code')
    rows = []
    for r in ['Ag', 'Sf']:
        P = summ[(summ.arm == 'P') & (summ.receptor == r)].set_index('ligand')
        A = aff[aff.receptor == r].groupby('ligand').log10_ic50_uM.mean() if not aff.empty else pd.Series(dtype=float)
        for m, ser, sign in [('ipTM', P.iptm, 1), ('ipSAE', P.ipsae, 1), ('log10 IC50', A, -1)]:
            idx = [l for l in c.index if l in ser.index and not pd.isna(ser[l])]
            for cm in ['score_center', 'lowest_energy']:
                if len(idx) >= 4:
                    rho, p = stats.spearmanr(sign * ser[idx], -c.loc[idx, cm])  # sinal: energia mais negativa = melhor
                    rows.append(dict(receptor=r, metrica=m, cluspro=cm, n=len(idx), rho=rho, p=p))
    d = pd.DataFrame(rows); d.to_csv(f'{RES}/cluspro_spearman.csv', index=False)
    if d.empty: return
    fig, ax = plt.subplots(figsize=(6, 3.2))
    lab = [f"{r.receptor}/{r.metrica}/{'SC' if r.cluspro == 'score_center' else 'LE'}" for r in d.itertuples()]
    ax.barh(range(len(d)), d.rho, color=['#1b9e77' if r_ == 'Ag' else '#d95f02' for r_ in d.receptor])
    ax.set_yticks(range(len(d))); ax.set_yticklabels(lab, fontsize=6.5); ax.axvline(0, c='k', lw=.6)
    ax.set_xlabel('Spearman ρ (positivo = concorda com o ClusPro), n=9–11'); ax.set_title('Concordância com o ClusPro (descritiva)', fontsize=9)
    fig.tight_layout(); fig.savefig(f'{FIG}/fig6_cluspro.png'); plt.close(fig)


def fig_pose_qc(poses, summ):
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.2))
    g = poses.groupby('arm').pose_qc.mean()
    ax[0].bar(g.index, g.values * 100, color='#4c72b0'); ax[0].set_ylabel('% de poses que passam no QC'); ax[0].set_ylim(0, 100)
    ax[0].set_title('QC por pose', fontsize=9)
    P = summ[summ.arm == 'P']
    ax[1].hist(P.top_cluster_frac.dropna(), bins=10, color='#4c72b0'); ax[1].set_xlabel('fração de poses no maior cluster'); ax[1].set_title('Consistência entre amostras', fontsize=9)
    ax[2].hist(poses[poses.arm == 'P'].pep_plddt.dropna(), bins=20, color='#4c72b0'); ax[2].set_xlabel('pLDDT do peptídeo'); ax[2].set_title('Confiança local do peptídeo', fontsize=9)
    fig.tight_layout(); fig.savefig(f'{FIG}/fig7_qc_poses.png'); plt.close(fig)


if __name__ == '__main__':
    fig_qc_receptors()
    sf = os.path.join(RES, 'complex_summary.csv')
    if os.path.exists(sf):
        summ = pd.read_csv(sf); poses = pd.read_csv(os.path.join(RES, 'poses.csv'))
        try:
            aff = pd.read_csv(os.path.join(RES, 'affinity.csv'))
        except Exception:
            aff = pd.DataFrame()
        crys = pd.read_csv(os.path.join(RES, 'crystal_controls.csv')) if os.path.exists(os.path.join(RES, 'crystal_controls.csv')) else None
        fig_heat_interaction(summ); fig_affinity(summ, aff); fig_ag_vs_sf(summ, aff); fig_controls(summ, crys)
        fig_cluspro(summ, aff); fig_pose_qc(poses, summ)
    print('figuras em', FIG)
