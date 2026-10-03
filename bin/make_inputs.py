"""Gera os YAMLs do Boltz-2 (receptor com MSA em cache; peptideo sem MSA).

Braços (nunca misturar no mesmo `boltz predict`: bug conhecido do passo de afinidade):
  P  = peptideo como cadeia de proteina B (interacao: ipTM/PAE/geometria)
  S  = peptideo como ligante SMILES (zwitterion pH ~8,2) + properties.affinity (afinidade Boltz-2)
Receptores: Ag, Sf, Bov (nativos) e Ag_D189G, Sf_D189G (controle negativo do receptor); Pig p/ controles cristalograficos.
O MSA dos mutantes e o do nativo com a linha-query trocada (MSA-WT: teste de memorizacao do sitio).
"""
import os, re, sys
sys.path.insert(0, os.path.dirname(__file__))
from pep_smiles import peptide_smiles

ROOT = os.path.expanduser('~/GOREs-boltz')
rec = dict(l.split('\t') for l in open(f'{ROOT}/data/private/receptors.tsv').read().splitlines()[1:])
lig = {}
for l in open(f'{ROOT}/data/private/ligands.tsv').read().splitlines()[1:]:
    c = l.split('\t'); lig[c[0]] = (c[1], c[2])
BEN = 'NC(=[NH2+])c1ccccc1'  # benzamidina catiônica (amidínio pKa ~11-12)


def d189_and_more(seq):
    m = re.search(r'GDSGG', seq)
    s195 = m.start() + 2  # S de GDSGG (G193 D194 S195)
    return s195 - 6  # índice 0-based de Asp189 (Ser195-6)


def mut_msa(wt_csv, out_csv, wt_seq, mut_seq):
    lines = open(wt_csv).read().splitlines()
    assert lines[1].split(',', 1)[1] == wt_seq, 'primeira linha do MSA nao e a query'
    lines[1] = lines[1].split(',', 1)[0] + ',' + mut_seq
    open(out_csv, 'w').write('\n'.join(lines) + '\n')


states = {}  # nome -> (seq, msa_csv)
for r in ('Ag', 'Sf', 'Bov', 'Pig'):
    states[r] = (rec[r], f'{ROOT}/msa/receptors_wt/{r}_0.csv')
os.makedirs(f'{ROOT}/msa/mut', exist_ok=True)
for r in ('Ag', 'Sf'):
    seq = rec[r]; i = d189_and_more(seq)
    assert seq[i] == 'D', (r, i, seq[i])
    mseq = seq[:i] + 'G' + seq[i + 1:]
    out = f'{ROOT}/msa/mut/{r}_D189G.csv'
    mut_msa(f'{ROOT}/msa/receptors_wt/{r}_0.csv', out, seq, mseq)
    states[f'{r}_D189G'] = (mseq, out)


def yaml_P(rname, seq, msa, pep):
    return (f"version: 1\nsequences:\n- protein:\n    id: A\n    sequence: {seq}\n    msa: {msa}\n"
            f"- protein:\n    id: B\n    sequence: {pep}\n    msa: empty\n")


def yaml_S(rname, seq, msa, smi):
    return (f"version: 1\nsequences:\n- protein:\n    id: A\n    sequence: {seq}\n    msa: {msa}\n"
            f"- ligand:\n    id: B\n    smiles: '{smi}'\nproperties:\n- affinity:\n    binder: B\n")


plan = []  # (receptor_state, ligand_code)
main = [c for c, (s, g) in lig.items() if g in ('serie_L', 'GORE', 'neg')]
for r in ('Ag', 'Sf', 'Bov', 'Ag_D189G', 'Sf_D189G'):
    for c in main:
        plan.append((r, c))
    plan.append((r, 'BEN'))
for c in ('X1', 'X2', 'X3', 'BEN'):
    plan.append(('Pig', c))

for arm in ('P', 'S'):
    d = f'{ROOT}/inputs/{arm}'
    os.makedirs(d, exist_ok=True)
    for r, c in plan:
        seq, msa = states[r]
        if c == 'BEN':
            if arm == 'P':
                continue  # benzamidina não é peptídeo
            y = yaml_S(r, seq, msa, BEN)
        else:
            pep = lig[c][0]
            y = yaml_P(r, seq, msa, pep) if arm == 'P' else yaml_S(r, seq, msa, peptide_smiles(pep))
        open(f'{d}/{r}__{c}.yaml', 'w').write(y)
    print(arm, len(os.listdir(d)), 'yamls')
open(f'{ROOT}/inputs/plan.tsv', 'w').write('receptor\tligand\n' + '\n'.join(f'{r}\t{c}' for r, c in plan) + '\n')
