"""YAMLs do braço S (afinidade) em pH 10, derivados dos YAMLs de pH 8,2 já usados na semente 101.

O Boltz-2 NÃO modela pH: receptor entra só como sequência+MSA. Logo, em pH 10 só o estado de ionização do LIGANTE
(SMILES) muda; efeitos do lado do receptor (ex.: Ile16 N-terminal, pKa ~9) só aparecem via PROPKA pós-hoc.
O braço P (peptídeo como cadeia de proteína) é invariante ao pH por construção e não é refeito.

Conjuntos (cada um em inputs/<nome>/, mesmos receptores/MSA dos YAMLs de 8,2; só a linha `smiles:` muda):
  S_pH10    : N-terminal NH2 neutro; Lys+ (76%); Arg+; Asp/Glu/C-term -; Tyr neutra (56%)   -> todos os peptídeos
  S_pH10Kn  : como S_pH10, com Lys épsilon neutra (24%)                                     -> só ligantes com K
  S_pH10Ym  : como S_pH10, com Tyr fenolato (44%)                                           -> só ligantes com Y
Benzamidina fica de fora (amidínio pKa ~11-12, SMILES idêntico em 8,2 e 10).
Saída: só códigos de ligante; sequências nunca são impressas.
"""
import os, re, sys, glob
sys.path.insert(0, os.path.dirname(__file__))
from pep_smiles import peptide_smiles
from rdkit import Chem

ROOT = os.path.expanduser('~/GOREs-boltz')
lig = {}
for l in open(f'{ROOT}/data/private/ligands.tsv').read().splitlines()[1:]:
    c = l.split('\t'); lig[c[0]] = c[1]

SETS = {
    'S_pH10':   (dict(nterm_charged=False), lambda s: True),
    'S_pH10Kn': (dict(nterm_charged=False, lys_charged=False), lambda s: 'K' in s),
    'S_pH10Ym': (dict(nterm_charged=False, tyr_anion=True), lambda s: 'Y' in s),
}
SMI = re.compile(r"^(\s+smiles: )'(.*)'$", re.M)


def charge(smi):
    return Chem.GetFormalCharge(Chem.MolFromSmiles(smi))


n_chk = 0
for name in SETS:
    os.makedirs(f'{ROOT}/inputs/{name}', exist_ok=True)
counts = {k: 0 for k in SETS}
same = {}
for f in sorted(glob.glob(f'{ROOT}/inputs/S/*.yaml')):
    base = os.path.basename(f)[:-5]
    rec, code = base.split('__')
    if code == 'BEN':
        continue
    txt = open(f).read()
    m = SMI.search(txt); assert m, f
    seq = lig[code]
    # regressão: o padrão regenerado tem de reproduzir exatamente o SMILES usado na semente 101
    assert peptide_smiles(seq) == m.group(2), f'SMILES padrao mudou em {base}'
    n_chk += 1
    q0 = charge(m.group(2))
    for name, (kw, sel) in SETS.items():
        if not sel(seq):
            continue
        smi = peptide_smiles(seq, **kw)
        if smi == m.group(2):  # ex.: X1-X3 (N-terminal acetilado, sem Lys/Tyr titulavel): identico a 8,2, nao se refaz
            same.setdefault(name, set()).add(code); continue
        open(f'{ROOT}/inputs/{name}/{base}.yaml', 'w').write(txt[:m.start()] + m.group(1) + f"'{smi}'" + txt[m.end():])
        counts[name] += 1
        if rec == 'Ag':  # uma linha por ligante (carga 8,2 -> variante), sem sequência
            print(f'{name:9s} {code:4s} carga {q0:+d} -> {charge(smi):+d}')
print('SMILES identico a 8,2 (pulados):', {k: sorted(v) for k, v in same.items()})
print('regressao 8,2 ok em', n_chk, 'YAMLs; gerados:', counts)
