"""SMILES zwitteriônico de peptídeo linear (L-aminoácidos) para o braço de afinidade do Boltz-2.

Estado de protonação em pH ~8,2 (modelo, o Boltz não titula):
  N-terminal NH3+ (pKa ~8: mistura; escolhido o zwitterion padrão), C-terminal COO-,
  Lys NH3+, Arg guanidínio+, Asp/Glu carboxilato. His neutra, Tyr neutra.
"""
from rdkit import Chem


def peptide_smiles(seq: str, nterm_charged: bool = True, lys_charged: bool = True, tyr_anion: bool = False) -> str:
    """Padrão (pH 8,2) inalterado. Microestados para pH 10: nterm_charged=False (NH2 neutro; pKa alfa-amina ~8),
    lys_charged=False (Lys épsilon neutra; pKa ~10,5, ~24% a pH 10), tyr_anion=True (fenolato; pKa ~10,1, ~44% a pH 10)."""
    m = Chem.MolFromSequence(seq, flavor=0)  # flavor 0: peptídeo com terminais livres (NH2 / COOH)
    rw = Chem.RWMol(m)
    patt = {
        'nterm': Chem.MolFromSmarts('[NX3;H2;!$(NC=O)][CX4]C(=O)'),
        'cterm': Chem.MolFromSmarts('C(=O)[OX2H1]'),
        'lys': Chem.MolFromSmarts('[NX3;H2][CH2][CH2][CH2][CH2][CX4]'),
        'arg': Chem.MolFromSmarts('[NX3;H1]([CH2])C(=[NX2;H1])[NX3;H2]'),
        'acid_side': Chem.MolFromSmarts('[CX4][CH2]C(=O)[OX2H1]'),
        'tyr_oh': Chem.MolFromSmarts('c[OX2H1]'),
    }
    def first_atoms(p, idx):
        return [mt[idx] for mt in rw.GetSubstructMatches(p)]
    # ordem: N-terminal = primeiro match do N-terminal (o que tem apenas 1 vizinho C alfa e H2)
    nt = first_atoms(patt['nterm'], 0)
    # seleciona o N-term verdadeiro: N ligado ao C-alfa do primeiro resíduo (menor índice atômico)
    if nt and nterm_charged:
        a = rw.GetAtomWithIdx(min(nt))
        a.SetFormalCharge(1); a.SetNumExplicitHs(3); a.SetNoImplicit(True)
    ct = first_atoms(patt['cterm'], 2)
    if ct:
        a = rw.GetAtomWithIdx(max(ct))  # C-terminal = último OH de ácido a nível de cadeia principal
        a.SetFormalCharge(-1); a.SetNumExplicitHs(0); a.SetNoImplicit(True)
    for i in (first_atoms(patt['lys'], 0) if lys_charged else []):
        a = rw.GetAtomWithIdx(i); a.SetFormalCharge(1); a.SetNumExplicitHs(3); a.SetNoImplicit(True)
    for i in first_atoms(patt['arg'], 3):  # =NH -> =NH2+
        a = rw.GetAtomWithIdx(i); a.SetFormalCharge(1); a.SetNumExplicitHs(2); a.SetNoImplicit(True)
    for i in first_atoms(patt['acid_side'], 4):
        a = rw.GetAtomWithIdx(i); a.SetFormalCharge(-1); a.SetNumExplicitHs(0); a.SetNoImplicit(True)
    if tyr_anion:
        for i in first_atoms(patt['tyr_oh'], 1):
            a = rw.GetAtomWithIdx(i); a.SetFormalCharge(-1); a.SetNumExplicitHs(0); a.SetNoImplicit(True)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    return Chem.MolToSmiles(mol)


if __name__ == '__main__':
    import sys
    from rdkit.Chem import Descriptors
    for s in sys.argv[1:]:
        smi = peptide_smiles(s)
        m = Chem.MolFromSmiles(smi)
        print(s, smi, 'carga', Chem.GetFormalCharge(m), 'heavy', m.GetNumHeavyAtoms())
