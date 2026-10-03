"""Boltz-2 com standardize() desativado (o ligante SMILES entra COM as cargas formais do YAML).

Uso identico ao `boltz`: python bin/boltz_nostd.py predict <yamls> ...   (EXPLORATORIO: o modulo de afinidade foi treinado com
ligante neutro/padronizado; entrada carregada esta fora da distribuicao de treino e NAO tem validacao.)
O patch vale so neste processo; o Boltz instalado (e o lote em andamento) nao e alterado.
"""
from rdkit.Chem import AllChem, Descriptors
AllChem.Descriptors = Descriptors  # o Boltz usa AllChem.Descriptors.MolWt; no caminho padrao o standardize() ja importava isso
import boltz.data.parse.schema as S
S.standardize = lambda smiles: smiles
from boltz.main import cli
if __name__ == '__main__':
    cli()
