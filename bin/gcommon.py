"""Funções comuns: numeração de tripsina (via alinhamento ao 1BTY), geometria, Kabsch, ipSAE/LIS."""
import os, re, json, glob
import numpy as np
from Bio.PDB import PDBParser
from Bio import Align

ROOT = os.path.expanduser('~/GOREs-boltz')
THREE = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H',
         'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
         'TYR': 'Y', 'VAL': 'V'}
S1_SET = [189, 190, 191, 192, 213, 214, 215, 216, 217, 219, 220, 226]  # bolsão S1 (numeração de tripsina)
CATAL = [57, 102, 195]
OXY = [193, 195]
S2P = [151]  # Tyr151: S'2 (hipótese do projeto análise-alostérica)


def load(path):
    return PDBParser(QUIET=True).get_structure('x', path)[0]


def chain_seq(chain):
    return ''.join(THREE.get(r.get_resname(), 'X') for r in chain if r.id[0] == ' ')


def aligner():
    a = Align.PairwiseAligner()
    a.substitution_matrix = Align.substitution_matrices.load('BLOSUM62')
    a.open_gap_score = -10; a.extend_gap_score = -0.5
    a.mode = 'global'
    return a


def ref_1bty():
    """Sequência e numeração (numeração de tripsina) da cadeia A do 1BTY."""
    m = load(f'{ROOT}/refs/1BTY.pdb')
    ch = m['A']
    res = [r for r in ch if r.id[0] == ' ']
    return m, ''.join(THREE.get(r.get_resname(), 'X') for r in res), [r.id[1] for r in res], res


def number_map(seq):
    """índice 0-based da sequência -> numeração de tripsina (1BTY) por alinhamento global; None se sem par."""
    _, rseq, rnum, _ = ref_1bty()
    al = aligner().align(rseq, seq)[0]
    out = {}
    for (a0, a1), (b0, b1) in zip(*al.aligned):
        for k in range(a1 - a0):
            out[b0 + k] = rnum[a0 + k]
    return out


def kabsch(P, Q):
    """R,t que levam P (Nx3) sobre Q; retorna R,t (P@R.T + t)."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    return R, qc - pc @ R.T


def rmsd(A, B):
    return float(np.sqrt(((A - B) ** 2).sum(1).mean()))


def receptor_residues(model, chain='A'):
    return [r for r in model[chain] if r.id[0] == ' ']


def ca_coords(residues):
    return np.array([r['CA'].coord for r in residues])


def dist(a, b):
    return float(np.linalg.norm(a.coord - b.coord))


def tm_score_from_aligned(P, Q, L):
    """TM-score dado o pareamento (mesma ordem), superposição Kabsch (aproximação; não otimiza)."""
    d0 = 1.24 * (max(L, 19) - 15) ** (1 / 3) - 1.8
    R, t = kabsch(P, Q)
    d = np.linalg.norm(P @ R.T + t - Q, axis=1)
    return float((1 / (1 + (d / d0) ** 2)).sum() / L)


def ipsae_lis(pae, nA, nB, cutoff=10.0):
    """ipSAE (variante d0res, corte de PAE) e LIS entre cadeia A (tokens 0..nA-1) e B (nA..nA+nB-1)."""
    ab = pae[:nA, nA:nA + nB]; ba = pae[nA:nA + nB, :nA]
    def one(M):
        best = 0.0
        for i in range(M.shape[0]):
            v = M[i] < cutoff
            n0 = int(v.sum())
            if n0 == 0:
                continue
            d0 = 1.24 * (max(n0, 27) - 15) ** (1 / 3) - 1.8
            s = float((1 / (1 + (M[i][v] / d0) ** 2)).mean())
            best = max(best, s)
        return best
    ipsae = max(one(ab), one(ba))
    inter = np.concatenate([ab.ravel(), ba.ravel()])
    sel = inter[inter <= 12.0]
    lis = float(((12.0 - sel) / 12.0).mean()) if sel.size else 0.0
    return ipsae, lis
