# GORE-peptides: peptídeos inibidores de tripsina avaliados com Boltz-2

Este repositório reúne o código, as tabelas e as figuras de uma triagem computacional de peptídeos curtos que inibem tripsinas de insetos-praga. A ferramenta usada foi o **Boltz-2**, um modelo de IA que prevê a estrutura de um complexo proteína–ligante e estima a afinidade.

**As sequências dos peptídeos não estão aqui.** Parte delas está sob proteção de patente. Todo o material usa códigos (L01–L11, G1, G2, N1–N4, X1–X3). O arquivo que liga código a sequência fica fora do repositório (`data/private/`, no `.gitignore`).

Autoria: Eulálio (UFV). Status: fase 1 concluída em 03/10/2026, sem validação experimental.

---

## 1. O problema, em poucas linhas

Lagartas como *Anticarsia gemmatalis* e *Spodoptera frugiperda* digerem proteínas com **tripsinas** no intestino médio. Peptídeos pequenos que se encaixam no sítio ativo da tripsina e a bloqueiam podem servir de base para controle de pragas.

Pergunta do projeto: **o Boltz-2 consegue dizer quais peptídeos se ligam melhor às tripsinas de *A. gemmatalis* (Ag), *S. frugiperda* (Sf) e à tripsina bovina (Bov, usada como referência de mamífero)?**

## 2. Como a tripsina reconhece o peptídeo (o mínimo para ler os resultados)

- O **bolsão S1** (resíduos 189–192, 213–220, 226, numeração de tripsina) é fundo e tem um **Asp189** no fundo, com carga negativa. Ele atrai cadeias laterais positivas, Lys ou Arg, do peptídeo. Essa posição do peptídeo chama-se **P1**.
- A **tríade catalítica** (His57, Asp102, Ser195) é quem corta a proteína no sítio ativo.
- Um bom inibidor coloca o P1 básico no S1 e mantém a tríade intacta.

Por isso o projeto inclui um **controle de receptor, o D189G**: troca-se o Asp189 por Gly em Ag e Sf. Sem o Asp189, o P1 básico não deveria mais "entrar" no S1. Se o modelo continua dando nota alta, ele não está enxergando a física do sítio.

## 3. O que foi feito

```
 sequências (privadas)          receptores (5)
        │                             │
        ▼                             ▼
   Boltz-2  ──►  braço P (peptídeo como cadeia)  ──► estrutura, ipTM, ipSAE, contatos
        │
        └──────►  braço S (peptídeo como SMILES)  ──► afinidade (IC50 previsto, ΔG)
                                   │
                  2 sementes (101, 102) • QC pré-registrado • controles
```

| Elemento | Conteúdo |
|---|---|
| Receptores | Ag, Sf, Bov, Ag D189G, Sf D189G (e tripsina suína, só para conferir estruturas de cristal) |
| Ligantes | 13 alvos (L01–L11, G1, G2) + 4 controles negativos (N1–N4) + benzamidina |
| Controles positivos de estrutura | X1–X3, tripeptídeos com estrutura cristalográfica (PDB 5XW8, 5XW9, 5XWJ) |
| Repetição | 2 sementes aleatórias por complexo |
| Estados de ionização | neutro (padrão do Boltz-2, validado) e carregado em pH 8,2 e 10 (exploratório) |

**Escopo deliberadamente estreito:** só Boltz-2 (interação e afinidade), ranking por energia e mapa de resíduos a menos de 5 Å. Sem dinâmica molecular, sem outros programas de predição.

## 4. Métricas, explicadas

- **ipTM e ipSAE:** confiança do modelo na *interface* entre receptor e peptídeo (0 a 1). O ipTM sai alto até para ligantes ruins, porque peptídeos curtos distorcem essa nota; o **ipSAE separa melhor** bons e maus ligantes.
- **Probabilidade de ligante:** o quanto o Boltz-2 acha que a molécula realmente se liga. Separou bem alvos de controles negativos.
- **IC50 previsto → ΔG:** o Boltz-2 devolve log₁₀ IC50 em µM. Convertemos em ΔG = −1,364 · (6 − log₁₀ IC50) kcal/mol para ranquear. **É uma conversão do valor previsto, não energia livre calculada.**
- **Mapa a menos de 5 Å:** para cada pose, quais resíduos da tripsina têm algum átomo a menos de 5 Å do peptídeo.
- **Reprodutibilidade:** a mesma entrada com duas sementes dá quase o mesmo número. Isso mede *precisão*, não *acerto*.

## 5. Resultados principais

1. **A estrutura é boa.** 59 de 60 modelos de receptor, 877 de 880 poses do braço P e 564 de 564 do braço S passam no controle de qualidade. Os cristais X1–X3 são reproduzidos com erro de 0,6 a 0,9 Å (mas são anteriores ao treino do Boltz-2, então testam o procedimento, não a capacidade de generalizar).
2. **O modelo põe o P1 básico no S1** nos receptores nativos e o tira no D189G. A afinidade prevista piora no D189G (Δ log₁₀ IC50 ≈ +0,7, nas duas sementes), mas **a nota de confiança estrutural quase não muda**.
3. **Ranking por energia (neutro):** L03, G2 e L01 no topo; L05 e L11 no fim. O meio do ranking muda muito de um receptor para outro.
4. **A afinidade do Boltz-2 não serve para ordenar a potência real destes peptídeos.** Dos dois peptídeos com Ki experimental conhecido, a ordem de potência só sai correta em 1 dos 3 receptores, e os valores previstos erram o Ki em 1 a 4 ordens de grandeza. É a segunda vez que isso acontece neste projeto.
5. **A incerteza que importa** é a discordância entre os dois modelos internos do Boltz-2 (cerca de 1,2 kcal/mol), muito maior que a variação entre sementes. Diferenças de energia menores que isso entre dois ligantes não sustentam conclusão.
6. **Efeito do pH (8,2 × 10), exploratório:** a ordem dos ligantes se mantém e o ΔG médio quase não muda. Isso **não** mostra que o pH é irrelevante para a tripsina: o Boltz-2 não representa pH, só a carga que se escreve no ligante, e declarar a carga já desloca a previsão mais que o ruído entre sementes.
7. **Insetos × bovina:** a afinidade prevista é mais forte nas tripsinas de inseto, mas sem dado experimental de seletividade isso é uma hipótese.

Detalhes, tabelas e limitações: **[`RELATORIO_FINAL_GORE_Boltz2.md`](RELATORIO_FINAL_GORE_Boltz2.md)** (também em `.docx`).

## 6. Figuras

| Figura | O que mostra |
|---|---|
| `fig1` | Qualidade dos modelos de receptor |
| `fig2` | Interação: ipTM, ipSAE e ocupação do S1 |
| `fig3` | Afinidade prevista (log₁₀ IC50 e probabilidade de ligante) |
| `fig4` | Ag × Sf |
| `fig5` | Controles: D189G e cristais |
| `fig6` | Comparação com o ClusPro |
| `fig7` | Qualidade por pose |
| `fig8` | **Ranking por energia** nos três receptores |
| `fig9` | **Mapa de resíduos a menos de 5 Å** |
| `fig10` | ΔG em pH 8,2 × pH 10 (exploratório) |
| `fig11` | Deslocamento neutro × carregado frente ao ruído |

Todas em `results/figuras/`.

## 7. Estrutura do repositório

```
RELATORIO_FINAL_GORE_Boltz2.md/.docx   relatório completo
bin/                                   scripts (preparo, execução, análise, figuras)
results/                               tabelas e CSVs agregados; results/figuras/ = fig1–fig11
refs/                                  estruturas de referência públicas (1BTY, 5XW8, 5XW9, 5XWJ)
qc_thresholds.yaml                     limiares de QC fixados antes de ver os resultados
data/cluspro.tsv                       pontuações do ClusPro por código de ligante
```

Não incluídos de propósito: sequências dos peptídeos, a tabela de poses (traz o resíduo P1 de cada peptídeo) e as sequências dos receptores de inseto.

## 8. Como reproduzir

Requer GPU com CUDA, Boltz 2.2.1, Python 3.11 com Biopython, RDKit, pandas, SciPy e matplotlib, além do TM-align.

1. Crie `data/private/ligands.tsv` com colunas `code`, `seq`, `group`, `note` (uma linha por ligante; **não versione**).
2. Gere as entradas: `bin/make_inputs.py` (e `make_inputs_ph.py` / `make_inputs_charged.py` para os conjuntos de ionização).
3. Rode o Boltz-2 com `bin/run_boltz.sh <pasta_de_entradas> <pasta_de_saída> <semente>`; para ligante carregado, `bin/run_boltz_nostd.sh`. Os scripts `01–05` orquestram o lote completo.
4. Analise, nesta ordem:
   ```bash
   python bin/analyze.py && python bin/analyze_affinity.py && python bin/compare_crystal.py && python bin/make_tables.py
   python bin/ranking_energia.py
   python bin/contatos_5A.py && python bin/contatos_5A.py --set C82 && python bin/contatos_5A.py --set C10
   python bin/analise_ph.py
   python bin/figures.py && python bin/figures_final.py
   ```
   Os scripts esperam as saídas do Boltz-2 em `~/GOREs-boltz` (ajuste `ROOT` em `bin/gcommon.py`).

Comando-base do Boltz-2: `--model boltz2 --diffusion_samples {5 (peptídeo como cadeia) | 3 (SMILES)} --recycling_steps 3 --sampling_steps 200 --use_potentials --no_kernels`, sementes 101 e 102.

## 9. Limitações

- **Sem validação experimental** dos resultados novos.
- Só dois peptídeos têm Ki experimental conhecido para checar a ordem de potência. O Ki e o pH do ensaio de outros controles não foram obtidos (texto de origem sem acesso aberto).
- Os controles cristalográficos são anteriores ao corte de treino do Boltz-2 (01/06/2023); não há controle posterior ao corte.
- O Boltz-2 não titula nem modela pH. Os conjuntos "carregados" usam uma alteração em processo próprio que desativa a neutralização do SMILES; é uma entrada fora da distribuição de treino do módulo de afinidade.
- Poucos controles negativos (4), poucas sementes (2) e poses por complexo limitam a estatística.

## 10. Citações e ferramentas

- Boltz-2 (modelo de estrutura e afinidade); ipSAE e LIS para interfaces (referências completas a conferir); ClusPro (docking).
- Estruturas cristalográficas de tripeptídeos: Saikhedkar et al., 2018, *Insect Biochem. Mol. Biol.* (PMID 29486250).
- Confira as referências completas no relatório antes de citar em manuscrito.
