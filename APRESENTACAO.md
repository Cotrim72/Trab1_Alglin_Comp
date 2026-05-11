---
marp: true
theme: default
paginate: true
header: 'COC473 — Álgebra Linear Computacional — Trabalho 1'
footer: '1º Semestre 2026'
---

# Trabalho 1
## Solução de Sistemas de Equações Lineares

**COC473 — Álgebra Linear Computacional**
Prof. José Luis Drummond Alves

Leonardo Peres Albertazzi Drummond
*(adicionar demais integrantes do grupo)*

ECI / Politécnica — UFRJ — 1º Sem. 2026

---

# Agenda

1. Problema e objetivos
2. Métodos implementados
3. Organização do código
4. Matrizes do *Matrix Market*
5. Resultados
6. Análise
7. Conclusões

---

# 1. O problema

Resolver sistemas lineares quadrados:

$$A\,\mathbf{x} = \mathbf{b}, \qquad A \in \mathbb{R}^{n\times n}$$

**Especificação do trabalho:**

1. Métodos **diretos** — sem `numpy`/`scipy`, verificar correção, medir tempo
2. Métodos **iterativos** — tolerância, máx. iterações, log + gráfico `log(R) × iter`, tempo
3. Usar matrizes de **repositórios públicos**
4. Comparar com bibliotecas (`numpy`/`scipy`)

---

# 2. Métodos implementados

| Categoria   | Método                  | Spec? | Slide referência   |
|-------------|--------------------------|:------:|--------------------|
| Direto      | Eliminação gaussiana     | ✓     | `ALC_MC_02` s.16   |
| Direto      | Fatoração LU             | ✓     | `ALC_MC_03` s.23   |
| Direto      | **Cholesky** *(extra)*   | —     | `ALC_MC_03` s.33   |
| Iterativo   | Jacobi                   | ✓     | `ALC_MC_04` s.11   |
| Iterativo   | Gauss-Seidel             | ✓     | `ALC_MC_04` s.16   |

Cada método tem também uma **versão com biblioteca** (`numpy.linalg.solve`, `scipy.linalg.lu_factor`, `scipy.linalg.cho_factor`) para comparação.

---

# 3. Métodos diretos — Gauss e LU

**Eliminação Gaussiana:** triangulariza `A` via combinações de linhas + retro-substituição.

**Fatoração LU:** `A = L U` (L triangular inferior unitária, U triangular superior). Resolve `Ly = b` (substituição para frente) e `Ux = y` (retro-substituição).

Custo: ambos $\sim \tfrac{2n^3}{3}$ operações.

```python
# Núcleo do LU (sistema_linear.py)
for i in range(n):
    for j in range(i+1, n):
        mult = U[j][i] / U[i][i]
        L[j][i] = mult
        for k in range(i, n):
            U[j][k] -= mult * U[i][k]
```

---

# 3. Métodos diretos — Cholesky (extensão)

Para `A` **simétrica positiva definida**: `A = L Lᵀ`, custo $\sim \tfrac{n^3}{3}$ — **metade** de Gauss/LU.

$$\ell_{ii} = \sqrt{a_{ii} - \sum_{k<i} \ell_{ik}^2}, \qquad
\ell_{ji} = \frac{a_{ji} - \sum_{k<i} \ell_{jk}\ell_{ik}}{\ell_{ii}}$$

Verifica simetria e positividade do radicando; se falham, retorna mensagem de erro.

Resolve `Ly = b` e `Lᵀx = y`, reaproveitando as substituições já implementadas.

---

# 4. Métodos iterativos

**Jacobi** *(slide 11)*:

$$x_i^{k+1} = \frac{b_i - \sum_{j \ne i} a_{ij}\,x_j^{k}}{a_{ii}}$$

**Gauss-Seidel** *(slide 16)* — usa valores já atualizados:

$$x_i^{k+1} = \frac{b_i - \sum_{j<i} a_{ij}\,x_j^{k+1} - \sum_{j>i} a_{ij}\,x_j^{k}}{a_{ii}}$$

**Critério de parada** (slides 11 e 16):

$$R = \frac{\| X^{k+1} - X^{k} \|_2}{\| X^{k+1} \|_2} \le t$$

A cada iteração: grava `(iter, R, x)` em log e gera gráfico `log(R) × iter` ao final.

---

# 5. Organização do código

```
Trab1_Alglin_Comp/
├── main.py                  # loop sobre 6 matrizes
├── sistema_linear.py        # classe SistemaLinear (todos os métodos)
├── testes.py                # orquestra: executa, mede, escreve
└── Dependencias/
    ├── utilidades.py        # vetores, ler_matriz_market, ...
    ├── log.py               # Log(iteracao, residuo, x)
    ├── results.py           # escreve .txt
    └── sumario.py           # gera tabela comparativa
```

`main.py` itera sobre uma lista de matrizes; cada uma gera sua pasta de resultados; ao final, sumário consolidado em **Markdown**.

---

# 6. Matrizes do Matrix Market

Seis matrizes reais, cobrindo dois domínios (engenharia estrutural e sistemas de potência) e três ordens de grandeza:

| Matriz   | n     | cond(A)    | Tipo                 |
|----------|------:|-----------:|----------------------|
| bcsstk01 | 48    | 8,8 × 10⁵  | rigidez estrutural   |
| bcsstk03 | 112   | 6,8 × 10⁶  | rigidez estrutural   |
| bcsstk04 | 132   | 2,3 × 10⁶  | rigidez estrutural   |
| 494_bus  | 494   | 2,4 × 10⁶  | sistema de potência  |
| m1       | 675   | 7,7 × 10⁶  | SPD diag. dominante  |
| 1138_bus | 1138  | 8,6 × 10⁶  | sistema de potência  |

Todas **SPD** (autovalores positivos) e **mal-condicionadas**.

---

# 7. Resultados — Tempo dos diretos (s)

| Matriz   |  Gauss | LU    | Cholesky | numpy | scipy LU | scipy Cholesky |
|----------|-------:|------:|---------:|------:|---------:|---------------:|
| bcsstk01 |  0,003 | 0,003 | 0,003    | 0,0002| 0,0003   | 0,0003         |
| bcsstk04 |  0,059 | 0,057 | 0,040    | 0,024 | 0,0009   | 0,0011         |
| 494_bus  |  3,19  | 2,95  | 1,75     | 0,018 | 0,018    | 0,011          |
| m1       |  8,18  | 7,60  | **4,47** | 0,057 | 0,193    | 0,019          |
| 1138_bus | 39,32  | 36,56 |**21,59** | 0,062 | 0,479    | 0,069          |

**Erro `‖A·x − b‖`** ficou na ordem de **10⁻¹³ a 10⁻⁹** em todas as combinações — solução essencialmente exata.

---

# 7. Resultados — Cholesky confirma a teoria

Teoria: Cholesky usa $n^3/3$ ops; Gauss/LU usam $2n^3/3$ — razão $\tfrac{1}{2}$.

Medido (puros):

| Matriz   | n    | Cholesky / Gauss |
|----------|-----:|------------------|
| bcsstk03 | 112  | 0,63             |
| bcsstk04 | 132  | 0,67             |
| 494_bus  | 494  | 0,55             |
| m1       | 675  | 0,55             |
| **1138_bus** | **1138** | **0,549** |

Para `n` grande, a razão converge ao valor teórico `0,5` (constantes de baixa ordem ficam diluídas).

---

# 7. Resultados — Iterativos

| Matriz   | Jacobi (50 iter)         | Gauss-Seidel (50 iter)    |
|----------|--------------------------|---------------------------|
| bcsstk01 | erro 1,77 × 10⁹           | erro 2,01 × 10⁵            |
| **bcsstk03** | **erro 6,26 × 10²³ (!!)** | erro 3,38 × 10⁸        |
| bcsstk04 | erro 1,10 × 10¹²          | erro 7,51 × 10²            |
| 494_bus  | erro 24,6                 | erro 42,5                  |
| **m1**       | **17 iter — R<t** (erro 12,5) | **13 iter — R<t** (erro 12,4) |
| 1138_bus | erro 35,7                 | erro 55,5                  |

Nenhuma matriz é diagonal dominante ⇒ convergência **não é garantida**. Em 5 das 6, divergem; em m1, "estagnam" no critério mas longe da solução real.

---

# 8. Análise — Sistemas mal-condicionados

Slides 22-25 de `ALC_MC_04` introduzem o problema: pequena variação em `b` causa grande variação em `x` se `det(A) ≈ 0`.

Em **m1** vimos isso na prática:
- O critério `R = ‖Δx‖/‖x‖` cai abaixo de `10⁻³` em 13–17 iterações.
- Mas `‖A·x − b‖ ≈ 12,5` (vs. `10⁻⁹` dos diretos).
- O iterado "estaciona" num platô, sem chegar à solução verdadeira.

**Lição:** o critério da aula é coerente com a definição, mas pode mascarar a falha em sistemas mal-condicionados. Métodos diretos são a escolha segura nesses casos.

---

# 8. Análise — Bibliotecas vs Python puro

Para `1138_bus` (a maior matriz):

| Implementação              | Tempo (s) | Ganho   |
|----------------------------|----------:|--------:|
| Eliminação gaussiana puro  | 39,3      | 1×      |
| `np.linalg.solve`          | 0,062     | **632×** |
| Cholesky puro              | 21,6      | 1×      |
| `scipy.linalg.cho_solve`   | 0,069     | **313×** |

Erro de aproximação **praticamente o mesmo** (~10⁻⁹) — a diferença é puramente de eficiência (BLAS/LAPACK compilado vs Python interpretado).

---

# 9. Gráficos de convergência

Para cada execução iterativa foi gerado um `.png` com `log₁₀(R) × iteração`, como exigido pelo item 2 da spec. **12 gráficos no total**.

- **Casos divergentes** (bcsstk01-04, 494_bus, 1138_bus): curva **crescente** — confirma a divergência observada nas tabelas.
- **m1**: queda rápida nas primeiras iterações + platô — visualiza a "estagnação enganosa" mencionada antes.

Arquivos em `Resultados/<matriz>/Resultados/{jacobi,gauss_seidel}.png`.

---

# 10. Conclusões

1. **Spec atendida**: 4 métodos puros + 3 com biblioteca + Cholesky extra; verificação, tempo, log, gráfico — tudo automático.

2. **Cholesky confirma teoria**: na matriz maior, **0,549 ×** o tempo de Gauss/LU — exatamente a razão `½` prevista.

3. **Iterativos clássicos falham** nessas matrizes reais (nenhuma diagonal dominante). Em 5/6 divergem; em 1, satisfaz o critério mas não a solução verdadeira.

4. **Bibliotecas** são **100-600× mais rápidas** — mesma precisão, eficiência incomparável.

5. Para matrizes SPD mal-condicionadas como as testadas, **métodos diretos (especialmente Cholesky) são a escolha prática**.

---

# Como reproduzir

```bash
python main.py
```

- Lê as 6 matrizes em `Matrizes de Teste/`
- Roda os 8 métodos em cada uma
- Gera `Resultados/<matriz>/{Resultados,Bibliotecas}/*.txt` + `.png`
- Gera `Resultados/sumario.txt` (Markdown comparativo)

**Tempo total:** ~110 s na máquina de teste.

Para adicionar uma matriz: dropar `.mtx` em `Matrizes de Teste/` e incluir na lista em [main.py](main.py).

---

# Obrigado!

**Perguntas?**

Documentos disponíveis no repositório:
- [RELATORIO.md](RELATORIO.md) — relatório completo
- [Resultados/sumario.txt](Resultados/sumario.txt) — tabela comparativa
- [Matrizes de Teste/](Matrizes%20de%20Teste/) — matrizes utilizadas
- Código fonte: [sistema_linear.py](sistema_linear.py), [testes.py](testes.py), [main.py](main.py)
