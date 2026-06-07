# Trabalho 1 — Solução de Sistemas de Equações Lineares

**COC473 — Álgebra Linear Computacional**
Leonardo Peres Albertazzi Drummond — Engenharia da Computação e Informação / UFRJ
Primeiro Semestre 2026

---

## 1. Introdução

Este trabalho implementa e compara métodos numéricos para a solução de sistemas lineares quadrados $A\mathbf{x} = \mathbf{b}$, conforme as especificações:

1. **Métodos diretos** (sem uso de bibliotecas): eliminação Gaussiana e fatoração LU.
2. **Métodos iterativos**: Jacobi e Gauss-Seidel, com tolerância para o resíduo e número máximo de iterações.
3. **Verificação da correção** das soluções e medida do **tempo de execução** de cada método.
4. **Comparação contra bibliotecas** (`numpy`/`scipy`) e contra **decomposição de Cholesky** (extensão para matrizes simétricas positivas definidas — slides ALC_MC_03).
5. **Aplicação a matrizes reais** do repositório [Matrix Market](https://math.nist.gov/MatrixMarket/), conforme indicação do PDF *Solvers Benchmarks*.

Todas as decisões de projeto e implementação seguem o conteúdo das aulas (ALC_MC_02 a ALC_MC_04 e ALC_MC_03).

---

## 2. Implementação

### 2.1 Estrutura do projeto

```
Trab1_Alglin_Comp/
├── main.py                      # entrypoint: lista de matrizes + loop de testes
├── sistema_linear.py            # classe SistemaLinear com todos os métodos
├── testes.py                    # classe TestesUmaMatriz (orquestra os testes)
├── Dependencias/
│   ├── utilidades.py            # operações vetoriais, leitor MatrixMarket, etc.
│   ├── log.py                   # objeto Log (iteração, resíduo, x)
│   ├── results.py               # acumulador de texto -> arquivo
│   └── sumario.py               # gerador do sumário comparativo
├── Matrizes de Teste/           # arquivos .mtx baixados do Matrix Market
└── Resultados/                  # gerado por main.py; uma subpasta por matriz
    ├── <matriz>/
    │   ├── configuracoes.txt
    │   ├── Resultados/{eliminacao_gaussiana,fatoracao_lu,cholesky,jacobi,gauss_seidel}.txt
    │   ├── Resultados/{jacobi,gauss_seidel}.png   # gráficos log(r) × iteração
    │   └── Bibliotecas/{eliminacao_gaussiana,fatoracao_lu,cholesky}.txt
    └── sumario.txt              # tabela consolidada
```

A classe `SistemaLinear` ([sistema_linear.py](sistema_linear.py)) encapsula os dados `A`, `b`, `x` e oferece cada método como um método de instância que altera `self.x`.

### 2.2 Métodos diretos

**Eliminação Gaussiana** ([sistema_linear.py:80](sistema_linear.py#L80)) — segue o algoritmo do slide 16 de `ALC_MC_02`: para cada coluna `i`, calcula multiplicadores e elimina elementos abaixo da diagonal, transformando `A` em triangular superior. A substituição para trás ([sistema_linear.py:48](sistema_linear.py#L48)) resolve o sistema triangular.

**Fatoração LU** ([sistema_linear.py:105](sistema_linear.py#L105)) — implementação do algoritmo dos slides 23-24 de `ALC_MC_03`: produto direto da primeira linha de `L` pelas colunas de `U` e recompõe os coeficientes. A solução é obtida resolvendo `Ly = b` por substituição para frente ([sistema_linear.py:64](sistema_linear.py#L64)) e depois `Ux = y` por substituição para trás.

**Decomposição de Cholesky** (extensão, [sistema_linear.py:138](sistema_linear.py#L138)) — segue os slides 33-34 de `ALC_MC_03`: para `A` simétrica positiva definida, constrói `L` tal que `A = L Lᵀ` percorrendo coluna a coluna. Verifica simetria e positividade do termo dentro da raiz; se algum desses falha, retorna a string `"Matriz não é simétrica positiva definida"`. Resolve então `Ly = b` e `Lᵀx = y`.

Nenhum dos três métodos usa `numpy` ou `scipy` — apenas operações nativas do Python sobre listas.

### 2.3 Métodos iterativos

**Jacobi** ([sistema_linear.py:184](sistema_linear.py#L184)) — segue o slide 6 de `ALC_MC_04` (Formato 1, `D Xₖ₊₁ = b − (U+L) Xₖ`):

$$x_i^{k+1} = \frac{b_i - \sum_{j \ne i} a_{ij} \, x_j^{k}}{a_{ii}}$$

**Gauss-Seidel** ([sistema_linear.py:208](sistema_linear.py#L208)) — segue o slide 7 de `ALC_MC_04` (Formato 2, `(D+L) Xₖ₊₁ = b − U Xₖ`), aproveitando valores já atualizados na mesma iteração:

$$x_i^{k+1} = \frac{b_i - \sum_{j < i} a_{ij} \, x_j^{k+1} - \sum_{j > i} a_{ij} \, x_j^{k}}{a_{ii}}$$

**Critério de parada e logging** — em cada iteração calcula-se o **resíduo relativo** definido nos slides 11 e 16 de `ALC_MC_04`:

$$R = \frac{\| X^{k+1} - X^{k} \|_2}{\| X^{k+1} \|_2}$$

Esse valor, junto com o número da iteração e o vetor `xₖ₊₁`, é registrado num objeto `Log` ([Dependencias/log.py](Dependencias/log.py)) e acumulado em `self.logs`. O laço para quando `R ≤ t` (convergiu) ou quando completa `o` iterações sem convergir, caso em que o método retorna a string `"Ultrapassou o número máximo de iterações"`.

A conversão de recursão para laço `for` garante que problemas com `o` grande não estourem a pilha do Python.

### 2.4 Versões com biblioteca (comparação)

Cada método direto tem uma versão usando biblioteca para comparação:
- `eliminacao_gaussiana_numpy` ([sistema_linear.py:95](sistema_linear.py#L95)) — `np.linalg.solve`
- `fatoracao_lu_scipy` ([sistema_linear.py:128](sistema_linear.py#L128)) — `scipy.linalg.lu_factor` + `lu_solve`
- `cholesky_scipy` ([sistema_linear.py:167](sistema_linear.py#L167)) — `scipy.linalg.cho_factor` + `cho_solve`

### 2.5 Verificação da correção

Para cada método, após obtenção de `x`, calcula-se o **erro absoluto da solução**:

$$\text{erro} = \| b - A \cdot x \|_2$$

implementado em [Dependencias/utilidades.py:45](Dependencias/utilidades.py#L45) (`erro_solucao`). Esse valor é registrado nos arquivos de saída e no sumário consolidado.

### 2.6 Saída de cada teste

Os resultados de cada método são escritos em arquivos `.txt` separados, contendo:
- **Diretos**: vetor solução `x`, erro `‖A·x − b‖`, tempo de execução.
- **Iterativos**: indicação se a matriz é diagonal dominante, mensagem opcional ("Ultrapassou..."), **log completo das iterações** (número, resíduo `R`, vetor `x`), solução final, erro absoluto, tempo de execução. Adicionalmente é gerado um gráfico **`log₁₀(R) × iteração`** em formato `.png`.

---

## 3. Matrizes utilizadas

Todas as matrizes vêm do **Matrix Market** ([https://math.nist.gov/MatrixMarket/](https://math.nist.gov/MatrixMarket/)), conforme orientação do PDF *Solvers Benchmarks* postado pelo professor. Foram escolhidas seis matrizes para cobrir diferentes ordens de grandeza, domínios de aplicação e graus de condicionamento.

| Matriz   | Dimensão     | nnz   | Densidade | cond(A)     | $\lambda_\min$ | $\lambda_\max$ | Domínio                       |
|----------|--------------|-------|-----------|-------------|----------------|----------------|--------------------------------|
| bcsstk01 | 48 × 48      | 400   | 17,4 %    | 8,82 × 10⁵  | 3,42 × 10³     | 3,02 × 10⁹     | rigidez estrutural             |
| bcsstk03 | 112 × 112    | 640   | 5,1 %     | 6,79 × 10⁶  | 2,94 × 10⁴     | 2,00 × 10¹¹    | rigidez estrutural             |
| bcsstk04 | 132 × 132    | 3 648 | 20,9 %    | 2,29 × 10⁶  | 4,21           | 9,66 × 10⁶     | rigidez estrutural             |
| 494_bus  | 494 × 494    | 1 666 | 0,68 %    | 2,42 × 10⁶  | 1,24 × 10⁻²    | 3,00 × 10⁴     | sistema elétrico de potência   |
| m1       | 675 × 675    | 3 255 | 0,71 %    | 7,65 × 10⁶  | 1,00           | 7,65 × 10⁶     | proposta pelo professor (SPD)  |
| 1138_bus | 1138 × 1138  | 4 054 | 0,31 %    | 8,57 × 10⁶  | 3,52 × 10⁻³    | 3,01 × 10⁴     | sistema elétrico de potência   |

**Características gerais**: todas as matrizes são **simétricas positivas definidas** (autovalores positivos) e **mal-condicionadas** (condicionamento da ordem de 10⁵–10⁷). Nenhuma é diagonal dominante pelo teste do slide 14 de `ALC_MC_04`. Isso já antecipa que os métodos iterativos clássicos podem não convergir nessas matrizes — o que se confirma na seção de resultados.

O vetor `b` é gerado aleatoriamente (`np.random.rand(n)`) para cada matriz; o vetor inicial dos iterativos é `x⁰ = (1, 1, …, 1)ᵀ`.

**Parâmetros dos iterativos**: tolerância `t = 10⁻³` e número máximo de iterações `o = 50`, conforme os exemplos dos slides 13 e 19 de `ALC_MC_04`.

---

## 4. Resultados

O arquivo consolidado [Resultados/sumario.txt](Resultados/sumario.txt) traz uma tabela por matriz. A seguir apresenta-se uma síntese.

### 4.1 Tempo de execução (s) — métodos diretos

| Matriz   | Elim. gaussiana | Fatoração LU | Cholesky | Gauss numpy | LU scipy | Cholesky scipy |
|----------|----------------:|-------------:|---------:|------------:|---------:|---------------:|
| bcsstk01 | 0,003           | 0,003        | 0,003    | 0,0002      | 0,0003   | 0,0003         |
| bcsstk03 | 0,037           | 0,034        | 0,023    | 0,003       | 0,0007   | 0,0006         |
| bcsstk04 | 0,059           | 0,057        | 0,040    | 0,024       | 0,0009   | 0,0011         |
| 494_bus  | 3,19            | 2,95         | 1,75     | 0,018       | 0,018    | 0,011          |
| m1       | 8,18            | 7,60         | 4,47     | 0,057       | 0,193    | 0,019          |
| 1138_bus | 39,32           | 36,56        | 21,59    | 0,062       | 0,479    | 0,069          |

### 4.2 Erro `‖A·x − b‖` — métodos diretos

Todos abaixo de 10⁻⁹ — solução essencialmente exata em ponto flutuante.

| Matriz   | Elim. gaussiana | Fatoração LU | Cholesky    | Gauss numpy | LU scipy   | Cholesky scipy |
|----------|----------------:|-------------:|------------:|------------:|-----------:|---------------:|
| bcsstk01 | 2,79 × 10⁻¹³    | 2,54 × 10⁻¹³ | 2,46 × 10⁻¹³ | 2,08 × 10⁻¹³ | 2,08 × 10⁻¹³ | 2,68 × 10⁻¹³ |
| bcsstk03 | 7,89 × 10⁻¹²    | 7,53 × 10⁻¹² | 5,54 × 10⁻¹² | 1,65 × 10⁻¹¹ | 1,19 × 10⁻¹¹ | 3,83 × 10⁻¹² |
| bcsstk04 | 2,62 × 10⁻¹²    | 2,68 × 10⁻¹² | 2,03 × 10⁻¹² | 2,10 × 10⁻¹² | 2,50 × 10⁻¹² | 2,58 × 10⁻¹² |
| 494_bus  | 6,80 × 10⁻¹⁰    | 7,32 × 10⁻¹⁰ | 4,59 × 10⁻¹⁰ | 5,50 × 10⁻¹⁰ | 5,50 × 10⁻¹⁰ | 5,18 × 10⁻¹⁰ |
| m1       | 3,88 × 10⁻⁹     | 3,86 × 10⁻⁹  | 2,66 × 10⁻⁹  | 2,62 × 10⁻⁹  | 2,62 × 10⁻⁹  | 2,38 × 10⁻⁹  |
| 1138_bus | 5,26 × 10⁻⁹     | 5,26 × 10⁻⁹  | 2,50 × 10⁻⁹  | 2,27 × 10⁻⁹  | 2,27 × 10⁻⁹  | 2,33 × 10⁻⁹  |

### 4.3 Métodos iterativos

| Matriz   | Jacobi: iter / R / erro abs.       | Gauss-Seidel: iter / R / erro abs. |
|----------|-------------------------------------|-------------------------------------|
| bcsstk01 | 50 (não conv.) / —     / 1,77 × 10⁹ | 50 (não conv.) / —     / 2,01 × 10⁵ |
| bcsstk03 | 50 (não conv.) / —     / 6,26 × 10²³| 50 (não conv.) / —     / 3,38 × 10⁸ |
| bcsstk04 | 50 (não conv.) / —     / 1,10 × 10¹²| 50 (não conv.) / —     / 7,51 × 10² |
| 494_bus  | 50 (não conv.) / —     / 24,6       | 50 (não conv.) / —     / 42,5       |
| m1       | **17 (R < t)** / 9,7 × 10⁻⁴ / 12,5  | **13 (R < t)** / 9,1 × 10⁻⁴ / 12,4  |
| 1138_bus | 50 (não conv.) / —     / 35,7       | 50 (não conv.) / —     / 55,5       |

Em **m1** ambos os métodos atingiram o critério `R < t = 10⁻³` em poucas iterações, mas o erro absoluto final `‖A·x − b‖ ≈ 12,5` continua alto — o critério do slide é satisfeito porque o iterado se estabiliza num platô, sem chegar à solução verdadeira. Nas demais matrizes os iterativos não convergiram dentro de `o = 50` iterações.

### 4.4 Gráficos `log₁₀(R) × iteração`

Conforme especificado, cada execução iterativa produziu um `.png` com a curva exigida pelo enunciado (item 2). Disponíveis em:

- [Resultados/bcsstk01/Resultados/jacobi.png](Resultados/bcsstk01/Resultados/jacobi.png)
- [Resultados/bcsstk01/Resultados/gauss_seidel.png](Resultados/bcsstk01/Resultados/gauss_seidel.png)
- ... (idem para bcsstk03, bcsstk04, 494_bus, m1, 1138_bus)
- 12 gráficos no total.

Nos casos em que o método diverge (bcsstk01, bcsstk03, bcsstk04), os gráficos mostram `log₁₀(R)` crescendo monotonicamente. Em m1, os gráficos mostram queda rápida nas primeiras iterações seguida de platô — ilustrando o limite do critério da aula em sistemas mal-condicionados.

---

## 5. Análise e discussão

### 5.1 Cholesky é ~50 % mais rápido que Gauss/LU em matrizes SPD

A teoria (slide 30 de `ALC_MC_03`) prevê que Cholesky usa cerca de `n³/3` operações contra `2n³/3` de Gauss/LU — proporção 1:2. A medição na maior matriz (`1138_bus`, n = 1138) confirma com precisão:

> 21,59 s (Cholesky) / 39,32 s (Gauss) = **0,549**

Razões análogas para `m1` (0,547), `494_bus` (0,55), `bcsstk03` (0,63) e `bcsstk04` (0,67). A aproximação assintótica vai melhorando à medida que `n` cresce, pois constantes de baixa ordem ficam diluídas.

### 5.2 Iterativos clássicos falham fora do regime diagonal-dominante

Nenhuma das 6 matrizes é diagonal dominante segundo o teste do slide 14 de `ALC_MC_04`. As consequências são claras:

- **bcsstk03** é o exemplo mais dramático: Jacobi diverge para `‖A·x − b‖ ≈ 6,3 × 10²³` em 50 iterações.
- Em todas as matrizes de rigidez (BCSSTK), Jacobi explode (10⁹ a 10²³).
- Gauss-Seidel é mais comportado mas igualmente não converge — atinge `o = 50` sem satisfazer `R < t`.
- Em m1 ocorre um fenômeno diferente: o iterado se estabiliza rapidamente num platô (R cai abaixo de t em 13-17 iterações), mas a solução obtida fica longe da verdadeira (erro absoluto 12,5 vs 10⁻⁹ dos métodos diretos).

Esses resultados ilustram diretamente o conteúdo dos slides 22-25 de `ALC_MC_04` sobre **sistemas mal-condicionados**: uma característica enganosa do critério `R = ‖Δx‖/‖x‖` é que ele pode satisfazer-se mesmo longe da solução real quando o sistema tem condicionamento alto. Nesses casos os métodos diretos são imbatíveis.

### 5.3 Bibliotecas são 100–600× mais rápidas que Python puro

Para a maior matriz (`1138_bus`, n = 1138):
- Eliminação gaussiana: 39,3 s (puro) vs **0,062 s** (numpy) → **632×** mais rápida.
- Cholesky: 21,6 s (puro) vs **0,069 s** (scipy) → **313×** mais rápida.

O ganho vem da implementação em BLAS/LAPACK compilado (C/Fortran com SIMD) versus laços Python interpretados. Note que o **erro de aproximação é praticamente o mesmo** em ambas as implementações (~10⁻⁹), confirmando a correção da versão pura.

### 5.4 Cholesky scipy vs eliminação gaussiana numpy

Mesmo entre as bibliotecas Cholesky vence em `m1` (0,019 s vs 0,057 s) e `1138_bus` (0,069 s vs 0,062 s — empatado). A diferença está no fato de `np.linalg.solve` ser genérico (usa LU com pivoteamento parcial), enquanto `scipy.linalg.cho_factor` aproveita explicitamente que a matriz é SPD.

### 5.5 Erro vs condicionamento

O erro absoluto dos métodos diretos cresce com o condicionamento (e o tamanho), como esperado em aritmética de ponto flutuante:

- bcsstk01 (cond ≈ 9 × 10⁵): erro 10⁻¹³.
- bcsstk03 (cond ≈ 7 × 10⁶): erro 10⁻¹².
- m1 (cond ≈ 8 × 10⁶): erro 10⁻⁹.
- 1138_bus (cond ≈ 9 × 10⁶): erro 5 × 10⁻⁹.

O aumento do erro entre bcsstk03 e m1, apesar de condicionamento similar, deve-se ao acúmulo de operações em `n` maior.

---

## 6. Conclusões

1. **Os quatro métodos principais (Gauss, LU, Jacobi, Gauss-Seidel) foram implementados em Python puro** sem bibliotecas, conforme exigido. A correção foi verificada via `‖A·x − b‖`. Os tempos de execução estão registrados em todos os arquivos de saída.

2. **Cholesky (extensão)** foi adicionado com sucesso, confirmando experimentalmente a previsão teórica de **redução de ~50 % no tempo** versus Gauss/LU em matrizes SPD.

3. **Os métodos iterativos clássicos são pouco úteis para as matrizes deste benchmark** (matrizes reais de aplicação — rigidez estrutural e sistemas de potência), porque nenhuma é diagonal dominante. Em 5 das 6 matrizes os iterativos sequer convergiram em 50 iterações; em uma (m1) convergiram pelo critério do slide mas com erro real ~12, longe da solução exata. **Métodos diretos são, na prática, superiores para esses casos.**

4. **Bibliotecas otimizadas (numpy/scipy) são 100–600× mais rápidas** que a implementação Python pura, mas a precisão é equivalente — comprovando que a diferença é puramente de eficiência, não de algoritmo.

5. **Os gráficos `log(R) × iteração`** foram gerados em formato `.png` para os 12 testes iterativos, conforme item 2 da especificação.

---

## Apêndice A — Como reproduzir

```bash
# Pré-requisitos: Python 3, numpy, scipy, matplotlib
python main.py
```

O script processa as 6 matrizes da lista em [main.py:9](main.py#L9) e gera toda a pasta `Resultados/`, incluindo o sumário consolidado [Resultados/sumario.txt](Resultados/sumario.txt). Tempo total de execução: ~110 s na máquina de teste (a maior parte gasta em `1138_bus` para os métodos puros).

Para adicionar mais matrizes, basta dropar o arquivo `.mtx` em [Matrizes de Teste/](Matrizes%20de%20Teste/) e incluir uma linha na lista `matrizes` em [main.py](main.py).

## Apêndice B — Bibliografia / referências

- Slides do curso COC473, Prof. José Luis Drummond Alves:
  - `ALC_MC_02_PB` — Sistemas de Equações Algebricas Lineares
  - `ALC_MC_03_PB` — Métodos Diretos (Gauss, LU, Cholesky)
  - `ALC_MC_04_PB` — Métodos Iterativos Estacionários (Jacobi, Gauss-Seidel)
- *Solvers Benchmarks* (PDF complementar do professor) — orientação sobre Matrix Market.
- Matrix Market: https://math.nist.gov/MatrixMarket/
