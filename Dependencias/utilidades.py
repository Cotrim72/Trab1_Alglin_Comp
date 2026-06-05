import math
import numpy as np

def soma_vetor(a: list, b: list):
    'Retorna a + b. Ambos devem ter o mesmo tamanho.'
    res = [0]*len(a)
    for i in range(len(a)):
        res[i] = a[i] + b[i]
    return res

def sub_vetor(a: list, b: list):
    'Retorna a - b. Ambos devem ter o mesmo tamanho.'
    res = [0]*len(a)
    for i in range(len(a)):
        res[i] = a[i] - b[i]
    return res

def modulo_vetor(a: list):
    soma = 0.0
    for i in range(len(a)):
        val = a[i]
        if math.isinf(val) or math.isnan(val) or abs(val) > 1e154:
            return float('inf') 
        soma += val ** 2
    return soma ** 0.5

def prod_vetor_escalar(a: list, k):
    'Retorna o produto entre o vetor a e o escalar k.'
    res = [0]*len(a)
    for i in range(len(a)):
        res[i] = a[i]*k

    return res

def coluna(A: list[list], j: int):
    'Retorna uma cópia da j-ésima coluna de A.'
    return [A[i][j] for i in range(len(A))]

def prod_matriz_vetor(A: list[list], y: list):
    'Retorna o produto entre A e o vetor y.'
    res = [0]*len(A)
    for i in range(len(A)):
        for j in range(len(y)):
            res[i] += A[i][j]*y[j]

    return res

def erro_solucao(A: list[list], x: list, b: list):
    'Retorna uma métrica para o erro de x enquanto solução do sistema A x = b. Exatamente, retorna o módulo de b - A x; quanto mais próximo de 0, melhor.'
    Ax = prod_matriz_vetor(A, x)
    erro = sub_vetor(b, Ax)
    return modulo_vetor(erro)

def erro_solução_b(A: list[list], x: list, b: list):
    erro = erro_solucao(A,x,b)
    return erro/modulo_vetor(b)   

def gera_matriz_diagonal_dominante(A:list[list]):
    'A é uma matriz quadrada. Altera A e depois a retorna, tornando-a diagonal dominante.'
    for i in range(len(A)):
        soma = 0
        for j in range(len(A[0])):
            if i != j:
                soma += abs(A[i][j])
        A[i][i] = abs(A[i][i]) + soma
    return A

def check_matriz_diagonal_dominante(A:list[list]):
    'A é uma matriz quadrada. Retorna uma string indicando se A é diagonal dominante.'
    for i in range(len(A)):
        soma = 0
        for j in range(len(A[0])):
            if i != j:
                soma += abs(A[i][j])
        if abs(A[i][i]) < soma:
            return "Não é uma matriz diagonal dominante"
    return "A matriz é diagonal dominante"

def ler_matriz_market(caminho: str) -> list[list[float]]:
    '''
    Lê uma matriz no formato MatrixMarket (coordenada) e retorna uma matriz densa como list[list[float]].
    Suporta os tipos de simetria 'general' e 'symmetric' (no caso simétrico, espelha a parte triangular armazenada).
    '''
    with open(caminho, 'r') as f:
        linhas = f.readlines()

    # Cabeçalho: %%MatrixMarket matrix coordinate real <symmetry>
    cabecalho = linhas[0].strip().lower().split()
    if len(cabecalho) < 5 or cabecalho[0] != '%%matrixmarket' or cabecalho[2] != 'coordinate':
        raise ValueError(f'Cabeçalho MatrixMarket inválido ou não suportado: {linhas[0].strip()}')
    simetrica = cabecalho[4] in ('symmetric', 'skew-symmetric')
    anti_simetrica = cabecalho[4] == 'skew-symmetric'

    # Pular comentários e linhas em branco até a linha de dimensões
    i = 1
    while i < len(linhas) and (linhas[i].startswith('%') or linhas[i].strip() == ''):
        i += 1

    partes = linhas[i].split()
    n_linhas = int(partes[0])
    n_colunas = int(partes[1])
    nnz = int(partes[2])

    A = [[0.0]*n_colunas for _ in range(n_linhas)]

    # Entradas: 'linha coluna valor' com índices 1-based
    entradas_lidas = 0
    k = i + 1
    while entradas_lidas < nnz and k < len(linhas):
        linha_txt = linhas[k].strip()
        k += 1
        if linha_txt == '' or linha_txt.startswith('%'):
            continue
        partes = linha_txt.split()
        linha = int(partes[0]) - 1
        coluna = int(partes[1]) - 1
        valor = float(partes[2])
        A[linha][coluna] = valor
        if simetrica and linha != coluna:
            A[coluna][linha] = -valor if anti_simetrica else valor
        entradas_lidas += 1

    return A

def erro_vetor_solucao(x_aleatorio_inicial: list, x_calculado: list):
    dif = sub_vetor(x_aleatorio_inicial, x_calculado)
    return modulo_vetor(dif)

def numero_condicionamento(A:list[list]):
    cond = np.linalg.cond(A)
    return cond

def gerar_matriz_exemplo(n:int):
    A = np.random.uniform(-1, 1, (n, n)) 
    A = (A + A.T) / 2    
    np.fill_diagonal(A, 1.0)                    
    return A.tolist()   

def aplicar_tau(A:list[list],tau:float):
    A_copia = [linha[:] for linha in A]
    for i in range(len(A_copia)):
        for j in range(len(A_copia[0])):
            if i != j: 
                if abs(A_copia[i][j]) > tau:
                    A_copia[i][j] = 0.0
    return A_copia