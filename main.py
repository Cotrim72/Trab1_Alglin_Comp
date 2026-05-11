import numpy as np

from testes import TestesUmaMatriz
from Dependencias.utilidades import gera_matriz_diagonal_dominante, ler_matriz_market

n = 3
x_inicial = [1]*n
t = 0.001
o = 50

t = TestesUmaMatriz(
    gera_matriz_diagonal_dominante(np.random.rand(n, n).tolist()),
    np.random.rand(n).tolist(),
    t,
    o,
    x_inicial,
    'Resultados3x3'
)

A_m1 = ler_matriz_market("./Matrizes de Teste/m1.txt")
n_m1 = len(A_m1)
x_inicial_m1 = [1]*n_m1

tr = TestesUmaMatriz(A_m1, np.random.rand(n_m1).tolist(), 0.001, o, x_inicial_m1, 'Resultados')

if __name__ == "__main__":
    # t.executar()
    tr.executar()
