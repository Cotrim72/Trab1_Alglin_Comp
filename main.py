import numpy as np

from testes import TestesUmaMatriz
from Dependencias.utilidades import ler_matriz_market
from Dependencias.sumario import gerar_sumario

t = 0.001
o = 50

matrizes = [
    ('bcsstk01', './Matrizes de Teste/bcsstk01.mtx'),
    ('bcsstk03', './Matrizes de Teste/bcsstk03.mtx'),
    ('bcsstk04', './Matrizes de Teste/bcsstk04.mtx'),
    ('494_bus',  './Matrizes de Teste/494_bus.mtx'),
    ('m1',       './Matrizes de Teste/m1.txt'),
    ('1138_bus', './Matrizes de Teste/1138_bus.mtx'),
]

if __name__ == "__main__":
    todos_resultados = []
    for nome, caminho in matrizes:
        A = ler_matriz_market(caminho)
        n = len(A)
        b = np.random.rand(n).tolist()
        x_inicial = [1]*n

        tr = TestesUmaMatriz(A, b, t, o, x_inicial, f'Resultados/{nome}')
        tr.executar()
        todos_resultados.append((nome, n, tr.resultados))

    gerar_sumario(todos_resultados, 'Resultados/sumario.txt')
