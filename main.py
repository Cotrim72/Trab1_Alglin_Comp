import numpy as np

from testes import TestesUmaMatriz
from Dependencias.utilidades import ler_matriz_market, prod_matriz_vetor, gerar_matriz_exemplo
from Dependencias.sumario import gerar_sumario

t = 10**(-8)
o = 2000

#matrizes = [
    #('bcsstk01', './Matrizes de Teste/bcsstk01.mtx'),
    #('bcsstk03', './Matrizes de Teste/bcsstk03.mtx'),
    #('bcsstk04', './Matrizes de Teste/bcsstk04.mtx'),
    #('494_bus',  './Matrizes de Teste/494_bus.mtx'),
    #('m1',       './Matrizes de Teste/m1.txt'),
    #('1138_bus', './Matrizes de Teste/1138_bus.mtx'),
#]

#matrizes = [
    #('bcsstm02', './Matrizes de Teste/bcsstm02.mtx'),
    #('gr_30_30', './Matrizes de Teste/gr_30_30.mtx'),
    #('bcsstm22', './Matrizes de Teste/bcsstm22.mtx'),
    #('bfw398b', './Matrizes de Teste/bfw398b.mtx'),
    #('bfw62b', './Matrizes de Teste/bfw62b.mtx'),
    #('dwb512', './Matrizes de Teste/dwb512.mtx'),
#]

#('Matriz Teste', gerar_matriz_exemplo()),

matrizes = [
    ('Matriz tau = 0.01', gerar_matriz_exemplo(0.01)),
    ('Matriz tau = 0.05', gerar_matriz_exemplo(0.05)),
    ('Matriz tau = 0.1', gerar_matriz_exemplo(0.1)),
    ('Matriz tau = 0.2', gerar_matriz_exemplo(0.2)),
    ('Matriz tau = 0.5', gerar_matriz_exemplo(0.2)),
    ('Matriz tau = 1', gerar_matriz_exemplo(1)),
]


if __name__ == "__main__":
    todos_resultados = []
    for nome, matriz in matrizes:
        #A = ler_matriz_market(matriz)
        A = matriz
        n = len(A)
        vetor_x = np.random.rand(n).tolist()
        b = prod_matriz_vetor(A,vetor_x)
        x_inicial = [1]*n

        tr = TestesUmaMatriz(A, b, t, o, x_inicial, f'Resultados/{nome}', vetor_x)
        tr.executar()
        todos_resultados.append((nome, n, tr.resultados))

    gerar_sumario(todos_resultados, 'Resultados/sumario.txt')
