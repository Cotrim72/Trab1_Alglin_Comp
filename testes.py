import time
from math import log10
import matplotlib.pyplot as plt
import numpy as np

from Dependencias.results import Results
from Dependencias.utilidades import erro_solucao, check_matriz_diagonal_dominante, erro_vetor_solucao, numero_condicionamento, erro_solução_b
from sistema_linear import SistemaLinear

class TestesUmaMatriz:
    'Objetivo: executar todos os métodos em uma única matriz A e um único vetor B, e escrever os resultados em arquivos'

    def __init__(self, A: list[list], b: list, t: float, o: int, x_inicial: list, caminho_saida: str, x_calculado: list):
        self.A = A
        self.b = b

        self.t = t
        self.o = o
        self.x_inicial = x_inicial

        self.caminho_saida = caminho_saida
        self.x_calculado = x_calculado

        # Cada elemento: {'metodo', 'erro', 'tempo', 'iteracoes' (None se direto), 'mensagem' (None se ok)}
        self.resultados: list[dict] = []

    def sistema(self):
        'Retorna o sistema linear que vamos resolver.'
        return SistemaLinear(self.A, self.b, self.x_inicial)

    def escrever_configuracoes(self):
        # Configurações
        r = Results()
        r.write('Sistema A x = b')

        r.skipline()
        s = self.sistema()
        r.write(f'Matriz A ({s.tamanho()} x {s.tamanho()}):')
        r.write(self.sistema())
        cond = numero_condicionamento(self.A)
        r.skipline()
        r.write(f"κ(A) = {cond:.4e}")

        r.skipline()
        r.write('Vetor b:')
        r.write(self.b)

        r.skipline()
        r.write('Configurações dos métodos iterativos:')
        r.write(f'Tolerância t: {self.t}')
        r.write(f'Número máximo de iterações o: {self.o}')

        r.generate_file(f'{self.caminho_saida}/configuracoes.txt')

    def teste_metodo_direto(self, metodo: str, caminho_saida: str):
        r = Results()
        s = self.sistema()

        inicio = time.time()
        res = getattr(s, metodo)() # Chamando o método
        final = time.time()

        if res != None:
            r.write(res)
            r.skipline()

        x = s.x
        r.write('Solução x:')
        r.write(x)

        r.skipline()
        r.write('X_aleatório_inicial:')
        r.write(self.x_calculado)

        r.skipline()
        erro = erro_solucao(self.A, x, self.b)
        erro_x = erro_vetor_solucao(x,self.x_calculado)
        erro_b = erro_solução_b(self.A, x, self.b)
        r.write(f'Erro da solução (|Ax - b|): {erro}')
        r.write(f'Erro (|Ax - b|)/|b| {erro_b}')
        r.write(f'Erro no vetor x: {erro_x}')
        r.skipline()
        r.write(f'Tempo de execução: {final - inicio}')

        r.generate_file(caminho_saida)

        self.resultados.append({
            'metodo': metodo,
            'erro': erro,
            'tempo': final - inicio,
            'iteracoes': None,
            'mensagem': res,
        })

    def teste_metodo_iterativo(self, metodo: str, caminho_saida: str):
        r = Results()
        s = self.sistema()

        r.write(f'{check_matriz_diagonal_dominante(self.A)}')
        r.skipline()
        
        inicio = time.time()
        res = getattr(s, metodo)(self.t, self.o) # Chamando o método
        final = time.time()

        if res != None:
            r.write(res)
            r.skipline()
            
        x = s.x
        logs = s.logs

        r.write('Logs intermediários:')
        for l in logs: r.write(l)

        r.skipline()
        r.write('Solução x:')
        r.write(x)

        r.skipline()
        r.write('X_aleatório_inicial:')
        r.write(self.x_calculado)

        r.skipline()
        erro = erro_solucao(self.A, x, self.b)
        erro_x = erro_vetor_solucao(x,self.x_calculado)
        erro_b = erro_solução_b(self.A, x, self.b)
        r.write(f'Erro da solução (|Ax - b|): {erro}')
        r.write(f'Erro (|Ax - b|)/|b| {erro_b}')
        r.write(f'Erro no vetor x: {erro_x}')
        r.skipline()
        r.write(f'Tempo de execução: {final - inicio}')

        r.generate_file(caminho_saida)

        self.resultados.append({
            'metodo': metodo,
            'erro': erro,
            'tempo': final - inicio,
            'iteracoes': len(logs),
            'mensagem': res,
        })

        # Gerar gráfico do log10 do resíduo relativo em função da iteração
        iteracoes = [l.iteracao for l in logs if l.residuo > 0]
        log_residuos = [log10(l.residuo) for l in logs if l.residuo > 0]

        plt.figure()
        plt.plot(iteracoes, log_residuos)
        plt.xlabel('Iteração')
        plt.ylabel('log10(R) = log10(||x_novo - x|| / ||x_novo||)')
        plt.title(f'Convergência - {metodo}')
        
        nome_grafico = caminho_saida.replace('.txt', '.png')
        plt.savefig(nome_grafico)
        plt.close()

    def executar(self):
        pasta = self.caminho_saida

        self.escrever_configuracoes()

        self.teste_metodo_direto('eliminacao_gaussiana', f'{pasta}/Resultados/eliminacao_gaussiana.txt')
        self.teste_metodo_direto('fatoracao_lu', f'{pasta}/Resultados/fatoracao_lu.txt')
        self.teste_metodo_direto('cholesky', f'{pasta}/Resultados/cholesky.txt')

        self.teste_metodo_iterativo('jacobi', f'{pasta}/Resultados/jacobi.txt')
        self.teste_metodo_iterativo('gauss_seidel', f'{pasta}/Resultados/gauss_seidel.txt')

        self.teste_metodo_direto('eliminacao_gaussiana_numpy', f'{pasta}/Bibliotecas/eliminacao_gaussiana.txt')
        self.teste_metodo_direto('fatoracao_lu_scipy', f'{pasta}/Bibliotecas/fatoracao_lu.txt')
        self.teste_metodo_direto('cholesky_scipy', f'{pasta}/Bibliotecas/cholesky.txt')
