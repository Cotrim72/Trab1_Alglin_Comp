import numpy as np
from scipy.linalg import lu_factor, lu_solve, cho_factor, cho_solve
from scipy.sparse.linalg import gmres

from Dependencias.utilidades import sub_vetor, modulo_vetor, soma_vetor, prod_vetor_escalar
from Dependencias.log import Log

class SistemaLinear:
    'Objeto mutável representando um sistema linear A x = b, em que A é uma matriz quadrada.'

    def __init__(self, A: list[list], b: list, x_inicial: list = None):
        '''
        x_inicial é o ponto de partida para os métodos iterativos. Os atributos não devem ser alterados de fora da classe após a instanciação.
        '''
        self.A = [linha[:] for linha in A]
        self.b = b[:]

        if x_inicial == None:
            self.x = [1]*self.tamanho()
        else:
            self.x = x_inicial[:]

        self.logs: list[Log] = [] # Registros para os métodos iterativos

        # Garantir que a matriz é quadrada e que, em geral, as dimensões são compatíveis
        n = len(self.A)
        assert n >= 1
        for l in self.A: assert len(l) == n
        assert len(self.b) == n
        assert len(self.x) == n

    def __repr__(self):
        'Retorna uma string representando A.'
        s = ''
        s += '['
        for i in range(len(self.A)):
            s += f'{self.A[i]}'
            if i < len(self.A) - 1:
                s += '\n'

        s += ']'
        return s

    def tamanho(self):
        'Retorna a quantidade de linhas de A.'
        return len(self.A)

    def substituicao_para_tras(self):
        'Altera x, realizando a substituição para trás.'
        b = self.b

        res = [0]*len(b)
        res[len(res)-1] = b[len(b)-1]/(self.A[len(self.A)-1][len(self.A[0])-1])
        i = len(b) - 2
        while i >= 0:
            soma = 0
            for j in range(i+1, len(b)):
                soma += (self.A[i][j]*res[j])
            res[i] = (b[i] - soma)/self.A[i][i]
            i -= 1
        
        self.x = res

    def substituiçao_para_frente(self):
        'Altera x, realizando a substituição para frente.'
        b = self.b

        res = [0]*len(b)
        res[0] = b[0]/(self.A[0][0])
        i = 1
        while i < len(b):
            soma = 0
            for j in range(0, i):
                soma += (self.A[i][j]*res[j])
            res[i] = (b[i] - soma)/self.A[i][i]
            i += 1
        
        self.x = res

    def eliminacao_gaussiana(self):
        'Altera A e x, resolvendo o sistema por eliminação gaussiana.'
        b = self.b

        n = len(b)
        for i in range(n):
            pivo = self.A[i][i]
            for j in range(i+1, n):
                multiplicador = self.A[j][i]/pivo
                for k in range(i, n):
                    self.A[j][k] -= self.A[i][k]*multiplicador
                b[j] -= multiplicador*b[i]

        self.substituicao_para_tras()

    def eliminacao_gaussiana_numpy(self):
        'Altera x, realizando a eliminação gaussiana com a biblioteca numpy.'
        b = self.b

        A = np.array(self.A)
        B = np.array(b)
        x = np.linalg.solve(A, B)
        
        self.x = x

    def fatoracao_lu(self):
        'Altera x, resolvendo o sistema com a fatoração LU.'
        b = self.b
        n = len(self.A)

        L = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
        U = [line[:] for line in self.A]

        for i in range(n):
            for j in range(i + 1, n):
                multiplicador = U[j][i] / U[i][i]
                L[j][i] = multiplicador
                for k in range(i, n):
                    U[j][k] -= multiplicador * U[i][k]
        
        fl = SistemaLinear(L, b)
        fl.substituiçao_para_frente()
        y = fl.x

        fu = SistemaLinear(U, y)
        fu.substituicao_para_tras()
        self.x = fu.x

    def fatoracao_lu_scipy(self):
        b = self.b

        A = np.array(self.A)
        B = np.array(b)
        lu, piv = lu_factor(A)
        x = lu_solve((lu, piv), b)

        self.x = x

    def cholesky(self):
        'Altera x, resolvendo o sistema pela decomposição de Cholesky A = L L^T. Aplica-se apenas a matrizes simétricas positivas definidas; retorna mensagem de erro caso contrário.'
        n = self.tamanho()

        for i in range(n):
            for j in range(i + 1, n):
                if self.A[i][j] != self.A[j][i]:
                    return "Matriz não é simétrica positiva definida"

        L = [[0.0]*n for _ in range(n)]
        for i in range(n):
            soma_diag = self.A[i][i] - sum(L[i][k]**2 for k in range(i))
            if soma_diag <= 0:
                return "Matriz não é simétrica positiva definida"
            L[i][i] = soma_diag ** 0.5
            for j in range(i + 1, n):
                soma = self.A[j][i] - sum(L[j][k]*L[i][k] for k in range(i))
                L[j][i] = soma / L[i][i]

        LT = [[L[j][i] for j in range(n)] for i in range(n)]

        fl = SistemaLinear(L, self.b)
        fl.substituiçao_para_frente()
        y = fl.x

        fu = SistemaLinear(LT, y)
        fu.substituicao_para_tras()
        self.x = fu.x

    def cholesky_scipy(self):
        'Altera x, resolvendo o sistema com a decomposição de Cholesky da biblioteca scipy.'
        A = np.array(self.A)
        b = np.array(self.b)
        if not np.array_equal(A, A.T):
            return "Matriz não é simétrica positiva definida"
        try:
            c, low = cho_factor(A)
        except np.linalg.LinAlgError:
            return "Matriz não é simétrica positiva definida"
        self.x = cho_solve((c, low), b)

    def jacobi(self, t: float, o: int):
        'Altera x, resolvendo o sistema pelo método de Jacobi. Para quando o resíduo R = ||x_novo - x|| / ||x_novo|| cai abaixo de t ou após o iterações.'
        b = self.b

        for n in range(1, o + 1):
            x = self.x
            x_novo = [0]*len(x)
            for i in range(len(x)):
                soma = 0
                for j in range(len(x)):
                    if i != j:
                        soma += self.A[i][j]*x[j]
                x_novo[i] = (b[i] - soma)/self.A[i][i]
            residuo = modulo_vetor(sub_vetor(x_novo, x)) / modulo_vetor(x_novo)

            self.logs.append(Log(n, residuo, x_novo))
            self.x = x_novo

            if residuo <= t:
                return

        return "Ultrapassou o número máximo de iterações"
    
    def gauss_seidel(self, t: float, o: int):
        'Altera x, resolvendo o sistema pelo método de Gauss-Seidel. Para quando o resíduo R = ||x_novo - x|| / ||x_novo|| cai abaixo de t ou após o iterações.'
        b = self.b

        for n in range(1, o + 1):
            x = self.x
            x_novo = [0]*len(x)
            for i in range(len(x)):
                soma_new = 0
                soma_old = 0
                for j in range(i):
                    soma_new += self.A[i][j]*x_novo[j]
                for j in range(i+1,len(x)):
                    soma_old += self.A[i][j]*x[j]
                x_novo[i] = (b[i] - soma_new - soma_old)/self.A[i][i]
            residuo = modulo_vetor(sub_vetor(x_novo, x)) / modulo_vetor(x_novo)

            self.logs.append(Log(n, residuo, x_novo))
            self.x = x_novo

            if residuo <= t:
                return

        return "Ultrapassou o número máximo de iterações"