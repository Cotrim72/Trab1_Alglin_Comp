class Log:
    def __init__(self, iteracao: int, residuo: float, x: list):
        self.iteracao = iteracao
        self.residuo = residuo
        self.x = x

    def __repr__(self):
        return f'Iteração {self.iteracao}: Resíduo relativo = {self.residuo}, x = {self.x}'