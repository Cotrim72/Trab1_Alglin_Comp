from Dependencias.results import Results

def gerar_sumario(todos_resultados, caminho_saida):
    '''
    todos_resultados: lista de tuplas (nome_matriz, n, resultados), onde
        resultados é a lista de dicts produzida por TestesUmaMatriz.resultados.
    Escreve uma tabela comparativa em formato Markdown.
    '''
    r = Results()
    r.write('# Sumário comparativo dos solvers')
    r.skipline()

    for nome, n, resultados in todos_resultados:
        r.write(f'## {nome} ({n} x {n})')
        r.skipline()
        r.write('| Método                     | Erro (|Ax - b|)/|b| | Tempo (s)  | Notas                |')
        r.write('|----------------------------|---------------------|------------|----------------------|')
        for res in resultados:
            metodo = res['metodo']
            erro = f"{res['erro']:.3e}"
            tempo = f"{res['tempo']:.4f}"
            notas = []
            if res['iteracoes'] is not None:
                notas.append(f"{res['iteracoes']} iter")
            if res['mensagem']:
                notas.append(res['mensagem'])
            nota = '; '.join(notas)
            r.write(f'| {metodo:<26} | {erro:<15} | {tempo:<10} | {nota:<20} |')
        r.skipline()

    r.generate_file(caminho_saida)
