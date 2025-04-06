def contar_nucleotideos(sequencia):
    contagem =  {
        "numA": 0,
        "numC": 0,
        "numG": 0,
        "numT": 0,
        "total": 0
    }

    for letra in sequencia:
        if letra == 'A':
            contagem["numA"] += 1
        elif letra == 'C':
            contagem["numC"] += 1
        elif letra == 'G':
            contagem["numG"] += 1
        elif letra == 'T':
            contagem["numT"] += 1
        contagem["total"] += 1

    return contagem

def import_fasta(caminho):
    vetor_dicionarios = list()
    dicionario = dict()
    numSequencia = 0

    with open(caminho, "r") as arquivo:
        for linha in arquivo:
            linha = linha.strip()
            if linha[0] == ">":
                if numSequencia > 0 and dicionario:
                    contagem = contar_nucleotideos(dicionario["sequencia"])
                    dicionario.update(contagem)
                    vetor_dicionarios.append(dicionario) #Insere um dicionári já preenchido

                #Cria um novo dicionário
                dicionario = {
                    "nome" : linha[1:].strip(),
                    "sequencia" : "",
                }
                numSequencia += 1


            else:
                if dicionario:
                    dicionario["sequencia"] += linha
    
    #Inserir quando o dicionário é único ou quando ele é o último
    if dicionario:
        contagem = contar_nucleotideos(dicionario["sequencia"])
        dicionario.update(contagem)
        vetor_dicionarios.append(dicionario)


    return vetor_dicionarios

def exportar_fasta(lista, caminho):
    with open(caminho, 'w') as arquivo:
        arquivo.write("Resultado:\n\n")
        for item in lista:
            total = item["numA"] + item["numC"] + item["numG"] + item["numT"]
            
            arquivo.write(f'{item["sequencia"]}\n\n')
            arquivo.write("Solução\n")
            arquivo.write("-------------------\n")
            arquivo.write(f'Nome do gene: {item["nome"]}\n')
            arquivo.write(f'A: {item["numA"]}\n')
            arquivo.write(f'C: {item["numC"]}\n')
            arquivo.write(f'G: {item["numG"]}\n')
            arquivo.write(f'T: {item["numT"]}\n')
            arquivo.write(f'Total: {total} nucleotídeos\n\n')


caminho_busca = "teste1.fasta"
caminho_saida = "resultado1.txt"
lista = import_fasta(caminho_busca)
exportar_fasta(lista, caminho_saida)

    






