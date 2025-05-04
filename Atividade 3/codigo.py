# Importação das bibliotecas necessárias
from Bio import SeqIO   # Biblioteca para manipulação de arquivos FASTA
import matplotlib.pyplot as plt # Biblioteca para criação de gráficos
import math
import csv  # Biblioteca para leitura e escrita de arquivos CSV
import os   # Biblioteca para manipulação de diretórios e arquivos no sistema operacional

def contar(seq):
    """
    Função que conta a quantidade de nucleotídeos A, T, C e G em uma sequência,
    calcula o conteúdo GC (%) e a temperatura de melting (Tm) baseada no conteúdo GC.
    Retorna um dicionário com esses valores.
    """
    Na_conc = 0.1  # Concentração de Na+, 100mM = 0.1M (valor usado na fórmula de Tm, caso queira trocar a fórmula)

    # Inicializa o dicionário de contagem
    contagem = {
        'A': 0,
        'T': 0,
        'C': 0,
        'G': 0,
        'total': 0,
        'GC': 0,
        'tm': None
    }

    # Conta a ocorrência de cada nucleotídeo na sequência
    for letra in seq:
        if letra == 'A':
            contagem['A'] += 1
        elif letra == 'T':
            contagem['T'] += 1
        elif letra == 'C':
            contagem['C'] += 1
        elif letra == 'G':
            contagem['G'] += 1
        contagem['total'] += 1  # Incírementa o total de nucleotdeos

    # Calcula o conteúdo GC (%) e a Temperatura de Melting (Tm) se houver sequência
    if contagem['total'] > 0:
        contagem['GC'] = ((contagem['G'] + contagem['C']) / contagem['total']) * 100
        contagem['tm'] = 64.9 + 0.41 * contagem['GC'] - (500 / contagem['total'])  # Fórmula simplificada para Tm
        # Fórmula alternativa comentada:
        # contagem['tm'] = 81.5 + 16.6*math.log10(Na_conc) + 0.41*contagem['GC'] - (500/contagem['total'])
    else:
        contagem['GC'] = 0  # Evita divisão por zero

    return contagem  # Retorna o dicionário com os resultados

def geraGrafico(lista):
    """
    Função que gera e salva um gráfico de dispersão do conteúdo GC (%) 
    versus a Temperatura de Melting (Tm) das sequências analisadas.
    """
    # Obtém listas de GC e Tm a partir da lista de sequências
    gc_values = [item.gc for item in lista]
    tm_values = [item.tm for item in lista]

    plt.figure(figsize=(8, 6))
    plt.scatter(tm_values, gc_values, color='red', marker='x')

    # Define rótulos, título e grade do gráfico
    plt.xlabel('Temperatura de Melting (°C)')
    plt.ylabel('Conteúdo GC (%)')
    plt.title('Conteúdo GC vs Temperatura de Melting\nIgor de Souza Monteiro')
    plt.grid(True)

    # Salva o gráfico em arquivo e exibe na tela
    plt.savefig('saidas/grafico_GC_vs_Tm.png')
    plt.show()

    print(f"Arquivo PNG gerado com sucesso em {grafico_png}")

# Define o caminho do arquivo FASTA de entrada
arquivo_fasta = "entradas/Ecoli_Sakai_cds_from_genomic.fna"

# Cria a lista que armazenará as sequências com suas informações adicionais
lista = list()

# Leitura do arquivo FASTA
for sequencia in SeqIO.parse(arquivo_fasta, "fasta"):
    # Calcula as informações para a sequência atual
    contagem = contar(sequencia)

    # Atribui os valores obtidos como atributos do objeto sequencia
    sequencia.numA = contagem['A']
    sequencia.numT = contagem['T']
    sequencia.numC = contagem['C']
    sequencia.numG = contagem['G']
    sequencia.total = contagem['total']
    sequencia.gc = contagem['GC']
    sequencia.tm = contagem['tm']

    # Adiciona o objeto sequencia à lista de resultados
    lista.append(sequencia)

# Verifica e cria o diretório de saída caso não exista
if not os.path.exists('saidas'):
    os.makedirs('saidas')

# Define os caminhos dos arquivos de saída
dados_da_sequencia_csv = 'saidas/dados_das_sequencias.csv'
conteudo_gc_csv = 'saidas/conteudo_GC.csv'
temperatura_x_gc_csv = 'saidas/caminho_temperatura_x_gc.csv'
grafico_png = 'saidas/grafico_GC_vs_Tm.png'

# Gera o arquivo CSV com as quantidades de nucleotídeos e tamanho total das sequências
with open(dados_da_sequencia_csv, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['ID', 'A', 'T', 'C', 'G', 'Total'])  # Cabeçalho do arquivo

    for item in lista:
        writer.writerow([
            item.id,
            item.numA,
            item.numT,
            item.numC,
            item.numG,
            item.total
        ])

    print(f"Arquivo CSV gerado com sucesso em {dados_da_sequencia_csv}")

# Gera o arquivo CSV com o conteúdo GC (%) de cada sequência
with open(conteudo_gc_csv, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['ID', 'GC'])

    for item in lista:
        writer.writerow([
            item.id,
            f"{item.gc / 100:.3f}"  # Salva GC como fração (0 a 1) com 3 casas decimais
        ])

    print(f"Arquivo CSV gerado com sucesso em {conteudo_gc_csv}")

# Gera o arquivo CSV com os valores de Tm e GC para posterior análise gráfica
with open(temperatura_x_gc_csv, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Temperatura de Melting (Cº)', 'Conteúdo GC (%)'])

    for item in lista:
        writer.writerow([
            f'{item.tm:.3f}',
            f'{item.gc:.2f}'
        ])

    print(f"Arquivo CSV gerado com sucesso em {temperatura_x_gc_csv}")

# Chama a função para gerar o gráfico
geraGrafico(lista)
