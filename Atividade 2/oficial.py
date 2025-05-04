# Atividade 2
from Bio import SeqIO
import csv
import os

os.system("cls" if os.name == "nt" else "clear")

print('Iniciando processo.')

# Caminhos dos arquivos
caminho_genoma_completo = "entradas/Ecoli Sakai sequence.fasta"
caminho_plasmídeo_1 = "entradas/Ecoli Sakai Plasmid 1 sequence (1).fasta"
caminho_plasmídeo_2 = "entradas/Ecoli Sakai plasmid 2 sequence (2).fasta"

lista_caminhos = [caminho_genoma_completo, caminho_plasmídeo_1, caminho_plasmídeo_2]

# Criar pasta de saída
os.makedirs("saidas", exist_ok=True)

for caminho in lista_caminhos:
    count = 0

    for registro in SeqIO.parse(caminho, "fasta"):
        seq_str = str(registro.seq).upper()
        numA = 0
        numC = 0
        numG = 0
        numT = 0
        total = 0
        count += 1

        for nucleotideo in seq_str:
            if nucleotideo == "A":
                numA += 1
            elif nucleotideo == "T":
                numT += 1
            elif nucleotideo == "C":
                numC += 1
            elif nucleotideo == "G":
                numG += 1
            total += 1

        nome_saida = os.path.basename(caminho).replace(".fasta", "").replace(" ", "_").replace("(", "").replace(")", "")
        nome_saida = f'{count}_{nome_saida}'
        
        caminho_saida = f"saidas/{nome_saida}_dados.csv"

        with open(caminho_saida, mode="w", newline='', encoding='utf-8') as file:
            writer = csv.writer(file)
            writer.writerow([f"Nome: {registro.name}", f"ID: {registro.id}", f"Descrição: {registro.description}"])
            writer.writerow(["A", "C", "G", "T", "Total"])
            writer.writerow([numA, numC, numG, numT, total])

print('Programa finalizado.')
