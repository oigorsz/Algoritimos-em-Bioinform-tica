import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def criaFasta(saida, lista):
    with open(saida, mode='w') as output:
        SeqIO.write(lista, output, "fasta")
    print(f'Arquivo {saida} criado com sucesso.')

def traduz_proteina(frame_rna, tabela):
    proteina = ''
    for i in range(0, len(frame_rna) - 2, 3):
        codon = str(frame_rna[i:i+3])
        proteina += tabela.get(codon, '*')  # Pega '*' se não encontrar
    return proteina

def criaFrames(frame):
    for i in range (1, 7):
        with open(f'saídas/frame{i}.fasta', mode = "w") as output:
            SeqIO.write(frame[i], output, "fasta")
        print(f"Arquivo do Frame {i} gerado com sucesso.")

# Verifica se existe a pasta de saídas
if not os.path.isdir("saídas"):
    os.mkdir("saídas")

# Dicionário código genético RNA → proteína (1 letra)
codon_table = {
    'UUU':'F', 'UUC':'F', 'UUA':'L', 'UUG':'L',
    'UCU':'S', 'UCC':'S', 'UCA':'S', 'UCG':'S',
    'UAU':'Y', 'UAC':'Y', 'UAA':'*', 'UAG':'*',
    'UGU':'C', 'UGC':'C', 'UGA':'*', 'UGG':'W',
    'CUU':'L', 'CUC':'L', 'CUA':'L', 'CUG':'L',
    'CCU':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    'CAU':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
    'CGU':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
    'AUU':'I', 'AUC':'I', 'AUA':'I', 'AUG':'M',
    'ACU':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    'AAU':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
    'AGU':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
    'GUU':'V', 'GUC':'V', 'GUA':'V', 'GUG':'V',
    'GCU':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
    'GAU':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
    'GGU':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'
}

# Criando um dicionário com listas para cada frame
frame = {i: [] for i in range(1, 7)}

entrada = "entradas/Ecoli_Sakai_cds_from_genomic.fna"
rna_fasta = "saídas/Sakai_RNA.fasta"

lista = list()

# Lê as sequências do arquivo de entrada
for sequencia in SeqIO.parse(entrada, "fasta"):
    sequencia_rna = str(sequencia.seq).replace('T', 'U')
    rna_seq = Seq(sequencia_rna)

    # Atualiza informações para fasta de RNA
    sequencia.seq = rna_seq
    sequencia.id = sequencia.id + "_RNA"
    sequencia.description = sequencia.description + " RNA"
    lista.append(sequencia)

    # Frames 1 a 3
    for i in range(3):
        frame_rna = rna_seq[i:]
        proteina = traduz_proteina(frame_rna, codon_table)
        frame[i + 1].append(SeqRecord(Seq(proteina),
                                      id=f"{sequencia.id}_Frame{i+1}",
                                      description=f"Proteína Frame {i+1} de {sequencia.id}"))

    # Frames 4 a 6 (na fita complementar reversa)
    rev_rna_seq = rna_seq.reverse_complement()
    for i in range(3):
        frame_rna = rev_rna_seq[i:]
        proteina = traduz_proteina(frame_rna, codon_table)
        frame[i + 4].append(SeqRecord(Seq(proteina),
                                      id=f"{sequencia.id}_Frame{i+4}",
                                      description=f"Proteína Frame {i+4} de {sequencia.id}"))

#Cria o arquivo fasta de RNA
criaFasta(rna_fasta, lista)

# Cria os arquivos fasta para cada frame
criaFrames(frame)

