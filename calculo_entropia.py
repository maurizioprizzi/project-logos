import math
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt

def shannon_entropy(dna_sequence):
    """
    Calcula a Entropia de Shannon para uma string de DNA.
    H(X) = - sum(p(x) * log2(p(x)))
    """
    tamanho = len(dna_sequence)
    if tamanho == 0:
        return 0
    
    # Conta frequência de cada base (A, C, T, G)
    contagem = {base: dna_sequence.count(base) for base in set(dna_sequence)}
    
    entropia = 0
    for base in contagem:
        p_x = contagem[base] / tamanho
        entropia += - p_x * math.log2(p_x)
        
    return entropia

def analisar_entropia_janela(arquivo, nome_especie, tamanho_janela=100):
    """
    Aplica uma janela deslizante para ver a 'densidade de informação' ao longo do gene.
    """
    print(f"--- Calculando Entropia: {nome_especie} ---")
    registro = SeqIO.read(arquivo, "fasta")
    seq = str(registro.seq)
    
    entropias = []
    
    # Janela Deslizante (Sliding Window)
    for i in range(len(seq) - tamanho_janela):
        sub_seq = seq[i : i + tamanho_janela]
        H = shannon_entropy(sub_seq)
        entropias.append(H)
        
    media_entropia = sum(entropias) / len(entropias)
    print(f"Entropia Média: {media_entropia:.4f} bits")
    
    # Retorna os dados para plotagem futura
    return entropias

if __name__ == "__main__":
    print(">>> PROJETO LOGOS: ANÁLISE DE ENTROPIA DE SHANNON <<<\n")
    
    # Calcula
    ent_humano = analisar_entropia_janela("humano_foxp2.fasta", "Humano")
    ent_chimp = analisar_entropia_janela("chimp_foxp2.fasta", "Chimpanzé")
    
    print("\n[ANÁLISE]:")
    print("A Entropia Máxima teórica (caos total) seria 2.0 bits.")
    print("Valores muito baixos indicam alta repetição/ordem.")
    print("Valores intermediários indicam complexidade linguística.")