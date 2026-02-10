import math
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

def shannon_entropy(dna_sequence):
    """Calcula a Entropia de Shannon (H) para uma string de DNA."""
    tamanho = len(dna_sequence)
    if tamanho == 0: return 0
    contagem = {base: dna_sequence.count(base) for base in set(dna_sequence)}
    entropia = sum(-(count/tamanho) * math.log2(count/tamanho) for count in contagem.values())
    return entropia

def gerar_dados_janela(arquivo_fasta, janela=100, passo=20):
    """Gera lista de entropias usando uma janela deslizante."""
    registro = SeqIO.read(arquivo_fasta, "fasta")
    seq = str(registro.seq).upper()
    
    entropias = []
    posicoes = []
    
    # Janela Deslizante
    for i in range(0, len(seq) - janela, passo):
        sub_seq = seq[i : i + janela]
        H = shannon_entropy(sub_seq)
        entropias.append(H)
        posicoes.append(i)
        
    return posicoes, entropias

def plotar_comparacao():
    print(">>> Gerando Gráfico de Entropia... <<<")
    
    # 1. Processar dados
    x_hum, y_hum = gerar_dados_janela("humano_foxp2.fasta", janela=100)
    x_chimp, y_chimp = gerar_dados_janela("chimp_foxp2.fasta", janela=100)
    
    # 2. Configurar o Gráfico (Estilo Científico)
    plt.figure(figsize=(14, 7))
    plt.style.use('ggplot') # Estilo limpo
    
    # Plot Humano (Azul)
    plt.plot(x_hum, y_hum, label='Homo sapiens (Humano)', color='#1f77b4', linewidth=1.5, alpha=0.9)
    
    # Plot Chimpanzé (Vermelho)
    plt.plot(x_chimp, y_chimp, label='Pan troglodytes (Chimpanzé)', color='#d62728', linewidth=1.5, alpha=0.8)
    
    # 3. Decoração
    plt.title('Densidade de Informação: Complexidade Algorítmica do Gene FOXP2', fontsize=16)
    plt.ylabel('Entropia de Shannon (Bits)', fontsize=12)
    plt.xlabel('Posição no Genoma (Pares de Base)', fontsize=12)
    plt.axhline(y=2.0, color='black', linestyle='--', label='Máxima Aleatoriedade (Caos)')
    
    # Anotação Científica
    plt.text(3000, 1.2, "Região exclusiva humana\n(Extensão Não-Codificante)", 
             fontsize=10, bbox=dict(facecolor='white', alpha=0.8))

    plt.legend(loc='lower right')
    plt.grid(True)
    
    # 4. Salvar
    nome_imagem = "comparacao_entropia_foxp2.png"
    plt.savefig(nome_imagem, dpi=300)
    print(f"[SUCESSO] Gráfico salvo como: {nome_imagem}")

if __name__ == "__main__":
    plotar_comparacao()