import zlib
import matplotlib.pyplot as plt
from Bio import SeqIO

# --- FUNÇÃO 1: A MÁQUINA DE COMPRESSÃO ---
def calcular_complexidade_zlib(sequencia_dna):
    # Converte para bytes (linguagem de máquina)
    bytes_seq = str(sequencia_dna).encode('utf-8')
    # Comprime no nível máximo (9)
    comprimido = zlib.compress(bytes_seq, level=9)
    
    tamanho_original = len(bytes_seq)
    tamanho_comprimido = len(comprimido)
    ratio = tamanho_comprimido / tamanho_original
    
    return tamanho_original, tamanho_comprimido, ratio

# --- FUNÇÃO 2: O GERENTE DE ARQUIVOS ---
def analisar_arquivo(caminho_arquivo, especie):
    print(f"--- Processando {especie} ---")
    try:
        registro = SeqIO.read(caminho_arquivo, "fasta")
        seq = registro.seq
        orig, comp, ratio = calcular_complexidade_zlib(seq)
        
        print(f"   Original: {orig} bytes")
        print(f"   Comprimido: {comp} bytes")
        print(f"   Ratio: {ratio:.4f}") # (Ex: 0.29 significa 29% do tamanho original)
        print("-" * 30)
        return ratio
    except FileNotFoundError:
        print(f"[ERRO] Arquivo {caminho_arquivo} não encontrado!")
        return 0

# --- FUNÇÃO 3: O DESENHISTA ---
def plotar_comparacao(ratio_humano, ratio_chimp):
    print(">>> Gerando Gráfico... <<<")
    etiquetas = ['Humano', 'Chimpanzé']
    valores = [ratio_humano, ratio_chimp]
    cores = ['#1f77b4', '#d62728'] # Azul vs Vermelho
    
    plt.figure(figsize=(8, 5))
    barras = plt.bar(etiquetas, valores, color=cores, width=0.5)
    
    plt.title('Complexidade Algorítmica (Quanto menor, mais estruturado)')
    plt.ylabel('Taxa de Compressão (Ratio)')
    plt.ylim(0, 0.5) # Foco na parte baixa do gráfico
    
    plt.savefig("comparacao_zlib.png", dpi=300)
    print("[SUCESSO] Gráfico salvo como 'comparacao_zlib.png'")

# --- O GATILHO (EXECUÇÃO) ---
if __name__ == "__main__":
    # Aqui é onde a mágica acontece. O Python lê isso e começa a trabalhar.
    print(">>> INICIANDO ANÁLISE AVANÇADA <<<\n")
    
    rh = analisar_arquivo("humano_foxp2.fasta", "Humano")
    rc = analisar_arquivo("chimp_foxp2.fasta", "Chimpanzé")
    
    if rh > 0 and rc > 0:
        plotar_comparacao(rh, rc)