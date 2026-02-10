import math
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO

def shannon_entropy(dna_sequence):
    """Calcula entropia de Shannon apenas para bases válidas (A, C, T, G)"""
    seq_limpa = ''.join([b for b in dna_sequence.upper() if b in 'ACGT'])
    
    if len(seq_limpa) == 0:
        return 0
    
    contagem = {'A': 0, 'C': 0, 'T': 0, 'G': 0}
    for base in seq_limpa:
        contagem[base] += 1
    
    entropia = 0
    for count in contagem.values():
        if count > 0:
            p = count / len(seq_limpa)
            entropia -= p * math.log2(p)
    
    return entropia

def gerar_dados_janela(arquivo_fasta, especie, janela=100, passo=20):
    """Gera entropias com janela deslizante e validações"""
    print(f"   Processando {especie}...")
    
    try:
        registro = SeqIO.read(arquivo_fasta, "fasta")
        seq = str(registro.seq).upper()
        
        # Remove bases ambíguas
        seq_limpa = ''.join([b for b in seq if b in 'ACGT'])
        removidos = len(seq) - len(seq_limpa)
        
        if removidos > 0:
            print(f"      └─ {removidos} bases ambíguas removidas")
        
        if len(seq_limpa) < janela:
            print(f"      └─ ❌ Sequência muito curta: {len(seq_limpa)} bp")
            return np.array([]), np.array([])
        
        entropias = []
        posicoes = []
        
        for i in range(0, len(seq_limpa) - janela, passo):
            sub_seq = seq_limpa[i:i + janela]
            H = shannon_entropy(sub_seq)
            entropias.append(H)
            posicoes.append(i)
        
        print(f"      └─ ✓ {len(entropias)} janelas analisadas")
        
        return np.array(posicoes), np.array(entropias)
        
    except FileNotFoundError:
        print(f"      └─ ❌ Arquivo não encontrado: {arquivo_fasta}")
        return np.array([]), np.array([])
    except Exception as e:
        print(f"      └─ ❌ Erro: {e}")
        return np.array([]), np.array([])

def plotar_comparacao(janela=100, passo=20):
    """Gera gráfico comparativo de entropia entre espécies"""
    print("=" * 70)
    print("VISUALIZAÇÃO: PERFIL DE ENTROPIA".center(70))
    print("=" * 70 + "\n")
    
    # Processar dados
    x_hum, y_hum = gerar_dados_janela("humano_foxp2.fasta", "Humano", janela, passo)
    x_chimp, y_chimp = gerar_dados_janela("chimp_foxp2.fasta", "Chimpanzé", janela, passo)
    
    # Validações
    if len(x_hum) == 0 or len(x_chimp) == 0:
        print("\n❌ Dados insuficientes para plotagem\n")
        return
    
    # Configurar gráfico
    plt.figure(figsize=(15, 8))
    plt.style.use('seaborn-v0_8-darkgrid')  # Estilo mais moderno
    
    # Plot das curvas
    plt.plot(x_hum, y_hum, label='Homo sapiens', 
             color='#3498db', linewidth=2, alpha=0.8)
    plt.plot(x_chimp, y_chimp, label='Pan troglodytes', 
             color='#e74c3c', linewidth=2, alpha=0.8)
    
    # Linha de referência
    plt.axhline(y=2.0, color='black', linestyle='--', 
                linewidth=1.5, label='Máximo Teórico (Caos)', alpha=0.6)
    
    # Zonas de complexidade
    plt.axhspan(1.5, 1.7, alpha=0.05, color='orange', 
                label='Baixa Complexidade')
    plt.axhspan(1.9, 2.0, alpha=0.05, color='green', 
                label='Alta Complexidade')
    
    # Estatísticas no canto
    media_h = np.mean(y_hum)
    media_c = np.mean(y_chimp)
    desvio_h = np.std(y_hum)
    desvio_c = np.std(y_chimp)
    
    stats_texto = (f"ESTATÍSTICAS:\n"
                  f"─────────────\n"
                  f"Humano:\n"
                  f"  μ = {media_h:.4f}\n"
                  f"  σ = {desvio_h:.4f}\n"
                  f"\n"
                  f"Chimpanzé:\n"
                  f"  μ = {media_c:.4f}\n"
                  f"  σ = {desvio_c:.4f}\n"
                  f"\n"
                  f"Δμ = {abs(media_h - media_c):.4f}")
    
    plt.text(0.02, 0.98, stats_texto,
             transform=plt.gca().transAxes,
             verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', 
                      edgecolor='gray', alpha=0.9),
             fontsize=9, family='monospace')
    
    # Destacar região de maior divergência
    tamanho_min = min(len(y_hum), len(y_chimp))
    diff = np.abs(y_hum[:tamanho_min] - y_chimp[:tamanho_min])
    
    if len(diff) > 0:
        idx_max = np.argmax(diff)
        if diff[idx_max] > 0.05:  # Diferença mínima de 0.05 bits
            plt.annotate('Maior Divergência',
                        xy=(x_hum[idx_max], y_hum[idx_max]),
                        xytext=(x_hum[idx_max] + 200, y_hum[idx_max] + 0.15),
                        arrowprops=dict(arrowstyle='->', color='red', lw=1.5),
                        fontsize=10, color='red',
                        bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))
    
    # Decorações
    plt.title(f'Perfil de Entropia de Shannon - Gene FOXP2\n'
              f'(Janela: {janela} bp, Passo: {passo} bp)',
              fontsize=15, fontweight='bold', pad=20)
    plt.xlabel('Posição na Sequência (bp)', fontsize=12, fontweight='bold')
    plt.ylabel('Entropia de Shannon (bits)', fontsize=12, fontweight='bold')
    plt.legend(loc='lower right', fontsize=10, framealpha=0.9)
    plt.ylim(1.4, 2.05)
    plt.grid(True, alpha=0.3)
    
    # Salvar
    nome_arquivo = f"entropia_perfil_janela{janela}_passo{passo}.png"
    plt.savefig(nome_arquivo, dpi=300, bbox_inches='tight')
    
    print(f"\n✅ Gráfico salvo: {nome_arquivo}")
    print(f"\n📊 RESUMO:")
    print(f"   Humano:     μ={media_h:.4f}, σ={desvio_h:.4f}")
    print(f"   Chimpanzé:  μ={media_c:.4f}, σ={desvio_c:.4f}")
    print(f"   Diferença:  Δμ={abs(media_h - media_c):.4f} bits")
    
    # Interpretação
    if abs(media_h - media_c) < 0.01:
        print(f"\n💡 As sequências apresentam complexidade muito similar")
    elif media_h > media_c:
        print(f"\n💡 FOXP2 humano apresenta maior complexidade informacional")
    else:
        print(f"\n💡 FOXP2 de chimpanzé apresenta maior complexidade informacional")
    
    print("=" * 70 + "\n")

if __name__ == "__main__":
    plotar_comparacao(janela=100, passo=20)