import zlib
import statistics
import matplotlib.pyplot as plt
from Bio import SeqIO

def calcular_complexidade_zlib(sequencia_dna):
    """
    Calcula métricas de compressibilidade usando ZLIB.
    
    Interpretação:
    - Ratio BAIXO = Alta compressão = Muita redundância (padrões repetitivos)
    - Ratio ALTO = Baixa compressão = Pouca redundância (mais aleatório)
    """
    seq_str = str(sequencia_dna).upper()
    bytes_seq = seq_str.encode('utf-8')
    comprimido = zlib.compress(bytes_seq, level=9)
    
    tamanho_original = len(bytes_seq)
    tamanho_comprimido = len(comprimido)
    ratio = tamanho_comprimido / tamanho_original
    
    return {
        'tamanho_seq': len(seq_str),
        'original_bytes': tamanho_original,
        'comprimido_bytes': tamanho_comprimido,
        'ratio': ratio,
        'economia_percentual': (1 - ratio) * 100,
        'bits_por_base': (tamanho_comprimido * 8) / len(seq_str)
    }

def analisar_janelas(seq, tamanho_janela=200):
    """Analisa complexidade local em janelas deslizantes."""
    seq_str = str(seq).upper()
    seq_limpa = ''.join([n for n in seq_str if n in 'ACGT'])
    
    ratios = []
    for i in range(0, len(seq_limpa) - tamanho_janela, tamanho_janela // 2):
        janela = seq_limpa[i:i+tamanho_janela]
        if len(janela) == tamanho_janela:
            stats = calcular_complexidade_zlib(janela)
            ratios.append(stats['ratio'])
    
    if len(ratios) < 2:
        return None
    
    return {
        'media': statistics.mean(ratios),
        'desvio': statistics.stdev(ratios),
        'min': min(ratios),
        'max': max(ratios),
        'n_janelas': len(ratios)
    }

def analisar_arquivo(caminho_arquivo, especie):
    """Analisa arquivo FASTA e retorna estatísticas de compressão."""
    print(f"--- {especie} ---")
    
    try:
        registro = SeqIO.read(caminho_arquivo, "fasta")
        seq = str(registro.seq).upper()
        
        # Validação
        if len(seq) < 100:
            print(f"⚠️ Sequência muito curta ({len(seq)} bp)\n")
            return None
        
        # Remove ambiguidades
        seq_limpa = ''.join([n for n in seq if n in 'ACGT'])
        removidos = len(seq) - len(seq_limpa)
        
        if removidos > 0:
            print(f"   Bases ambíguas removidas: {removidos}")
        
        # Análise global
        stats_global = calcular_complexidade_zlib(seq_limpa)
        
        # Análise local (janelas)
        stats_janelas = analisar_janelas(seq_limpa)
        
        # Output
        print(f"   Tamanho: {len(seq_limpa):,} bp")
        print(f"   Ratio: {stats_global['ratio']:.4f}")
        print(f"   Economia: {stats_global['economia_percentual']:.1f}%")
        print(f"   Bits/base: {stats_global['bits_por_base']:.2f}")
        
        if stats_janelas:
            print(f"   Variação local: {stats_janelas['desvio']:.4f} "
                  f"(σ em {stats_janelas['n_janelas']} janelas)")
        
        print("-" * 50)
        
        return {**stats_global, 'janelas': stats_janelas, 'especie': especie}
        
    except FileNotFoundError:
        print(f"❌ Arquivo não encontrado: {caminho_arquivo}\n")
        return None
    except Exception as e:
        print(f"❌ Erro ao processar: {e}\n")
        return None

def plotar_comparacao(stats_humano, stats_chimp):
    """Gera visualização comparativa."""
    if not stats_humano or not stats_chimp:
        print("⚠️ Dados insuficientes para plotagem\n")
        return
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    especies = ['Humano', 'Chimpanzé']
    
    # Gráfico 1: Ratio de Compressão
    ratios = [stats_humano['ratio'], stats_chimp['ratio']]
    cores_ratio = ['#3498db', '#e74c3c']
    
    bars1 = ax1.bar(especies, ratios, color=cores_ratio, width=0.6, edgecolor='black')
    ax1.set_ylabel('Ratio de Compressão', fontsize=11)
    ax1.set_title('Compressibilidade\n(Menor = Mais Padrões Repetitivos)', fontsize=12, fontweight='bold')
    ax1.set_ylim(0, max(ratios) * 1.15)
    ax1.grid(axis='y', alpha=0.3, linestyle='--')
    
    for i, (bar, val) in enumerate(zip(bars1, ratios)):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2, height + 0.005,
                f'{val:.4f}', ha='center', va='bottom', fontweight='bold')
    
    # Gráfico 2: Economia de Espaço
    economias = [stats_humano['economia_percentual'], stats_chimp['economia_percentual']]
    cores_econ = ['#2ecc71', '#f39c12']
    
    bars2 = ax2.bar(especies, economias, color=cores_econ, width=0.6, edgecolor='black')
    ax2.set_ylabel('Economia (%)', fontsize=11)
    ax2.set_title('Percentual Comprimido', fontsize=12, fontweight='bold')
    ax2.set_ylim(0, 100)
    ax2.grid(axis='y', alpha=0.3, linestyle='--')
    
    for i, (bar, val) in enumerate(zip(bars2, economias)):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2, height + 2,
                f'{val:.1f}%', ha='center', va='bottom', fontweight='bold')
    
    plt.suptitle('Análise de Complexidade Algorítmica (ZLIB)', 
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig("comparacao_zlib.png", dpi=300, bbox_inches='tight')
    
    print("\n✅ Gráfico salvo: comparacao_zlib.png")
    
    # Interpretação
    diff = abs(stats_humano['ratio'] - stats_chimp['ratio'])
    mais_compressivel = 'Humano' if stats_humano['ratio'] < stats_chimp['ratio'] else 'Chimpanzé'
    
    print(f"\n📊 INTERPRETAÇÃO:")
    print(f"   {mais_compressivel} apresenta MAIOR compressibilidade")
    print(f"   Diferença: {diff:.4f} ({diff/min(ratios)*100:.1f}%)")
    print(f"   Isso sugere mais padrões repetitivos em {mais_compressivel}\n")

if __name__ == "__main__":
    print("=" * 70)
    print("ANÁLISE DE COMPLEXIDADE ALGORÍTMICA (ZLIB)".center(70))
    print("=" * 70 + "\n")
    
    stats_h = analisar_arquivo("humano_foxp2.fasta", "Humano")
    stats_c = analisar_arquivo("chimp_foxp2.fasta", "Chimpanzé")
    
    if stats_h and stats_c:
        plotar_comparacao(stats_h, stats_c)