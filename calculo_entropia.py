import math
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from scipy import stats as sp_stats

def shannon_entropy(dna_sequence):
    """
    Calcula a Entropia de Shannon para DNA (apenas A, C, T, G).
    H(X) = -Σ p(x) * log₂(p(x))
    
    Máximo teórico: 2.0 bits (distribuição uniforme)
    """
    seq_limpa = ''.join([b for b in dna_sequence.upper() if b in 'ACGT'])
    
    if len(seq_limpa) == 0:
        return 0
    
    contagem = {'A': 0, 'C': 0, 'T': 0, 'G': 0}
    for base in seq_limpa:
        contagem[base] += 1
    
    entropia = 0
    tamanho = len(seq_limpa)
    
    for base in contagem:
        if contagem[base] > 0:
            p_x = contagem[base] / tamanho
            entropia -= p_x * math.log2(p_x)
    
    return entropia

def analisar_entropia_janela(arquivo, nome_especie, tamanho_janela=100):
    """Calcula entropia com janela deslizante."""
    print(f"--- {nome_especie} ---")
    
    try:
        registro = SeqIO.read(arquivo, "fasta")
        seq = str(registro.seq).upper()
        
        # Remove bases ambíguas
        seq_limpa = ''.join([b for b in seq if b in 'ACGT'])
        removidos = len(seq) - len(seq_limpa)
        
        if removidos > 0:
            print(f"   Bases ambíguas removidas: {removidos}")
        
        if len(seq_limpa) < tamanho_janela:
            print(f"   ❌ Sequência muito curta ({len(seq_limpa)} bp)\n")
            return None
        
        # Entropia global
        entropia_global = shannon_entropy(seq_limpa)
        
        # Janela deslizante
        entropias = []
        for i in range(len(seq_limpa) - tamanho_janela + 1):
            sub_seq = seq_limpa[i:i + tamanho_janela]
            H = shannon_entropy(sub_seq)
            entropias.append(H)
        
        entropias = np.array(entropias)
        
        # Estatísticas
        print(f"   Tamanho: {len(seq_limpa):,} bp")
        print(f"   Entropia global: {entropia_global:.4f} bits")
        print(f"   Entropia média (janelas): {np.mean(entropias):.4f} bits")
        print(f"   Desvio padrão: {np.std(entropias):.4f}")
        print(f"   Range: [{np.min(entropias):.4f}, {np.max(entropias):.4f}]")
        print(f"   Janelas: {len(entropias):,}")
        print("-" * 50)
        
        return {
            'entropias': entropias,
            'global': entropia_global,
            'media': np.mean(entropias),
            'desvio': np.std(entropias),
            'min': np.min(entropias),
            'max': np.max(entropias),
            'especie': nome_especie,
            'tamanho_seq': len(seq_limpa)
        }
        
    except FileNotFoundError:
        print(f"   ❌ Arquivo não encontrado: {arquivo}\n")
        return None
    except Exception as e:
        print(f"   ❌ Erro: {e}\n")
        return None

def plotar_analise_completa(stats_h, stats_c, tamanho_janela):
    """Visualização completa da análise de entropia"""
    if not stats_h or not stats_c:
        print("⚠️ Dados insuficientes para plotagem")
        return
    
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)
    
    # 1. Perfil ao longo da sequência
    ax1 = fig.add_subplot(gs[0, :])
    ax1.plot(stats_h['entropias'], label='Humano', color='#3498db', alpha=0.7, linewidth=0.8)
    ax1.plot(stats_c['entropias'], label='Chimpanzé', color='#e74c3c', alpha=0.7, linewidth=0.8)
    ax1.axhline(y=2.0, color='black', linestyle='--', linewidth=1, 
                label='Máximo Teórico', alpha=0.5)
    ax1.set_xlabel(f'Posição (janelas de {tamanho_janela} bp)')
    ax1.set_ylabel('Entropia (bits)')
    ax1.set_title('Perfil de Entropia - Gene FOXP2', fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(1.5, 2.05)
    
    # 2. Histograma
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.hist(stats_h['entropias'], bins=50, alpha=0.6, color='#3498db', 
             label='Humano', edgecolor='black')
    ax2.hist(stats_c['entropias'], bins=50, alpha=0.6, color='#e74c3c', 
             label='Chimpanzé', edgecolor='black')
    ax2.axvline(stats_h['media'], color='#3498db', linestyle='--', linewidth=2)
    ax2.axvline(stats_c['media'], color='#e74c3c', linestyle='--', linewidth=2)
    ax2.set_xlabel('Entropia (bits)')
    ax2.set_ylabel('Frequência')
    ax2.set_title('Distribuição', fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3, axis='y')
    
    # 3. Box Plot
    ax3 = fig.add_subplot(gs[1, 1])
    bp = ax3.boxplot([stats_h['entropias'], stats_c['entropias']], 
                     labels=['Humano', 'Chimpanzé'], patch_artist=True, widths=0.6)
    for patch, cor in zip(bp['boxes'], ['#3498db', '#e74c3c']):
        patch.set_facecolor(cor)
        patch.set_alpha(0.6)
    ax3.set_ylabel('Entropia (bits)')
    ax3.set_title('Comparação Estatística', fontweight='bold')
    ax3.grid(True, alpha=0.3, axis='y')
    ax3.set_ylim(1.5, 2.05)
    
    # 4. Tabela
    ax4 = fig.add_subplot(gs[2, :])
    ax4.axis('off')
    dados = [
        ['Métrica', 'Humano', 'Chimpanzé', 'Δ'],
        ['Global', f"{stats_h['global']:.4f}", f"{stats_c['global']:.4f}", 
         f"{abs(stats_h['global']-stats_c['global']):.4f}"],
        ['Média', f"{stats_h['media']:.4f}", f"{stats_c['media']:.4f}", 
         f"{abs(stats_h['media']-stats_c['media']):.4f}"],
        ['Desvio', f"{stats_h['desvio']:.4f}", f"{stats_c['desvio']:.4f}", 
         f"{abs(stats_h['desvio']-stats_c['desvio']):.4f}"]
    ]
    table = ax4.table(cellText=dados, cellLoc='center', loc='center', 
                     colWidths=[0.3, 0.2, 0.2, 0.15])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2)
    for i in range(4):
        table[(0, i)].set_facecolor('#34495e')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    plt.suptitle('Análise de Entropia de Shannon', fontsize=14, fontweight='bold')
    plt.savefig("entropia_completa.png", dpi=300, bbox_inches='tight')
    print("\n✅ Gráfico salvo: entropia_completa.png")

def teste_estatistico(stats_h, stats_c):
    """Testa significância estatística"""
    if not stats_h or not stats_c:
        return
    
    t_stat, p_val = sp_stats.ttest_ind(stats_h['entropias'], stats_c['entropias'])
    
    print("\n" + "=" * 70)
    print("TESTE DE SIGNIFICÂNCIA".center(70))
    print("=" * 70)
    print(f"Teste t: t={t_stat:.4f}, p={p_val:.6f}")
    print(f"Resultado: {'SIGNIFICATIVO' if p_val < 0.05 else 'NÃO SIGNIFICATIVO'} (α=0.05)")
    print("=" * 70 + "\n")

if __name__ == "__main__":
    print("=" * 70)
    print("ANÁLISE DE ENTROPIA DE SHANNON".center(70))
    print("=" * 70 + "\n")
    
    TAMANHO_JANELA = 100
    
    stats_h = analisar_entropia_janela("humano_foxp2.fasta", "Humano", TAMANHO_JANELA)
    stats_c = analisar_entropia_janela("chimp_foxp2.fasta", "Chimpanzé", TAMANHO_JANELA)
    
    if stats_h and stats_c:
        plotar_analise_completa(stats_h, stats_c, TAMANHO_JANELA)
        teste_estatistico(stats_h, stats_c)
        
        print("\n📊 INTERPRETAÇÃO:")
        print(f"   Entropia próxima de 2.0 = Alta complexidade (quase aleatório)")
        print(f"   Entropia ~1.5-1.8 = Complexidade moderada (típico de genes)")
        print(f"   Entropia <1.5 = Baixa complexidade (muitas repetições)\n")