"""
PROJETO LOGOS - Script 03: Análise de Complexidade Algorítmica (Zlib)
=====================================================================
Compressibilidade como proxy para Complexidade de Kolmogorov.

CORREÇÕES vs. VERSÃO ANTERIOR:
1. Compara CDS (mesmo tamanho) em vez de mRNAs (tamanhos diferentes)
2. Análise por janelas não-sobrepostas (elimina viés de tamanho do DEFLATE)
3. Inclui teste de viés: comprime subamostras do mesmo tamanho
4. Reporta intervalos de confiança, não apenas médias

REQUISITOS:
    pip install matplotlib numpy scipy

USO:
    python 03_analise_compressao.py
"""

import os
import zlib
import math
import random
import statistics
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats as sp_stats

# ============================================================================
# FUNÇÕES UTILITÁRIAS
# ============================================================================

def ler_fasta_simples(caminho):
    """Lê FASTA simples."""
    header = ""
    seq_parts = []
    with open(caminho, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                header = line[1:]
            else:
                seq_parts.append(line.upper())
    return header, "".join(seq_parts)


def limpar_sequencia(seq):
    """Remove bases ambíguas."""
    return "".join(b for b in seq if b in "ACGT")


# ============================================================================
# COMPRESSÃO
# ============================================================================

def comprimir_sequencia(seq_str):
    """
    Calcula métricas de compressão Zlib para uma sequência de DNA.
    
    Interpretação:
    - Ratio BAIXO = Alta compressão = Mais padrões/redundância
    - Ratio ALTO = Baixa compressão = Mais aleatório
    
    NOTA: O ratio é influenciado pelo tamanho da sequência.
    Sequências mais longas tendem a ter ratios menores porque
    o DEFLATE constrói dicionários mais eficientes.
    """
    bytes_seq = seq_str.encode("utf-8")
    comprimido = zlib.compress(bytes_seq, level=9)
    
    tam_original = len(bytes_seq)
    tam_comprimido = len(comprimido)
    
    # Subtrair overhead do header zlib (~11 bytes)
    ZLIB_OVERHEAD = 11
    tam_dados = max(tam_comprimido - ZLIB_OVERHEAD, 1)
    
    return {
        "original_bytes": tam_original,
        "comprimido_bytes": tam_comprimido,
        "ratio": tam_comprimido / tam_original,
        "ratio_corrigido": tam_dados / tam_original,  # Sem overhead
        "economia_pct": (1 - tam_comprimido / tam_original) * 100,
        "bits_por_base": (tam_dados * 8) / len(seq_str),
    }


def analise_janelas_nao_sobrepostas(seq, tamanho_janela=200):
    """
    Analisa compressibilidade em janelas NÃO SOBREPOSTAS.
    
    Vantagens sobre janelas sobrepostas:
    1. Cada janela é independente (válido para testes estatísticos)
    2. Tamanho fixo elimina viés do DEFLATE
    3. N amostral real, não inflado
    """
    ratios = []
    bits_por_base = []
    
    for i in range(0, len(seq) - tamanho_janela + 1, tamanho_janela):
        janela = seq[i:i + tamanho_janela]
        if len(janela) == tamanho_janela:
            stats = comprimir_sequencia(janela)
            ratios.append(stats["ratio"])
            bits_por_base.append(stats["bits_por_base"])
    
    if len(ratios) < 3:
        return None
    
    return {
        "ratios": np.array(ratios),
        "bits_por_base": np.array(bits_por_base),
        "n_janelas": len(ratios),
        "media_ratio": np.mean(ratios),
        "desvio_ratio": np.std(ratios, ddof=1),
        "media_bpb": np.mean(bits_por_base),
        "desvio_bpb": np.std(bits_por_base, ddof=1),
    }


def teste_vies_tamanho(sequencias, tamanhos_teste=[200, 500, 1000, 2000]):
    """
    Demonstra o viés de tamanho do DEFLATE.
    Para cada sequência, comprime subamostras de tamanhos diferentes.
    
    Se o ratio muda significativamente com o tamanho, 
    comparações globais entre sequências de tamanhos diferentes são inválidas.
    """
    random.seed(42)
    print(f"\n{'='*70}")
    print("  TESTE DE VIÉS DE TAMANHO DO DEFLATE")
    print(f"{'='*70}")
    print(f"  {'Espécie':<15} | ", end="")
    for t in tamanhos_teste:
        print(f" {t:>5} bp |", end="")
    print()
    print(f"  {'-'*15}-+-" + "-+-".join(["-"*8 for _ in tamanhos_teste]) + "-|")
    
    for nome, seq in sequencias:
        print(f"  {nome:<15} | ", end="")
        for t in tamanhos_teste:
            if len(seq) >= t:
                # Média de 5 subamostras aleatórias
                ratios = []
                for _ in range(5):
                    start = random.randint(0, len(seq) - t)
                    sub = seq[start:start + t]
                    stats = comprimir_sequencia(sub)
                    ratios.append(stats["ratio"])
                media = statistics.mean(ratios)
                print(f" {media:.4f}  |", end="")
            else:
                print(f"   N/A   |", end="")
        print()
    
    print(f"\n  ℹ️  Se o ratio diminui com tamanho crescente, comparações")
    print(f"     globais entre sequências de tamanhos diferentes são enviesadas.")
    print(f"     Solução: usar janelas de tamanho fixo.\n")


# ============================================================================
# ANÁLISE PRINCIPAL
# ============================================================================

def analisar_arquivo(caminho, especie, tamanho_janela=200):
    """Analisa compressão de um arquivo FASTA."""
    if not os.path.exists(caminho):
        print(f"  ❌ Não encontrado: {caminho}")
        return None
    
    header, seq = ler_fasta_simples(caminho)
    seq = limpar_sequencia(seq)
    
    if len(seq) < tamanho_janela:
        print(f"  ⚠️ {especie}: sequência muito curta ({len(seq)} bp)")
        return None
    
    # Compressão global
    global_stats = comprimir_sequencia(seq)
    
    # Compressão por janelas
    janelas_stats = analise_janelas_nao_sobrepostas(seq, tamanho_janela)
    
    print(f"  {especie:15s}: {len(seq):>6,} bp | "
          f"Global: {global_stats['ratio']:.4f} | "
          f"Janelas (μ): {janelas_stats['media_ratio']:.4f} ± "
          f"{janelas_stats['desvio_ratio']:.4f} "
          f"(n={janelas_stats['n_janelas']})")
    
    return {
        "especie": especie,
        "tamanho": len(seq),
        "global": global_stats,
        "janelas": janelas_stats,
        "seq": seq,  # Para testes adicionais
    }


def teste_estatistico(resultados, tamanho_janela):
    """
    Teste estatístico robusto entre pares de espécies.
    Usa janelas não-sobrepostas (amostras independentes).
    Reporta tanto p-value quanto tamanho de efeito (Cohen's d).
    """
    print(f"\n{'='*70}")
    print("  TESTES ESTATÍSTICOS (janelas não-sobrepostas)")
    print(f"{'='*70}")
    
    validos = [r for r in resultados if r is not None]
    
    for i in range(len(validos)):
        for j in range(i + 1, len(validos)):
            a = validos[i]
            b = validos[j]
            
            ratios_a = a["janelas"]["ratios"]
            ratios_b = b["janelas"]["ratios"]
            
            # Teste t (Welch, não assume variâncias iguais)
            t_stat, p_val = sp_stats.ttest_ind(ratios_a, ratios_b, equal_var=False)
            
            # Tamanho de efeito: Cohen's d
            na, nb = len(ratios_a), len(ratios_b)
            pooled_std = math.sqrt(
                ((na - 1) * np.std(ratios_a, ddof=1)**2 + (nb - 1) * np.std(ratios_b, ddof=1)**2) 
                / (na + nb - 2)
            )
            cohens_d = abs(np.mean(ratios_a) - np.mean(ratios_b)) / pooled_std if pooled_std > 0 else 0
            
            # Interpretação do Cohen's d
            if cohens_d < 0.2:
                efeito = "NEGLIGÍVEL"
            elif cohens_d < 0.5:
                efeito = "PEQUENO"
            elif cohens_d < 0.8:
                efeito = "MÉDIO"
            else:
                efeito = "GRANDE"
            
            # Mann-Whitney U (não-paramétrico, como controle)
            u_stat, p_mw = sp_stats.mannwhitneyu(ratios_a, ratios_b, alternative="two-sided")
            
            print(f"\n  {a['especie']} vs {b['especie']}")
            print(f"  {'─'*50}")
            print(f"  Teste t (Welch):   t = {t_stat:>8.4f},  p = {p_val:.6f}")
            print(f"  Mann-Whitney U:    U = {u_stat:>8.0f},  p = {p_mw:.6f}")
            print(f"  Cohen's d:         {cohens_d:.4f} → Efeito {efeito}")
            print(f"  N janelas:         {len(ratios_a)} vs {len(ratios_b)}")
            print(f"  Δ média:           {abs(np.mean(ratios_a) - np.mean(ratios_b)):.4f}")


def plotar_resultados(resultados, tamanho_janela):
    """Gera visualizações comparativas."""
    validos = [r for r in resultados if r is not None]
    if len(validos) < 2:
        return
    
    cores = ["#3498db", "#e74c3c", "#2ecc71", "#9b59b6", "#f39c12"]
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    # 1. Box plot dos ratios por janela
    ax1 = axes[0]
    dados_box = [r["janelas"]["ratios"] for r in validos]
    nomes = [r["especie"] for r in validos]
    
    bp = ax1.boxplot(dados_box, labels=nomes, patch_artist=True, widths=0.6)
    for patch, cor in zip(bp["boxes"], cores[:len(validos)]):
        patch.set_facecolor(cor)
        patch.set_alpha(0.6)
    
    ax1.set_ylabel("Ratio de Compressão (Zlib)", fontsize=11)
    ax1.set_title(f"Distribuição por Janelas\n({tamanho_janela} bp, não-sobrepostas)",
                  fontsize=12, fontweight="bold")
    ax1.grid(True, alpha=0.3, axis="y")
    
    # 2. Comparação global vs janelas
    ax2 = axes[1]
    x = np.arange(len(validos))
    width = 0.35
    
    globais = [r["global"]["ratio"] for r in validos]
    janelas = [r["janelas"]["media_ratio"] for r in validos]
    erros = [r["janelas"]["desvio_ratio"] for r in validos]
    
    bars1 = ax2.bar(x - width/2, globais, width, label="Global", 
                    color=[cores[i] for i in range(len(validos))], alpha=0.4,
                    edgecolor="black")
    bars2 = ax2.bar(x + width/2, janelas, width, label=f"Média janelas",
                    color=[cores[i] for i in range(len(validos))], alpha=0.8,
                    edgecolor="black", yerr=erros, capsize=4)
    
    ax2.set_xticks(x)
    ax2.set_xticklabels(nomes)
    ax2.set_ylabel("Ratio de Compressão")
    ax2.set_title("Global vs Janelas\n(viés de tamanho visível no global)", 
                  fontsize=12, fontweight="bold")
    ax2.legend()
    ax2.grid(True, alpha=0.3, axis="y")
    
    # 3. Histograma sobreposto
    ax3 = axes[2]
    for idx, r in enumerate(validos):
        ax3.hist(r["janelas"]["ratios"], bins=20, alpha=0.5,
                color=cores[idx], label=r["especie"], edgecolor="black")
    
    ax3.set_xlabel("Ratio de Compressão")
    ax3.set_ylabel("Frequência")
    ax3.set_title(f"Distribuição de Compressibilidade\n(janelas de {tamanho_janela} bp)",
                  fontsize=12, fontweight="bold")
    ax3.legend()
    ax3.grid(True, alpha=0.3, axis="y")
    
    plt.suptitle("Análise de Complexidade Algorítmica — FOXP2 (CDS)",
                fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig("compressao_foxp2.png", dpi=300, bbox_inches="tight")
    print(f"\n  ✅ Gráfico salvo: compressao_foxp2.png")


# ============================================================================
# MAIN
# ============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("  ANÁLISE DE COMPLEXIDADE ALGORÍTMICA (ZLIB) — FOXP2 CDS")
    print("=" * 70)
    
    DATA_DIR = "."
    TAMANHO_JANELA = 200  # bp por janela
    
    especies = [
        ("Humano", "homo_sapiens"),
        ("Chimpanzé", "pan_troglodytes"),
        ("Camundongo", "mus_musculus"),
        ("Gorila", "gorilla_gorilla"),
        ("Macaco Rhesus", "macaca_mulatta"),
    ]
    
    # Análise principal (CDS)
    print(f"\n>>> Compressão da CDS (janelas de {TAMANHO_JANELA} bp, não-sobrepostas)")
    resultados = []
    sequencias_para_teste = []
    
    for nome, tag in especies:
        caminho = os.path.join(DATA_DIR, f"{tag}_foxp2_cds.fasta")
        r = analisar_arquivo(caminho, nome, TAMANHO_JANELA)
        if r:
            resultados.append(r)
            sequencias_para_teste.append((nome, r["seq"]))
    
    if len(resultados) >= 2:
        # Demonstrar viés de tamanho
        teste_vies_tamanho(sequencias_para_teste)
        
        # Testes estatísticos
        teste_estatistico(resultados, TAMANHO_JANELA)
        
        # Visualizações
        plotar_resultados(resultados, TAMANHO_JANELA)
    else:
        print("\n  ⚠️ Menos de 2 sequências disponíveis.")