"""
PROJETO LOGOS - Script 04: Análise de Entropia de Shannon
==========================================================
Mede a complexidade informacional do DNA via Entropia de Shannon.

CORREÇÕES vs. VERSÃO ANTERIOR:
1. Janelas NÃO-SOBREPOSTAS para independência estatística
2. Cohen's d (tamanho de efeito) junto com p-value
3. Teste de Levene para homogeneidade de variâncias
4. Análise sobre CDS (não mRNA) para comparação justa
5. Múltiplas espécies como controle filogenético

CONCEITOS:
- H(X) = -Σ p(x) · log₂(p(x))
- Máximo teórico para DNA: 2.0 bits (distribuição uniforme de A, C, G, T)
- Na prática, DNA real nunca atinge 2.0 devido ao viés composicional (GC%)
- Entropia ~1.9 é NORMAL para sequências codificantes

REQUISITOS:
    pip install matplotlib numpy scipy

USO:
    python 04_analise_entropia.py
"""

import os
import math
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
# ENTROPIA DE SHANNON
# ============================================================================

def shannon_entropy(dna_sequence):
    """
    Calcula Entropia de Shannon para uma sequência de DNA.
    
    H(X) = -Σ p(xᵢ) · log₂(p(xᵢ))
    
    Returns:
        float: entropia em bits por base (0 a 2.0)
    """
    seq = "".join(b for b in dna_sequence.upper() if b in "ACGT")
    n = len(seq)
    
    if n == 0:
        return 0.0
    
    contagem = {"A": 0, "C": 0, "G": 0, "T": 0}
    for base in seq:
        contagem[base] += 1
    
    entropia = 0.0
    for count in contagem.values():
        if count > 0:
            p = count / n
            entropia -= p * math.log2(p)
    
    return entropia


def entropia_esperada_por_gc(gc_fraction):
    """
    Calcula a entropia ESPERADA dado um conteúdo GC.
    
    DNA com GC% ≠ 50% terá entropia < 2.0 mesmo sendo "aleatório".
    Isso é crucial para interpretação correta.
    
    Assume: p(G) = p(C) = gc/2, p(A) = p(T) = (1-gc)/2
    """
    if gc_fraction <= 0 or gc_fraction >= 1:
        return 0.0
    
    p_gc = gc_fraction / 2  # p(G) = p(C)
    p_at = (1 - gc_fraction) / 2  # p(A) = p(T)
    
    h = 0.0
    for p in [p_gc, p_gc, p_at, p_at]:
        if p > 0:
            h -= p * math.log2(p)
    
    return h


# ============================================================================
# ANÁLISE POR JANELAS
# ============================================================================

def analisar_entropia_janelas(seq, tamanho_janela=200):
    """
    Calcula entropia em janelas NÃO-SOBREPOSTAS.
    
    Cada janela é independente → válido para testes estatísticos.
    """
    entropias = []
    posicoes = []
    
    for i in range(0, len(seq) - tamanho_janela + 1, tamanho_janela):
        janela = seq[i:i + tamanho_janela]
        if len(janela) == tamanho_janela:
            h = shannon_entropy(janela)
            entropias.append(h)
            posicoes.append(i)
    
    if len(entropias) < 3:
        return None
    
    return {
        "entropias": np.array(entropias),
        "posicoes": np.array(posicoes),
        "n_janelas": len(entropias),
        "media": np.mean(entropias),
        "desvio": np.std(entropias, ddof=1),
        "mediana": np.median(entropias),
        "min": np.min(entropias),
        "max": np.max(entropias),
    }


def analisar_entropia_deslizante(seq, tamanho_janela=100, passo=20):
    """
    Janela deslizante para VISUALIZAÇÃO do perfil ao longo da sequência.
    NÃO usar para testes estatísticos (amostras dependentes).
    """
    entropias = []
    posicoes = []
    
    for i in range(0, len(seq) - tamanho_janela + 1, passo):
        janela = seq[i:i + tamanho_janela]
        h = shannon_entropy(janela)
        entropias.append(h)
        posicoes.append(i + tamanho_janela // 2)  # Centro da janela
    
    return np.array(posicoes), np.array(entropias)


# ============================================================================
# ANÁLISE PRINCIPAL
# ============================================================================

def analisar_arquivo(caminho, especie, tamanho_janela=200):
    """Analisa entropia de um arquivo FASTA."""
    if not os.path.exists(caminho):
        return None
    
    header, seq = ler_fasta_simples(caminho)
    seq = limpar_sequencia(seq)
    
    if len(seq) < tamanho_janela:
        return None
    
    # Entropia global
    h_global = shannon_entropy(seq)
    
    # Entropia esperada pelo GC%
    gc = (seq.count("G") + seq.count("C")) / len(seq)
    h_esperada = entropia_esperada_por_gc(gc)
    
    # Análise por janelas (não-sobrepostas, para estatística)
    janelas = analisar_entropia_janelas(seq, tamanho_janela)
    
    # Perfil deslizante (para visualização)
    pos_vis, ent_vis = analisar_entropia_deslizante(seq, 100, 20)
    
    print(f"  {especie:15s}: {len(seq):>6,} bp | "
          f"H_global={h_global:.4f} | H_esperada(GC)={h_esperada:.4f} | "
          f"H_janelas(μ)={janelas['media']:.4f}±{janelas['desvio']:.4f} "
          f"(n={janelas['n_janelas']})")
    
    return {
        "especie": especie,
        "tamanho": len(seq),
        "gc_pct": gc * 100,
        "h_global": h_global,
        "h_esperada_gc": h_esperada,
        "h_excesso": h_global - h_esperada,  # Positivo = mais complexo que esperado
        "janelas": janelas,
        "perfil_pos": pos_vis,
        "perfil_ent": ent_vis,
        "seq": seq,
    }


def testes_estatisticos(resultados, tamanho_janela):
    """Testes estatísticos robustos com Cohen's d."""
    print(f"\n{'='*70}")
    print("  TESTES ESTATÍSTICOS — ENTROPIA")
    print(f"  (janelas de {tamanho_janela} bp, não-sobrepostas)")
    print(f"{'='*70}")
    
    validos = [r for r in resultados if r is not None]
    
    for i in range(len(validos)):
        for j in range(i + 1, len(validos)):
            a = validos[i]
            b = validos[j]
            
            ea = a["janelas"]["entropias"]
            eb = b["janelas"]["entropias"]
            
            # Teste de Levene (homogeneidade de variâncias)
            lev_stat, lev_p = sp_stats.levene(ea, eb)
            
            # Teste t de Welch (não assume variâncias iguais)
            t_stat, p_val = sp_stats.ttest_ind(ea, eb, equal_var=False)
            
            # Cohen's d
            na, nb = len(ea), len(eb)
            pooled_std = math.sqrt(
                ((na - 1) * np.std(ea, ddof=1)**2 + (nb - 1) * np.std(eb, ddof=1)**2) 
                / (na + nb - 2)
            )
            cohens_d = abs(np.mean(ea) - np.mean(eb)) / pooled_std if pooled_std > 0 else 0
            
            if cohens_d < 0.2:
                efeito = "NEGLIGÍVEL"
            elif cohens_d < 0.5:
                efeito = "PEQUENO"
            elif cohens_d < 0.8:
                efeito = "MÉDIO"
            else:
                efeito = "GRANDE"
            
            # Intervalo de confiança da diferença de médias
            diff = np.mean(ea) - np.mean(eb)
            se_diff = math.sqrt(np.var(ea, ddof=1)/len(ea) + np.var(eb, ddof=1)/len(eb))
            ci_95 = (diff - 1.96 * se_diff, diff + 1.96 * se_diff)
            
            print(f"\n  {a['especie']} vs {b['especie']}")
            print(f"  {'─'*55}")
            print(f"  Médias:         {np.mean(ea):.4f} vs {np.mean(eb):.4f}")
            print(f"  Diferença (Δ):  {diff:+.4f} bits")
            print(f"  IC 95%:         [{ci_95[0]:+.4f}, {ci_95[1]:+.4f}]")
            print(f"  Teste t:        t = {t_stat:>8.4f},  p = {p_val:.6f}")
            print(f"  Cohen's d:      {cohens_d:.4f} → Efeito {efeito}")
            print(f"  Levene:         W = {lev_stat:>8.4f},  p = {lev_p:.4f} "
                  f"({'variâncias ≠' if lev_p < 0.05 else 'variâncias ≈'})")
            print(f"  N janelas:      {len(ea)} vs {len(eb)}")


def tabela_contexto_gc(resultados):
    """
    Mostra que a entropia observada deve ser contextualizada pelo GC%.
    """
    validos = [r for r in resultados if r is not None]
    if not validos:
        return
    
    print(f"\n{'='*70}")
    print("  CONTEXTO: ENTROPIA vs GC%")
    print(f"{'='*70}")
    print(f"  {'Espécie':<15} | {'GC%':>5} | {'H_obs':>6} | "
          f"{'H_esp(GC)':>9} | {'Excesso':>8} | Interp.")
    print(f"  {'-'*15}-+-{'-'*5}-+-{'-'*6}-+-{'-'*9}-+-{'-'*8}-+--------")
    
    for r in validos:
        exc = r["h_excesso"]
        if abs(exc) < 0.01:
            interp = "≈ esperado"
        elif exc > 0:
            interp = "mais variado"
        else:
            interp = "mais restrito"
        
        print(f"  {r['especie']:<15} | {r['gc_pct']:>4.1f}% | "
              f"{r['h_global']:>6.4f} | {r['h_esperada_gc']:>9.4f} | "
              f"{exc:>+8.4f} | {interp}")
    
    print(f"\n  ℹ️  H_esperada(GC) é a entropia que uma sequência ALEATÓRIA")
    print(f"     com o mesmo GC% teria. Diferenças da esperada são mais")
    print(f"     informativas que diferenças absolutas entre espécies.")
    print(f"{'='*70}")


# ============================================================================
# VISUALIZAÇÕES
# ============================================================================

def plotar_analise(resultados, tamanho_janela):
    """Gera visualizações completas."""
    validos = [r for r in resultados if r is not None]
    if len(validos) < 2:
        return
    
    cores = ["#3498db", "#e74c3c", "#2ecc71", "#9b59b6", "#f39c12"]
    
    fig = plt.figure(figsize=(18, 12))
    gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.3)
    
    # 1. Perfil ao longo da sequência (janela deslizante, só visual)
    ax1 = fig.add_subplot(gs[0, :])
    for idx, r in enumerate(validos):
        ax1.plot(r["perfil_pos"], r["perfil_ent"],
                label=r["especie"], color=cores[idx], alpha=0.7, linewidth=1.2)
    
    ax1.axhline(y=2.0, color="black", linestyle="--", linewidth=1,
               label="Máximo teórico", alpha=0.5)
    ax1.set_xlabel("Posição na CDS (bp)")
    ax1.set_ylabel("Entropia (bits)")
    ax1.set_title("Perfil de Entropia — FOXP2 CDS\n(janela deslizante 100 bp, "
                  "passo 20 bp — apenas visualização)", fontweight="bold")
    ax1.legend(loc="lower right")
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(1.4, 2.05)
    
    # 2. Box plot (janelas não-sobrepostas)
    ax2 = fig.add_subplot(gs[1, 0])
    dados_box = [r["janelas"]["entropias"] for r in validos]
    nomes = [r["especie"] for r in validos]
    
    bp = ax2.boxplot(dados_box, labels=nomes, patch_artist=True, widths=0.6)
    for patch, cor in zip(bp["boxes"], cores[:len(validos)]):
        patch.set_facecolor(cor)
        patch.set_alpha(0.6)
    
    ax2.set_ylabel("Entropia (bits)")
    ax2.set_title(f"Distribuição\n({tamanho_janela} bp, não-sobrepostas)",
                  fontweight="bold")
    ax2.grid(True, alpha=0.3, axis="y")
    ax2.set_ylim(1.5, 2.05)
    
    # 3. Entropia observada vs esperada (contextualização por GC%)
    ax3 = fig.add_subplot(gs[1, 1])
    x_pos = np.arange(len(validos))
    width = 0.35
    
    h_obs = [r["h_global"] for r in validos]
    h_esp = [r["h_esperada_gc"] for r in validos]
    
    ax3.bar(x_pos - width/2, h_obs, width, label="Observada",
           color=[cores[i] for i in range(len(validos))], alpha=0.8, edgecolor="black")
    ax3.bar(x_pos + width/2, h_esp, width, label="Esperada (GC%)",
           color="lightgray", edgecolor="black")
    
    ax3.set_xticks(x_pos)
    ax3.set_xticklabels(nomes, rotation=15, ha="right")
    ax3.set_ylabel("Entropia (bits)")
    ax3.set_title("Observada vs Esperada (GC%)", fontweight="bold")
    ax3.legend()
    ax3.grid(True, alpha=0.3, axis="y")
    ax3.set_ylim(1.8, 2.0)
    
    # 4. Histograma sobreposto
    ax4 = fig.add_subplot(gs[1, 2])
    for idx, r in enumerate(validos):
        ax4.hist(r["janelas"]["entropias"], bins=15, alpha=0.5,
                color=cores[idx], label=r["especie"], edgecolor="black")
    
    ax4.set_xlabel("Entropia (bits)")
    ax4.set_ylabel("Frequência")
    ax4.set_title(f"Histograma de Entropias\n(janelas {tamanho_janela} bp)",
                  fontweight="bold")
    ax4.legend()
    ax4.grid(True, alpha=0.3, axis="y")
    
    plt.suptitle("Análise de Entropia de Shannon — FOXP2 (CDS)",
                fontsize=14, fontweight="bold", y=1.01)
    plt.savefig("entropia_foxp2.png", dpi=300, bbox_inches="tight")
    print(f"\n  ✅ Gráfico salvo: entropia_foxp2.png")


# ============================================================================
# MAIN
# ============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("  ANÁLISE DE ENTROPIA DE SHANNON — FOXP2 CDS")
    print("=" * 70)
    
    DATA_DIR = "."
    TAMANHO_JANELA = 200
    
    especies = [
        ("Humano", "homo_sapiens"),
        ("Chimpanzé", "pan_troglodytes"),
        ("Camundongo", "mus_musculus"),
        ("Gorila", "gorilla_gorilla"),
        ("Macaco Rhesus", "macaca_mulatta"),
    ]
    
    print(f"\n>>> Entropia da CDS (janelas de {TAMANHO_JANELA} bp, não-sobrepostas)")
    resultados = []
    
    for nome, tag in especies:
        caminho = os.path.join(DATA_DIR, f"{tag}_foxp2_cds.fasta")
        r = analisar_arquivo(caminho, nome, TAMANHO_JANELA)
        resultados.append(r)
    
    validos = [r for r in resultados if r is not None]
    
    if len(validos) >= 2:
        tabela_contexto_gc(resultados)
        testes_estatisticos(resultados, TAMANHO_JANELA)
        plotar_analise(resultados, TAMANHO_JANELA)
    else:
        print("\n  ⚠️ Menos de 2 sequências disponíveis.")