"""
PROJETO LOGOS - Script 02: Análise Básica
==========================================
Estatísticas descritivas das sequências com foco em comparabilidade.

Analisa tanto mRNA completo quanto CDS para demonstrar
por que a comparação deve ser feita sobre a CDS.

REQUISITOS:
    pip install biopython

USO:
    python 02_analise_basica.py
"""

import os
import json
import math
from collections import Counter

# ============================================================================
# FUNÇÕES UTILITÁRIAS (compartilhadas entre scripts)
# ============================================================================

def ler_fasta_simples(caminho):
    """Lê um arquivo FASTA simples (uma sequência) sem Biopython."""
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
    """Remove bases ambíguas, mantém apenas ACGT."""
    limpa = "".join(b for b in seq if b in "ACGT")
    removidos = len(seq) - len(limpa)
    return limpa, removidos


def composicao_bases(seq):
    """Retorna contagens e proporções de cada base."""
    contagens = Counter(seq)
    total = len(seq)
    proporcoes = {b: contagens.get(b, 0) / total for b in "ACGT"}
    return dict(contagens), proporcoes


def conteudo_gc(seq):
    """Calcula conteúdo GC como proporção."""
    c = Counter(seq)
    total = len(seq)
    if total == 0:
        return 0
    return (c.get("G", 0) + c.get("C", 0)) / total


def complexidade_dinucleotideos(seq):
    """
    Conta dinucleotídeos únicos observados / 16 possíveis.
    Nota: com seq > 1kb, quase sempre satura perto de 1.0.
    Mais útil para sequências curtas ou regiões específicas.
    """
    if len(seq) < 2:
        return 0
    dinucs = set(seq[i:i+2] for i in range(len(seq) - 1))
    return len(dinucs) / 16


# ============================================================================
# ANÁLISE PRINCIPAL
# ============================================================================

def analisar_sequencia(caminho, especie, tipo_seq="CDS"):
    """
    Analisa estatísticas básicas de uma sequência FASTA.
    
    Returns:
        dict com estatísticas
    """
    if not os.path.exists(caminho):
        print(f"  ❌ Arquivo não encontrado: {caminho}")
        return None
    
    header, seq = ler_fasta_simples(caminho)
    seq_limpa, removidos = limpar_sequencia(seq)
    
    if len(seq_limpa) == 0:
        print(f"  ⚠️ Sequência vazia em {caminho}")
        return None
    
    gc = conteudo_gc(seq_limpa)
    contagens, proporcoes = composicao_bases(seq_limpa)
    complex_di = complexidade_dinucleotideos(seq_limpa)
    
    return {
        "especie": especie,
        "tipo": tipo_seq,
        "arquivo": caminho,
        "tamanho": len(seq_limpa),
        "bases_removidas": removidos,
        "gc_pct": round(gc * 100, 2),
        "contagens": contagens,
        "proporcoes": {k: round(v, 4) for k, v in proporcoes.items()},
        "complexidade_dinuc": round(complex_di, 3),
    }


def imprimir_tabela(resultados, titulo):
    """Imprime tabela comparativa formatada."""
    validos = [r for r in resultados if r is not None]
    if not validos:
        return
    
    print(f"\n{'='*75}")
    print(f"  {titulo}")
    print(f"{'='*75}")
    print(f"  {'Espécie':<18} | {'Tamanho':>10} | {'GC%':>6} | "
          f"{'A%':>5} | {'C%':>5} | {'G%':>5} | {'T%':>5}")
    print(f"  {'-'*18}-+-{'-'*10}-+-{'-'*6}-+-"
          f"{'-'*5}-+-{'-'*5}-+-{'-'*5}-+-{'-'*5}")
    
    for r in validos:
        p = r["proporcoes"]
        print(f"  {r['especie']:<18} | {r['tamanho']:>10,} | "
              f"{r['gc_pct']:>5.1f}% | "
              f"{p.get('A',0)*100:>4.1f}% | {p.get('C',0)*100:>4.1f}% | "
              f"{p.get('G',0)*100:>4.1f}% | {p.get('T',0)*100:>4.1f}%")
    
    # Diferenças
    if len(validos) >= 2:
        tamanhos = [r["tamanho"] for r in validos]
        max_diff = max(tamanhos) - min(tamanhos)
        max_diff_pct = max_diff / min(tamanhos) * 100
        
        gcs = [r["gc_pct"] for r in validos]
        gc_diff = max(gcs) - min(gcs)
        
        print(f"  {'-'*75}")
        print(f"  Δ Tamanho: {max_diff:,} bp ({max_diff_pct:.1f}%)")
        print(f"  Δ GC%: {gc_diff:.2f} pp")
    
    print(f"{'='*75}")


if __name__ == "__main__":
    print("=" * 75)
    print("  PROJETO LOGOS — ANÁLISE BÁSICA DE SEQUÊNCIAS")
    print("=" * 75)
    
    DATA_DIR = "."
    
    # Definir espécies e arquivos
    especies = [
        ("Humano", "homo_sapiens"),
        ("Chimpanzé", "pan_troglodytes"),
        ("Camundongo", "mus_musculus"),
        ("Gorila", "gorilla_gorilla"),
        ("Macaco Rhesus", "macaca_mulatta"),
    ]
    
    # Análise de CDS
    resultados_cds = []
    print("\n>>> Analisando CDS (coding sequences)...")
    for nome, tag in especies:
        caminho = os.path.join(DATA_DIR, f"{tag}_foxp2_cds.fasta")
        if os.path.exists(caminho):
            r = analisar_sequencia(caminho, nome, "CDS")
            resultados_cds.append(r)
    
    imprimir_tabela(resultados_cds, "COMPARAÇÃO DE CDS — FOXP2")
    
    # Análise de mRNA
    resultados_mrna = []
    print("\n>>> Analisando mRNA completo...")
    for nome, tag in especies:
        caminho = os.path.join(DATA_DIR, f"{tag}_foxp2_mrna.fasta")
        if os.path.exists(caminho):
            r = analisar_sequencia(caminho, nome, "mRNA")
            resultados_mrna.append(r)
    
    imprimir_tabela(resultados_mrna, "COMPARAÇÃO DE mRNA COMPLETO — FOXP2")
    
    # Explicação pedagógica
    cds_validos = [r for r in resultados_cds if r is not None]
    mrna_validos = [r for r in resultados_mrna if r is not None]
    
    if len(cds_validos) >= 2 and len(mrna_validos) >= 2:
        cds_diff = max(r["tamanho"] for r in cds_validos) - min(r["tamanho"] for r in cds_validos)
        cds_diff_pct = cds_diff / min(r["tamanho"] for r in cds_validos) * 100
        
        mrna_diff = max(r["tamanho"] for r in mrna_validos) - min(r["tamanho"] for r in mrna_validos)
        mrna_diff_pct = mrna_diff / min(r["tamanho"] for r in mrna_validos) * 100
        
        print(f"\n{'='*75}")
        print(f"  NOTA METODOLÓGICA")
        print(f"{'='*75}")
        print(f"  Diferença de tamanho nas CDS:  {cds_diff_pct:.1f}%")
        print(f"  Diferença de tamanho nos mRNA: {mrna_diff_pct:.1f}%")
        print(f"")
        print(f"  As CDS codificam a mesma proteína (~715 aa) e diferem pouco")
        print(f"  em tamanho. As diferenças no mRNA refletem variação nas UTRs")
        print(f"  (regiões não traduzidas), que NÃO devem ser interpretadas")
        print(f"  como 'expansão de código'.")
        print(f"{'='*75}")