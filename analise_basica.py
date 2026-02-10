import os
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from collections import Counter

def analisar_sequencia(nome_arquivo, especie, verbose=True):
    """
    Analisa estatísticas de uma sequência FASTA.
    
    Retorna:
        dict com estatísticas ou None se houver erro
    """
    if verbose:
        print(f"--- {especie} ---")
    
    # Validações
    if not os.path.exists(nome_arquivo):
        print(f"❌ Arquivo não encontrado: {nome_arquivo}\n")
        return None
    
    try:
        registro = SeqIO.read(nome_arquivo, "fasta")
        seq = str(registro.seq).upper()  # Normaliza para maiúsculas
    except Exception as e:
        print(f"❌ Erro ao ler arquivo: {e}\n")
        return None
    
    # Estatísticas básicas
    tamanho = len(seq)
    if tamanho == 0:
        print(f"⚠️ Sequência vazia em {nome_arquivo}\n")
        return None
    
    gc = gc_fraction(seq) * 100
    
    # Composição nucleotídica
    contagens = {
        'A': seq.count('A'),
        'C': seq.count('C'),
        'G': seq.count('G'),
        'T': seq.count('T')
    }
    outros = tamanho - sum(contagens.values())
    
    # Métricas avançadas
    dinucs = [seq[i:i+2] for i in range(len(seq)-1) if i+1 < len(seq)]
    complexidade = len(set(dinucs)) / 16 if dinucs else 0
    
    # Display
    if verbose:
        print(f"Arquivo: {nome_arquivo}")
        print(f"Descrição: {registro.description[:70]}...")
        print(f"Tamanho: {tamanho:,} bp")
        print(f"Conteúdo GC: {gc:.2f}%")
        print(f"Complexidade: {complexidade:.1%}")
        print(f"\nComposição:")
        for base in 'ACGT':
            count = contagens[base]
            print(f"  {base}: {count:6,} ({count/tamanho:6.1%})")
        
        if outros > 0:
            print(f"  N/Outros: {outros:6,} ({outros/tamanho:6.1%})")
        
        print("-" * 50)
    
    # Retorna dados estruturados
    return {
        'especie': especie,
        'arquivo': nome_arquivo,
        'descricao': registro.description,
        'tamanho': tamanho,
        'gc_percentual': round(gc, 2),
        'complexidade': round(complexidade, 3),
        **contagens,
        'outros': outros
    }

def comparar_especies(stats_lista):
    """Exibe tabela comparativa entre espécies."""
    if not stats_lista or all(s is None for s in stats_lista):
        print("⚠️ Nenhuma estatística disponível para comparação\n")
        return
    
    # Filtra None
    stats_validas = [s for s in stats_lista if s is not None]
    
    print("\n" + "=" * 70)
    print("COMPARAÇÃO ENTRE ESPÉCIES".center(70))
    print("=" * 70)
    print(f"{'Espécie':<25} | {'Tamanho (bp)':>12} | {'GC%':>6} | {'Complex.':>8}")
    print("-" * 70)
    
    for s in stats_validas:
        print(f"{s['especie']:<25} | {s['tamanho']:>12,} | {s['gc_percentual']:>5.1f}% | {s['complexidade']:>7.1%}")
    
    # Diferenças
    if len(stats_validas) == 2:
        diff_tamanho = abs(stats_validas[0]['tamanho'] - stats_validas[1]['tamanho'])
        diff_gc = abs(stats_validas[0]['gc_percentual'] - stats_validas[1]['gc_percentual'])
        
        print("-" * 70)
        print(f"Δ Tamanho: {diff_tamanho:,} bp ({diff_tamanho/max(s['tamanho'] for s in stats_validas):.1%})")
        print(f"Δ GC%: {diff_gc:.2f} pontos percentuais")
    
    print("=" * 70 + "\n")

if __name__ == "__main__":
    print("=" * 70)
    print("ANÁLISE DE SEQUÊNCIAS - PROJETO LOGOS".center(70))
    print("=" * 70 + "\n")
    
    resultados = []
    
    # Analisa cada espécie
    resultados.append(
        analisar_sequencia("humano_foxp2.fasta", "Humano (Homo sapiens)")
    )
    resultados.append(
        analisar_sequencia("chimp_foxp2.fasta", "Chimpanzé (Pan troglodytes)")
    )
    
    # Comparação final
    comparar_especies(resultados)
    
    # Exemplo de uso dos dados retornados
    if all(resultados):
        print("💾 Dados disponíveis para análises adicionais")
        print(f"   Exemplo: {resultados[0]['especie']} tem {resultados[0]['A']:,} adeninas\n")