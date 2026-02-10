from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

def analisar_sequencia(nome_arquivo, especie):
    """
    Lê um arquivo FASTA e calcula estatísticas básicas:
    - Tamanho total
    - Conteúdo GC (Densidade de informação química)
    - Contagem de Nucleotídeos (A, C, T, G)
    """
    print(f"--- Analisando: {especie} ---")
    
    # Lê a primeira sequência do arquivo
    registro = SeqIO.read(nome_arquivo, "fasta")
    sequencia = registro.seq
    
    # Estatísticas
    tamanho = len(sequencia)
    gc = gc_fraction(sequencia) * 100  # Converte para porcentagem
    
    print(f"Arquivo: {nome_arquivo}")
    print(f"Tamanho: {tamanho} bp")
    print(f"Conteúdo GC: {gc:.2f}%")
    print(f"Composição:")
    print(f"  A: {sequencia.count('A')} ({sequencia.count('A')/tamanho:.1%})")
    print(f"  C: {sequencia.count('C')} ({sequencia.count('C')/tamanho:.1%})")
    print(f"  G: {sequencia.count('G')} ({sequencia.count('G')/tamanho:.1%})")
    print(f"  T: {sequencia.count('T')} ({sequencia.count('T')/tamanho:.1%})")
    print("-" * 30)

if __name__ == "__main__":
    print(">>> PROJETO LOGOS: ESTATÍSTICAS BÁSICAS <<<\n")
    
    analisar_sequencia("humano_foxp2.fasta", "Humano (Homo sapiens)")
    analisar_sequencia("chimp_foxp2.fasta", "Chimpanzé (Pan troglodytes)")