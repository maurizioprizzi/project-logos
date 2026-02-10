import os
from time import sleep
from Bio import Entrez, SeqIO
from urllib.error import HTTPError

# --- CONFIGURAÇÃO ---
Entrez.email = "maurizioprizzi@gmail.com"

# Validação de email
if not Entrez.email or "@" not in Entrez.email or "exemplo" in Entrez.email:
    raise ValueError("Configure um email real em Entrez.email!")

def buscar_e_baixar(organismo, gene, nome_arquivo, usar_refseq=True):
    """
    Busca e baixa sequências de genes do NCBI de forma dinâmica.
    
    Args:
        organismo: Nome científico (ex: "Homo sapiens")
        gene: Nome do gene (ex: "FOXP2")
        nome_arquivo: Caminho para salvar o FASTA
        usar_refseq: Se True, prioriza sequências RefSeq
    """
    print(f"--- {organismo} | {gene} ---")
    
    try:
        # Construir query com filtro opcional de RefSeq
        filtro_refseq = ' AND refseq[filter]' if usar_refseq else ''
        term = (f'"{organismo}"[Organism] AND "{gene}"[Gene Name] '
                f'AND biomol_mrna[PROP]{filtro_refseq}')
        
        # PASSO 1: Buscar ID
        print(f"   Buscando ID...")
        handle = Entrez.esearch(db="nucleotide", term=term, retmax=1, sort="relevance")
        resultado = Entrez.read(handle)
        handle.close()
        sleep(0.4)  # Rate limiting
        
        if not resultado["IdList"]:
            print(f"Nenhuma sequência encontrada\n")
            return False
        
        id_seq = resultado["IdList"][0]
        print(f"ID: {id_seq}")
        
        # PASSO 2: Baixar sequência
        print(f"   Baixando...")
        handle = Entrez.efetch(db="nucleotide", id=id_seq, 
                              rettype="fasta", retmode="text")
        registro = SeqIO.read(handle, "fasta")
        handle.close()
        sleep(0.4)
        
        # Validação
        if len(registro.seq) == 0:
            print(f"Sequência vazia!\n")
            return False
        
        # Salvar
        with open(nome_arquivo, "w") as f:
            SeqIO.write(registro, f, "fasta")
        
        print(f"Salvo: {nome_arquivo}")
        print(f"{len(registro.seq)} bp | {registro.description[:60]}...\n")
        return True
        
    except HTTPError as e:
        print(f"Erro de conexão: {e}\n")
        return False
    except RuntimeError as e:
        print(f"Erro no FASTA: {e}\n")
        return False
    except Exception as e:
        print(f"Erro inesperado: {e}\n")
        return False

if __name__ == "__main__":
    print("=== COLETA DE SEQUÊNCIAS NCBI ===\n")
    
    resultados = []
    resultados.append(buscar_e_baixar("Homo sapiens", "FOXP2", "humano_foxp2.fasta"))
    resultados.append(buscar_e_baixar("Pan troglodytes", "FOXP2", "chimp_foxp2.fasta"))
    
    # Relatório final
    sucessos = sum(resultados)
    print(f"Finalizado: {sucessos}/{len(resultados)} sequências baixadas com sucesso.")