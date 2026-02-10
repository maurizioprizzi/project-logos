import os
from Bio import Entrez, SeqIO

# --- CONFIGURAÇÃO ---
Entrez.email = "seu_email_aqui@gmail.com"  # <--- MANTENHA SEU EMAIL AQUI

def buscar_e_baixar(organismo, gene, nome_arquivo):
    """
    1. Busca o ID mais recente no NCBI (para não depender de versões fixas).
    2. Baixa a sequência e salva.
    """
    print(f"--- Processando: {organismo} | Gene: {gene} ---")
    
    try:
        # PASSO 1: A Busca Inteligente (Dynamic Search)
        # Query: "Organismo"[Organism] AND "Gene"[Gene Name] AND biomol_mrna[PROP]
        term = f'"{organismo}"[Organism] AND "{gene}"[Gene Name] AND biomol_mrna[PROP]'
        
        print(f"   [...] Buscando ID atualizado para {organismo}...")
        handle_search = Entrez.esearch(db="nucleotide", term=term, retmax=1, sort="relevance")
        record_search = Entrez.read(handle_search)
        handle_search.close()
        
        if not record_search["IdList"]:
            print(f"   [ERRO] Nenhum gene encontrado para {organismo}.")
            return

        # Pega o primeiro ID da lista (o mais relevante/recente)
        id_recente = record_search["IdList"][0]
        print(f"   [!] ID Encontrado: {id_recente}")

        # PASSO 2: O Download (Fetch)
        print(f"   [...] Baixando sequência...")
        handle_fetch = Entrez.efetch(db="nucleotide", id=id_recente, rettype="fasta", retmode="text")
        registro = SeqIO.read(handle_fetch, "fasta")
        handle_fetch.close()
        
        # Salva no disco
        with open(nome_arquivo, "w") as f:
            SeqIO.write(registro, f, "fasta")
            
        print(f"   [SUCESSO] Arquivo salvo: {nome_arquivo}")
        print(f"   Descrição: {registro.description}")
        print(f"   Tamanho: {len(registro.seq)} pares de base\n")
        
    except Exception as e:
        print(f"   [ERRO CRÍTICO]: {e}\n")

if __name__ == "__main__":
    print(">>> PROJETO LOGOS 2.0: COLETA DINÂMICA <<<\n")

    # Humano
    buscar_e_baixar("Homo sapiens", "FOXP2", "humano_foxp2.fasta")

    # Chimpanzé (Agora vai funcionar independente da versão)
    buscar_e_baixar("Pan troglodytes", "FOXP2", "chimp_foxp2.fasta")
    
    print("Processo finalizado.")