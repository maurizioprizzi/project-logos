"""
PROJETO LOGOS - Script 01: Coleta de Dados Controlada
=====================================================
Busca sequências FOXP2 do NCBI usando accession numbers específicos
para garantir comparação entre isoformas equivalentes.

ESTRATÉGIA:
- Usa accession numbers fixos (não busca por "relevância")
- Baixa tanto o mRNA completo quanto a CDS (coding sequence) separadamente
- Inclui múltiplas espécies como controle filogenético
- Valida que as proteínas codificadas têm tamanhos comparáveis

REQUISITOS:
    pip install biopython

USO:
    python 01_coleta_dados.py
"""

import os
import json
from time import sleep
from Bio import Entrez, SeqIO
from urllib.error import HTTPError

# --- CONFIGURAÇÃO ---
Entrez.email = "maurizioprizzi@gmail.com"  # OBRIGATÓRIO pelo NCBI

# Diretório de saída (mesmo diretório dos scripts)
OUTPUT_DIR = "."

# ============================================================================
# ACCESSION NUMBERS CONTROLADOS
# ============================================================================
# Cada entrada usa o accession number específico da isoforma canônica
# para garantir comparação justa entre espécies.
#
# Critérios de seleção:
# - Humano: NM_014491 = transcript variant 1 (isoforma I, 715 aa, canônica)
# - Chimpanzé: NM_001009020 = única isoforma NM_ curada (715 aa)
# - Camundongo: NM_053242 = isoforma canônica (711 aa)
# - Gorila/Macaco: XM_ (modelos preditos, melhor disponível)
#
# NOTA: A proteína FOXP2 tem 715 aa em humano/chimp/gorila, 711 em mouse.
# As diferenças proteicas são: 2 aa entre humano e chimp, 3 entre humano e mouse.
# ============================================================================

SEQUENCIAS = {
    "Homo_sapiens": {
        "accession": "NM_014491",
        "nome_comum": "Humano",
        "cor": "#3498db",
        "nota": "Isoforma I canônica (715 aa), RefSeq Select"
    },
    "Pan_troglodytes": {
        "accession": "NM_001009020",
        "nome_comum": "Chimpanzé",
        "cor": "#e74c3c",
        "nota": "Única isoforma NM_ curada (715 aa)"
    },
    "Mus_musculus": {
        "accession": "NM_053242",
        "nome_comum": "Camundongo",
        "cor": "#2ecc71",
        "nota": "Isoforma canônica (711 aa)"
    },
}

# Espécies adicionais — buscadas dinamicamente porque accession XM_ mudam
# entre versões de anotação do genoma. O script busca pelo gene FOXP2 e
# filtra pela CDS com proteína de ~715 aa.
SEQUENCIAS_EXTRAS_BUSCA = {
    "Gorilla_gorilla": {
        "nome_comum": "Gorila",
        "cor": "#9b59b6",
        "gene": "FOXP2",
        "proteina_esperada": (700, 720),  # Range de aa esperado
    },
    "Macaca_mulatta": {
        "nome_comum": "Macaco Rhesus",
        "cor": "#f39c12",
        "gene": "FOXP2",
        "proteina_esperada": (700, 720),
    },
}


def baixar_sequencia(accession, especie, nome_comum):
    """
    Baixa uma sequência específica do NCBI pelo accession number.
    Salva tanto o mRNA completo (GenBank) quanto a CDS extraída.
    
    Returns:
        dict com metadados ou None se falhou
    """
    print(f"\n{'='*60}")
    print(f"  {nome_comum} ({especie})")
    print(f"  Accession: {accession}")
    print(f"{'='*60}")
    
    try:
        # Baixar em formato GenBank (tem anotação de CDS)
        print(f"  Baixando GenBank...")
        handle = Entrez.efetch(
            db="nucleotide",
            id=accession,
            rettype="gb",
            retmode="text"
        )
        registro = SeqIO.read(handle, "genbank")
        handle.close()
        sleep(0.5)
        
        seq_mrna = str(registro.seq).upper()
        print(f"  mRNA completo: {len(seq_mrna):,} bp")
        print(f"  Descrição: {registro.description[:70]}...")
        
        # Extrair CDS
        cds_seq = None
        proteina_len = None
        
        for feature in registro.features:
            if feature.type == "CDS":
                cds_seq = str(feature.extract(registro.seq)).upper()
                # Verificar se há tradução anotada
                if "translation" in feature.qualifiers:
                    proteina_len = len(feature.qualifiers["translation"][0])
                break
        
        if cds_seq:
            print(f"  CDS extraída: {len(cds_seq):,} bp")
            if proteina_len:
                print(f"  Proteína: {proteina_len} aminoácidos")
        else:
            print(f"  ⚠️  CDS não encontrada no GenBank, usando mRNA completo")
            cds_seq = seq_mrna
        
        # Salvar arquivos
        tag = especie.lower().replace(" ", "_")
        
        # mRNA completo (FASTA)
        mrna_path = os.path.join(OUTPUT_DIR, f"{tag}_foxp2_mrna.fasta")
        with open(mrna_path, "w") as f:
            f.write(f">{accession} {especie} FOXP2 mRNA completo\n")
            for i in range(0, len(seq_mrna), 70):
                f.write(seq_mrna[i:i+70] + "\n")
        
        # CDS apenas (FASTA)
        cds_path = os.path.join(OUTPUT_DIR, f"{tag}_foxp2_cds.fasta")
        with open(cds_path, "w") as f:
            f.write(f">{accession}_CDS {especie} FOXP2 CDS\n")
            for i in range(0, len(cds_seq), 70):
                f.write(cds_seq[i:i+70] + "\n")
        
        print(f"  ✓ Salvos: {mrna_path}, {cds_path}")
        
        return {
            "especie": especie,
            "nome_comum": nome_comum,
            "accession": accession,
            "mrna_tamanho": len(seq_mrna),
            "cds_tamanho": len(cds_seq),
            "proteina_aa": proteina_len,
            "mrna_path": mrna_path,
            "cds_path": cds_path,
            "descricao": registro.description
        }
        
    except HTTPError as e:
        print(f"  ❌ Erro HTTP: {e}")
        print(f"     Verifique o accession number: {accession}")
        return None
    except Exception as e:
        print(f"  ❌ Erro: {e}")
        return None


def validar_comparabilidade(resultados):
    """
    Verifica se as sequências baixadas são comparáveis.
    Alerta sobre diferenças de tamanho inesperadas.
    """
    print(f"\n{'='*60}")
    print("  VALIDAÇÃO DE COMPARABILIDADE")
    print(f"{'='*60}")
    
    validos = [r for r in resultados if r is not None]
    
    if len(validos) < 2:
        print("  ⚠️  Menos de 2 sequências válidas. Comparação impossível.")
        return False
    
    # Comparar tamanhos de CDS
    cds_sizes = [(r["nome_comum"], r["cds_tamanho"]) for r in validos]
    prot_sizes = [(r["nome_comum"], r["proteina_aa"]) for r in validos if r["proteina_aa"]]
    
    print(f"\n  Tamanhos de CDS:")
    for nome, tam in cds_sizes:
        print(f"    {nome:15s}: {tam:,} bp")
    
    if prot_sizes:
        print(f"\n  Tamanhos de proteína:")
        for nome, tam in prot_sizes:
            print(f"    {nome:15s}: {tam} aa")
    
    # Verificar se as CDS têm tamanhos razoavelmente próximos
    tamanhos = [t for _, t in cds_sizes]
    max_diff = max(tamanhos) - min(tamanhos)
    max_diff_pct = max_diff / min(tamanhos) * 100
    
    if max_diff_pct > 10:
        print(f"\n  ⚠️  ALERTA: Diferença de CDS > 10% ({max_diff_pct:.1f}%)")
        print(f"     Verifique se as isoformas são realmente equivalentes!")
    else:
        print(f"\n  ✓ Diferença de CDS: {max_diff_pct:.1f}% (aceitável)")
    
    # Comparar tamanhos de mRNA (aqui diferenças são esperadas por UTR)
    mrna_sizes = [(r["nome_comum"], r["mrna_tamanho"]) for r in validos]
    print(f"\n  Tamanhos de mRNA (UTR incluída):")
    for nome, tam in mrna_sizes:
        print(f"    {nome:15s}: {tam:,} bp")
    
    mrna_tams = [t for _, t in mrna_sizes]
    mrna_diff_pct = (max(mrna_tams) - min(mrna_tams)) / min(mrna_tams) * 100
    if mrna_diff_pct > 50:
        print(f"\n  ℹ️  mRNA diferem em {mrna_diff_pct:.0f}% — esperado (UTRs variam)")
        print(f"     Use a CDS para comparações informacionais justas.")
    
    return True


def buscar_foxp2_dinamico(organismo, gene, nome_comum, proteina_range):
    """
    Busca FOXP2 dinamicamente para espécies sem accession NM_ estável.
    Filtra pela CDS com proteína no range esperado (~715 aa).
    
    Returns:
        dict com metadados ou None
    """
    print(f"\n{'='*60}")
    print(f"  {nome_comum} ({organismo}) — BUSCA DINÂMICA")
    print(f"{'='*60}")
    
    try:
        # Buscar por gene + organismo, mRNA RefSeq
        term = (f'"{organismo}"[Organism] AND "{gene}"[Gene Name] '
                f'AND biomol_mrna[PROP] AND refseq[filter]')
        
        print(f"  Buscando...")
        handle = Entrez.esearch(db="nucleotide", term=term, retmax=10, sort="relevance")
        resultado = Entrez.read(handle)
        handle.close()
        sleep(0.5)
        
        if not resultado["IdList"]:
            print(f"  ❌ Nenhuma sequência encontrada")
            return None
        
        print(f"  {len(resultado['IdList'])} candidatos encontrados")
        
        # Testar cada candidato até achar um com proteína no range correto
        for id_seq in resultado["IdList"]:
            handle = Entrez.efetch(db="nucleotide", id=id_seq,
                                  rettype="gb", retmode="text")
            registro = SeqIO.read(handle, "genbank")
            handle.close()
            sleep(0.5)
            
            # Verificar CDS
            for feature in registro.features:
                if feature.type == "CDS" and "translation" in feature.qualifiers:
                    prot_len = len(feature.qualifiers["translation"][0])
                    if proteina_range[0] <= prot_len <= proteina_range[1]:
                        # Encontramos a isoforma correta!
                        accession = registro.id
                        seq_mrna = str(registro.seq).upper()
                        cds_seq = str(feature.extract(registro.seq)).upper()
                        
                        print(f"  ✓ Encontrado: {accession}")
                        print(f"  mRNA: {len(seq_mrna):,} bp | CDS: {len(cds_seq):,} bp | "
                              f"Proteína: {prot_len} aa")
                        print(f"  {registro.description[:70]}...")
                        
                        # Salvar
                        tag = organismo.lower().replace(" ", "_")
                        
                        mrna_path = os.path.join(OUTPUT_DIR, f"{tag}_foxp2_mrna.fasta")
                        with open(mrna_path, "w") as f:
                            f.write(f">{accession} {organismo} FOXP2 mRNA\n")
                            for i in range(0, len(seq_mrna), 70):
                                f.write(seq_mrna[i:i+70] + "\n")
                        
                        cds_path = os.path.join(OUTPUT_DIR, f"{tag}_foxp2_cds.fasta")
                        with open(cds_path, "w") as f:
                            f.write(f">{accession}_CDS {organismo} FOXP2 CDS\n")
                            for i in range(0, len(cds_seq), 70):
                                f.write(cds_seq[i:i+70] + "\n")
                        
                        print(f"  ✓ Salvos: {mrna_path}, {cds_path}")
                        
                        return {
                            "especie": organismo,
                            "nome_comum": nome_comum,
                            "accession": accession,
                            "mrna_tamanho": len(seq_mrna),
                            "cds_tamanho": len(cds_seq),
                            "proteina_aa": prot_len,
                            "mrna_path": mrna_path,
                            "cds_path": cds_path,
                            "descricao": registro.description
                        }
            
            print(f"    {registro.id}: proteína fora do range, pulando...")
        
        print(f"  ❌ Nenhuma isoforma com ~715 aa encontrada")
        return None
        
    except Exception as e:
        print(f"  ❌ Erro: {e}")
        return None


if __name__ == "__main__":
    print("=" * 60)
    print("  PROJETO LOGOS — COLETA DE DADOS CONTROLADA")
    print("  Gene: FOXP2 | Fonte: NCBI RefSeq")
    print("=" * 60)
    
    resultados = []
    
    # Espécies principais (NM_ curadas — accession fixo)
    print("\n>>> ESPÉCIES PRINCIPAIS (RefSeq curadas)")
    for especie, info in SEQUENCIAS.items():
        r = baixar_sequencia(info["accession"], especie, info["nome_comum"])
        if r:
            resultados.append(r)
    
    # Espécies extras (busca dinâmica com validação de proteína)
    print("\n>>> ESPÉCIES ADICIONAIS (busca dinâmica)")
    for especie, info in SEQUENCIAS_EXTRAS_BUSCA.items():
        r = buscar_foxp2_dinamico(
            especie, info["gene"], info["nome_comum"], info["proteina_esperada"]
        )
        if r:
            resultados.append(r)
    
    # Validação
    validar_comparabilidade(resultados)
    
    # Salvar metadados
    meta_path = os.path.join(OUTPUT_DIR, "metadados.json")
    with open(meta_path, "w") as f:
        json.dump(resultados, f, indent=2, ensure_ascii=False)
    print(f"\n  💾 Metadados salvos: {meta_path}")
    
    # Resumo
    print(f"\n{'='*60}")
    print(f"  COLETA FINALIZADA: {len(resultados)} sequências")
    print(f"  Diretório: {OUTPUT_DIR}/")
    print(f"{'='*60}")