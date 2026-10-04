#!/usr/bin/env python3
import gzip
import re
import urllib.request
import urllib.parse
import json
import time
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def main():
    script_dir = Path(__file__).resolve().parent
    exp_dir = script_dir.parent
    goal_dir = exp_dir.parent.parent
    project_root = goal_dir.parent.parent

    input_csv = project_root / "data/tcga-brca/preprocessed/expression/gen_expr_processed.csv"
    output_fasta = goal_dir / "data/uniprot_proteome.fasta"
    results_dir = exp_dir / "results"
    reports_fig_dir = exp_dir / "reports/figures"
    cache_dir = results_dir / "cache"

    results_dir.mkdir(parents=True, exist_ok=True)
    reports_fig_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    # 1. Ingesta de genes desde el CSV del transcriptoma
    with open(input_csv, "r", encoding="utf-8") as f:
        header = f.readline().strip().split(",")
    genes = [col.strip() for col in header if col.strip() and col.strip() != "patient_id"]
    gene_set = set(genes)
    total_genes = len(genes)
    print(f"[Init] Ingestados {total_genes:,} identificadores HGNC desde {input_csv.name}")

    # 2. Descarga del proteoma de referencia de Homo sapiens (UP000005640)
    uniprot_urls = [
        ("UP000005640_9606.fasta.gz", "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz"),
        ("UP000005640_9606_additional.fasta.gz", "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606_additional.fasta.gz")
    ]

    gene_to_seqs = {}

    for fname, url in uniprot_urls:
        dest_path = cache_dir / fname
        if not dest_path.exists() or dest_path.stat().st_size == 0:
            print(f"[Download] Descargando {fname} desde mirror de UniProt...")
            urllib.request.urlretrieve(url, dest_path)
            print(f"[Download] {fname} descargado ({dest_path.stat().st_size / (1024*1024):.2f} MB)")
        else:
            print(f"[Cache] Utilizando archivo en cache {fname} ({dest_path.stat().st_size / (1024*1024):.2f} MB)")

        # Parseo de FASTA comprimido
        entries_count = 0
        with gzip.open(dest_path, "rt", encoding="utf-8") as gz_f:
            curr_header = None
            curr_seq = []
            for line in gz_f:
                line = line.strip()
                if line.startswith(">"):
                    if curr_header:
                        _process_entry(curr_header, "".join(curr_seq), gene_to_seqs)
                        entries_count += 1
                    curr_header = line
                    curr_seq = []
                else:
                    curr_seq.append(line)
            if curr_header:
                _process_entry(curr_header, "".join(curr_seq), gene_to_seqs)
                entries_count += 1
        print(f"[Parsed] {entries_count:,} entradas procesadas de {fname}")

    # 3. Fallback dirigido para genes no cubiertos directamente en proteoma principal
    unresolved = [g for g in genes if g not in gene_to_seqs]
    print(f"[Coverage] {len(genes) - len(unresolved):,} de {total_genes:,} genes resueltos en proteoma principal. {len(unresolved):,} pendientes de consulta.")

    if unresolved:
        print(f"[Fallback] Iniciando consultas dirigidas para {len(unresolved):,} genes con REST search...")
        recovered_fallback = 0
        chunk_size = 50
        for i in range(0, len(unresolved), chunk_size):
            chunk = unresolved[i:i+chunk_size]
            query_genes = " OR ".join([f"gene_exact:{g}" for g in chunk])
            query = f"(organism_id:9606) AND ({query_genes})"
            params = urllib.parse.urlencode({"query": query, "format": "fasta"})
            req_url = f"https://rest.uniprot.org/uniprotkb/search?{params}"

            for attempt in range(3):
                try:
                    req = urllib.request.Request(req_url, headers={"User-Agent": "WSI_RNA_multimodal_pipeline"})
                    with urllib.request.urlopen(req, timeout=20) as resp:
                        content = resp.read().decode("utf-8")
                        curr_header = None
                        curr_seq = []
                        for line in content.splitlines():
                            line = line.strip()
                            if line.startswith(">"):
                                if curr_header:
                                    _process_entry(curr_header, "".join(curr_seq), gene_to_seqs)
                                curr_header = line
                                curr_seq = []
                            else:
                                curr_seq.append(line)
                        if curr_header:
                            _process_entry(curr_header, "".join(curr_seq), gene_to_seqs)
                    break
                except Exception as e:
                    time.sleep(1.0 * (attempt + 1))
            time.sleep(0.3)
            
            # Chequeo de recuperados en este chunk
            chunk_recovered = sum(1 for g in chunk if g in gene_to_seqs)
            recovered_fallback += chunk_recovered
            if (i // chunk_size) % 10 == 0 or i + chunk_size >= len(unresolved):
                print(f"[Fallback Batch] Procesados {min(i+chunk_size, len(unresolved)):,}/{len(unresolved):,} | Recuperados: {recovered_fallback:,}")

    # 4. Selección determinista de la isoforma más larga (Lmax)
    mapped_records = []
    unmapped_genes = []
    fasta_lines = []

    for gene in sorted(genes):
        if gene in gene_to_seqs and gene_to_seqs[gene]:
            # Ordenar por longitud descendente; desempatar alfabéticamente por accession
            candidates = gene_to_seqs[gene]
            candidates.sort(key=lambda x: (len(x[1]), x[0]), reverse=True)
            best_acc, best_seq = candidates[0]
            best_len = len(best_seq)

            mapped_records.append({
                "gene_id": gene,
                "uniprot_acc": best_acc,
                "length_aa": best_len,
                "status": "mapped",
                "organism": "Homo_sapiens"
            })
            fasta_lines.append(f">{gene}|{best_acc}|{best_len}|Homo_sapiens\n{best_seq}\n")
        else:
            unmapped_genes.append(gene)

    # 5. Guardado del artefacto primario FASTA
    with open(output_fasta, "w", encoding="utf-8") as out_f:
        out_f.writelines(fasta_lines)
    print(f"[Artifact] FASTA escrito exitosamente en {output_fasta} ({len(mapped_records):,} secuencias)")

    # 6. Tablas cuantitativas de resultados
    df_mapped = pd.DataFrame(mapped_records)
    df_mapped.to_csv(results_dir / "mapping_stats.csv", index=False)

    with open(results_dir / "unmapped_genes.txt", "w", encoding="utf-8") as unmapped_f:
        for ug in sorted(unmapped_genes):
            unmapped_f.write(f"{ug}\n")

    n_mapped = len(mapped_records)
    n_unmapped = len(unmapped_genes)
    recovery_pct = (n_mapped / total_genes) * 100
    loss_pct = (n_unmapped / total_genes) * 100

    recovery_summary = [
        {"Metrica": "Total Genes Transcriptoma", "Conteo": total_genes, "Porcentaje": 100.0},
        {"Metrica": "Genes Mapeados a Proteina (Homo sapiens)", "Conteo": n_mapped, "Porcentaje": round(recovery_pct, 2)},
        {"Metrica": "Genes No Codificantes / No Mapeados", "Conteo": n_unmapped, "Porcentaje": round(loss_pct, 2)}
    ]
    df_recovery = pd.DataFrame(recovery_summary)
    df_recovery.to_csv(results_dir / "recovery_table.csv", index=False)
    print(f"[Metrics] Recuperacion: {n_mapped:,}/{total_genes:,} ({recovery_pct:.2f}%) | No mapeados: {n_unmapped:,} ({loss_pct:.2f}%)")

    # 7. Histograma de distribución de longitudes
    lengths = df_mapped["length_aa"].values
    med = float(np.median(lengths))
    mean_val = float(np.mean(lengths))
    q1 = float(np.percentile(lengths, 25))
    q3 = float(np.percentile(lengths, 75))
    p95 = float(np.percentile(lengths, 95))
    max_val = int(np.max(lengths))
    le_1022_pct = float(np.mean(lengths <= 1022) * 100)

    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    bins = np.linspace(0, min(max_val, 4000), 80)
    ax.hist(lengths, bins=bins, color="#2b5c8f", edgecolor="black", alpha=0.8, label="Isoformas más largas (Lmax)")

    ax.axvline(med, color="#d95f02", linestyle="--", linewidth=2, label=f"Mediana ({med:.0f} AA)")
    ax.axvline(mean_val, color="#7570b3", linestyle=":", linewidth=2, label=f"Media ({mean_val:.0f} AA)")
    ax.axvline(1022, color="#1b9e77", linestyle="-", linewidth=2.5, label=f"Límite ESM-2 (1022 AA, {le_1022_pct:.1f}%)")

    ax.set_title("Distribución de Longitudes de Secuencias Proteicas Humanas ($L_{max}$)", fontsize=14, fontweight="bold")
    ax.set_xlabel("Longitud de la secuencia (Aminoácidos - AA)", fontsize=12)
    ax.set_ylabel("Frecuencia (Número de genes)", fontsize=12)
    ax.legend(fontsize=11, loc="upper right")
    ax.grid(True, linestyle=":", alpha=0.5)

    stats_text = (
        f"Genes Mapeados: {n_mapped:,}\n"
        f"Min: {int(np.min(lengths))} AA\n"
        f"Q1: {q1:.0f} AA\n"
        f"Mediana: {med:.0f} AA\n"
        f"Q3: {q3:.0f} AA\n"
        f"P95: {p95:.0f} AA\n"
        f"Máx: {max_val:,} AA\n"
        f"$\\leq$ 1022 AA: {le_1022_pct:.2f}%"
    )
    ax.text(0.68, 0.45, stats_text, transform=ax.transAxes, fontsize=10,
            verticalalignment="top", bbox=dict(boxstyle="round,pad=0.5", facecolor="white", edgecolor="#cccccc", alpha=0.9))

    plot_path = reports_fig_dir / "sequence_length_histogram.png"
    plt.tight_layout()
    plt.savefig(plot_path)
    plt.close()
    print(f"[Plot] Histograma guardado en {plot_path}")
    print(f"[Done] exp_01 completado exitosamente.")

def _process_entry(header, seq, gene_to_seqs):
    # Encabezado UniProt: >db|UniqueIdentifier|EntryName ProteinName OS=Organism OX=Taxon [GN=GeneName] PE=... SV=...
    if "OS=Homo sapiens" not in header and "OX=9606" not in header:
        return
    gn_match = re.search(r"GN=([^\s]+)", header)
    if not gn_match:
        return
    gene_name = gn_match.group(1).strip()
    
    parts = header.split("|")
    acc = parts[1] if len(parts) > 1 else header[1:].split()[0]
    
    if gene_name not in gene_to_seqs:
        gene_to_seqs[gene_name] = []
    gene_to_seqs[gene_name].append((acc, seq))

if __name__ == "__main__":
    main()
