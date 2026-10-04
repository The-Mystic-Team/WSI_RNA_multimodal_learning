#!/usr/bin/env python3
import hashlib
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

def parse_fasta(fasta_path):
    records = []
    current_header = None
    current_seq = []
    with open(fasta_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header:
                    seq_str = "".join(current_seq)
                    records.append((current_header, seq_str))
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            records.append((current_header, "".join(current_seq)))
    return records

def format_fasta_seq(seq, width=80):
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

def main():
    script_dir = Path(__file__).resolve().parent
    exp_dir = script_dir.parent
    goal_dir = exp_dir.parent.parent
    project_root = goal_dir.parent.parent

    input_csv = project_root / "data/tcga-brca/preprocessed/expression/gen_expr_processed.csv"
    input_fasta = goal_dir / "data/uniprot_proteome.fasta"
    input_stats = exp_dir.parent / "exp_01/results/mapping_stats.csv"

    output_fasta = goal_dir / "data/filtered_proteome.fasta"
    results_dir = exp_dir / "results"
    reports_fig_dir = exp_dir / "reports/figures"

    results_dir.mkdir(parents=True, exist_ok=True)
    reports_fig_dir.mkdir(parents=True, exist_ok=True)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    # 1. Ingesta de genes de referencia desde el transcriptoma
    with open(input_csv, "r", encoding="utf-8") as f:
        header_cols = f.readline().strip().split(",")
    transcriptome_genes = [c.strip() for c in header_cols if c.strip() and c.strip() != "patient_id"]
    n_total = len(transcriptome_genes)
    print(f"[Phase A] Ingestados {n_total:,} genes desde el transcriptoma")

    # 2. Carga y parseo del proteoma de exp_01
    fasta_records = parse_fasta(input_fasta)
    n_mapped = len(fasta_records)
    print(f"[Phase A] Cargadas {n_mapped:,} secuencias de Homo sapiens desde {input_fasta.name}")

    # 3. Filtrado determinístico L <= 1022 AA
    passed_records = []
    excluded_by_length = []

    for header, seq in fasta_records:
        parts = header.split("|")
        gene_id = parts[0]
        uniprot_acc = parts[1] if len(parts) > 1 else "Unknown"
        length_aa = len(seq)

        if length_aa <= 1022:
            passed_records.append((gene_id, uniprot_acc, length_aa, seq))
        else:
            excluded_by_length.append({
                "gene_id": gene_id,
                "uniprot_acc": uniprot_acc,
                "length_aa": length_aa,
                "exclusion_reason": "exceeds_length_threshold"
            })

    n_passed = len(passed_records)
    n_exceeded = len(excluded_by_length)
    print(f"[Phase B] Filtro L <= 1022 AA completado: {n_passed:,} aprobados | {n_exceeded:,} exceden longitud")

    # 4. Escritura de filtered_proteome.fasta en data/ del goal
    sha256 = hashlib.sha256()
    with open(output_fasta, "w", encoding="utf-8") as f:
        for gene_id, acc, length_aa, seq in passed_records:
            entry = f">{gene_id}|{acc}|{length_aa}|Homo_sapiens\n{format_fasta_seq(seq)}\n"
            f.write(entry)
            sha256.update(entry.encode("utf-8"))

    fasta_hash = sha256.hexdigest()
    with open(results_dir / "filtered_proteome.sha256", "w", encoding="utf-8") as f:
        f.write(f"{fasta_hash}  filtered_proteome.fasta\n")
    print(f"[Phase B] Artefacto generado: {output_fasta.name} (SHA-256: {fasta_hash[:12]}...)")

    # 5. Carga de genes no mapeados/no codificantes
    unmapped_file = exp_dir.parent / "exp_01/results/unmapped_genes.txt"
    if unmapped_file.exists():
        with open(unmapped_file, "r", encoding="utf-8") as f:
            unmapped_genes = [line.strip() for line in f if line.strip()]
    else:
        mapped_set = {rec[0] for rec in fasta_records}
        unmapped_genes = [g for g in transcriptome_genes if g not in mapped_set]

    excluded_non_coding = [
        {
            "gene_id": g,
            "uniprot_acc": "None",
            "length_aa": 0,
            "exclusion_reason": "non_coding_or_unmapped"
        }
        for g in unmapped_genes
    ]
    n_non_coding = len(excluded_non_coding)

    # Consolidación de exclusiones
    all_exclusions = pd.DataFrame(excluded_non_coding + excluded_by_length)
    all_exclusions.sort_values(by=["exclusion_reason", "length_aa", "gene_id"], ascending=[True, False, True], inplace=True)
    all_exclusions.to_csv(results_dir / "excluded_genes_summary.csv", index=False)
    print(f"[Phase C] Detalle de exclusiones consolidado: {len(all_exclusions):,} genes descartados")

    # 6. Tablas de métricas cuantitativas
    pct_transcriptome_passed = (n_passed / n_total) * 100
    pct_proteome_passed = (n_passed / n_mapped) * 100
    pct_mapped = (n_mapped / n_total) * 100
    pct_non_coding = (n_non_coding / n_total) * 100
    pct_exceeded = (n_exceeded / n_total) * 100

    recovery_df = pd.DataFrame([
        {"Metrica": "Total Genes Transcriptoma", "Conteo": n_total, "Porcentaje_Transcriptoma": 100.0, "Porcentaje_Proteoma": "-"},
        {"Metrica": "Genes Mapeados a Proteina (Homo sapiens)", "Conteo": n_mapped, "Porcentaje_Transcriptoma": round(pct_mapped, 2), "Porcentaje_Proteoma": 100.0},
        {"Metrica": "Proteinas que Superan Filtro (L <= 1022 AA)", "Conteo": n_passed, "Porcentaje_Transcriptoma": round(pct_transcriptome_passed, 2), "Porcentaje_Proteoma": round(pct_proteome_passed, 2)}
    ])
    recovery_df.to_csv(results_dir / "recovery_table.csv", index=False)

    breakdown_df = pd.DataFrame([
        {"Categoria": "Genes No Codificantes / No Mapeados", "Conteo": n_non_coding, "Porcentaje": round(pct_non_coding, 2)},
        {"Categoria": "Secuencias Excluidas por Longitud (> 1022 AA)", "Conteo": n_exceeded, "Porcentaje": round(pct_exceeded, 2)},
        {"Categoria": "Aprobados en FASTA Curado (L <= 1022 AA)", "Conteo": n_passed, "Porcentaje": round(pct_transcriptome_passed, 2)},
        {"Categoria": "Total Transcriptoma", "Conteo": n_total, "Porcentaje": 100.0}
    ])
    breakdown_df.to_csv(results_dir / "filtration_breakdown.csv", index=False)
    print(f"[Phase D] Tablas generadas: recovery_table.csv y filtration_breakdown.csv")

    # 7. Diagrama de Cascada (Waterfall Chart)
    labels = [
        "Transcriptoma\nInicial",
        "No Codificantes\n/ Sin Mapeo",
        "Proteoma\nMapeado",
        "Longitud\n> 1022 AA",
        "FASTA Curado\nFinal"
    ]
    deltas = [n_total, -n_non_coding, 0, -n_exceeded, 0]
    totals = [n_total, n_total - n_non_coding, n_mapped, n_mapped - n_exceeded, n_passed]
    bottoms = [0, totals[1], 0, totals[3], 0]
    heights = [n_total, n_non_coding, n_mapped, n_exceeded, n_passed]
    colors = ["#2b5c8f", "#d9534f", "#6c757d", "#f0ad4e", "#28a745"]

    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    bars = ax.bar(labels, heights, bottom=bottoms, color=colors, width=0.55, edgecolor="#222222", linewidth=1.2)

    for i, bar in enumerate(bars):
        if i in [0, 2, 4]:
            val = totals[i]
            pct = (val / n_total) * 100
            txt = f"{val:,}\n({pct:.1f}%)"
            y_pos = val + 400
        else:
            val = heights[i]
            pct = (val / n_total) * 100
            txt = f"-{val:,}\n(-{pct:.1f}%)"
            y_pos = bottoms[i] + val / 2 - 400
        ax.text(bar.get_x() + bar.get_width() / 2, y_pos, txt, ha="center", va="bottom" if i in [0, 2, 4] else "center",
                fontsize=9.5, fontweight="bold", color="black" if i in [0, 2, 4] else "white")

    ax.set_ylim(0, 21500)
    ax.set_ylabel("Numero de Genes / Secuencias", fontsize=11, fontweight="bold")
    ax.set_title("exp_02: Cascada de Filtrado Determinístico del Transcriptoma a FASTA Curado (L <= 1022 AA)",
                 fontsize=12, fontweight="bold", pad=15)
    ax.grid(axis="y", linestyle="--", alpha=0.5)
    ax.set_axisbelow(True)

    fig.tight_layout()
    chart_path = reports_fig_dir / "filtering_waterfall.png"
    fig.savefig(chart_path)
    plt.close(fig)
    print(f"[Phase E] Visualización guardada en {chart_path.name}")
    print(f"[Summary] Completado con éxito: {n_passed:,}/{n_total:,} genes retenidos ({pct_transcriptome_passed:.2f}%)")

if __name__ == "__main__":
    main()
