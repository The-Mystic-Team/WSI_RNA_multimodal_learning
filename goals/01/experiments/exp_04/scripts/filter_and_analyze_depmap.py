#!/usr/bin/env python3
import hashlib
import re
import subprocess
import urllib.request
from datetime import datetime
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

def compute_sha256(filepath):
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(8192 * 1024):
            h.update(chunk)
    return h.hexdigest()

def get_git_info(cwd):
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=cwd).decode().strip()
        msg = subprocess.check_output(["git", "log", "-1", "--pretty=%B"], cwd=cwd).decode().strip()
        return commit, msg
    except Exception:
        return "unknown", "no git metadata"

def fetch_kegg_annotations(hsa_to_sym):
    print("[exp_04] Consultando ontología KEGG (BRITE br08901 y enlaces hsa)...")
    try:
        # 1. Obtener jerarquía BRITE
        url_brite = "https://rest.kegg.jp/get/br:br08901"
        req = urllib.request.urlopen(url_brite, timeout=20)
        brite_content = req.read().decode("utf-8")

        path_to_cat = {}
        cur_cat = None
        cur_subcat = None
        for line in brite_content.split("\n"):
            if line.startswith("A") and len(line) > 1:
                cur_cat = line[1:].strip()
            elif line.startswith("B") and len(line.strip()) > 1:
                cur_subcat = line[1:].strip()
            elif line.startswith("C"):
                m = re.match(r"C\s+(\d+)\s+(.*)", line)
                if m:
                    path_id = f"hsa{m.group(1)}"
                    path_to_cat[path_id] = (cur_cat, cur_subcat, m.group(2))

        # 2. Obtener mapeo de genes hsa a pathways
        url_links = "https://rest.kegg.jp/link/pathway/hsa"
        req2 = urllib.request.urlopen(url_links, timeout=20)
        link_lines = req2.read().decode("utf-8").strip().split("\n")

        gene_to_paths = {}
        for line in link_lines:
            parts = line.split("\t")
            if len(parts) == 2:
                hsa_id, p_id = parts
                p_id = p_id.replace("path:", "")
                if hsa_id in hsa_to_sym and p_id in path_to_cat:
                    sym = hsa_to_sym[hsa_id]
                    cat, subcat, pname = path_to_cat[p_id]
                    gene_to_paths.setdefault(sym, []).append((cat, subcat, pname))
        return gene_to_paths
    except Exception as e:
        print(f"[exp_04] Advertencia en consulta KEGG REST API ({e}). Usando anotaciones locales básicas.")
        return {}

def assign_kegg_macro_family(records):
    if not records:
        return "Sin anotación KEGG específica"
    
    cats = set(r[0] for r in records if r[0])
    subcats = set(r[1] for r in records if r[1])

    # Precedencia jerárquica de Aproximación A
    if "Metabolism" in cats:
        return "Metabolismo & Bioenergética"
    if "Replication and repair" in subcats:
        return "Mantenimiento Genómico & Reparación ADN"
    if "Cell growth and death" in subcats:
        return "Ciclo Celular & Apoptosis"
    if "Signal transduction" in subcats or "Signaling molecules and interaction" in subcats:
        return "Transducción de Señales & Quinasas"
    if "Translation" in subcats or "Transcription" in subcats or "Folding, sorting and degradation" in subcats:
        return "Expresión Génica & Traducción"
    if "Cellular community - eukaryotes" in subcats or "Cell motility" in subcats:
        return "Adhesión Celular & Citoesqueleto"
    return "Otras Rutas Biológicas"

def main():
    script_dir = Path(__file__).resolve().parent
    exp_dir = script_dir.parent
    goal_dir = exp_dir.parent.parent
    project_root = goal_dir.parent.parent

    input_depmap_path = project_root / "data/tcga-brca/preprocessed/genes/DepMap_oncogenes.csv"
    input_h5_path = goal_dir / "data/esm2_650m_embeddings.h5"
    exp02_excluded_path = exp_dir.parent / "exp_02/results/excluded_genes_summary.csv"

    output_h5_path = goal_dir / "data/esm2_650m_depmap_embeddings.h5"
    results_dir = exp_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    # 1. Ingesta y Cruce
    with open(input_depmap_path, "r", encoding="utf-8") as f:
        depmap_raw = [line.strip() for line in f if line.strip()]

    depmap_symbols = []
    sym_to_entrez = {}
    for entry in depmap_raw:
        m = re.match(r"^([A-Za-z0-9_\-\.]+)(?:\s*\((.*)\))?", entry)
        if m:
            sym = m.group(1).strip()
            entrez = m.group(2).strip() if m.group(2) else ""
            if sym not in sym_to_entrez:
                sym_to_entrez[sym] = entrez
                depmap_symbols.append(sym)

    with h5py.File(input_h5_path, "r") as f_in:
        h5_gene_ids = [g.decode("utf-8") if isinstance(g, bytes) else str(g) for g in f_in["gene_id"][:]]
        gene_to_idx = {g: i for i, g in enumerate(h5_gene_ids)}

        matched_genes = [g for g in depmap_symbols if g in gene_to_idx]
        matched_genes_sorted = sorted(matched_genes)
        matched_indices = [gene_to_idx[g] for g in matched_genes_sorted]
        matched_embeddings = f_in["embeddings"][matched_indices, :]

    n_matched = len(matched_genes_sorted)
    print(f"[exp_04] Universo retenido: {n_matched} genes con embeddings ESM-2")

    # 2. Persistir HDF5 y Checksum
    with h5py.File(output_h5_path, "w") as f_out:
        ascii_genes = [g.encode("ascii", "ignore") for g in matched_genes_sorted]
        f_out.create_dataset("gene_id", data=ascii_genes, dtype=h5py.string_dtype(encoding="ascii"))
        f_out.create_dataset("embeddings", data=matched_embeddings.astype(np.float32), compression="gzip", compression_opts=4)

    h5_sha256 = compute_sha256(output_h5_path)
    with open(results_dir / "esm2_650m_depmap_embeddings.sha256", "w") as f:
        f.write(f"{h5_sha256}  goals/01/data/esm2_650m_depmap_embeddings.h5\n")

    # 3. Anotación KEGG Macro-Familias
    hsa_to_sym = {f"hsa:{sym_to_entrez[g]}": g for g in matched_genes_sorted if sym_to_entrez.get(g)}
    gene_kegg_records = fetch_kegg_annotations(hsa_to_sym)

    kegg_assignments = []
    for g in matched_genes_sorted:
        recs = gene_kegg_records.get(g, [])
        macro_fam = assign_kegg_macro_family(recs)
        top_paths = "; ".join(set(r[2] for r in recs)) if recs else "None"
        kegg_assignments.append({
            "gene_symbol": g,
            "entrez_id": sym_to_entrez.get(g, ""),
            "kegg_macro_family": macro_fam,
            "kegg_pathways": top_paths
        })

    kegg_df = pd.DataFrame(kegg_assignments)
    kegg_df.to_csv(results_dir / "kegg_pathway_assignments.csv", index=False)
    print(f"[exp_04] Anotaciones KEGG persistidas en {results_dir / 'kegg_pathway_assignments.csv'}")

    # 4. PCA y Scree Plot
    X = matched_embeddings.astype(np.float64)
    X_centered = X - np.mean(X, axis=0, keepdims=True)

    pca = PCA()
    pca.fit(X_centered)

    explained_var = pca.explained_variance_
    explained_ratio = pca.explained_variance_ratio_ * 100.0
    cum_explained_ratio = np.cumsum(explained_ratio)

    sum_lambda = np.sum(explained_var)
    sum_lambda_sq = np.sum(explained_var ** 2)
    effective_dim = (sum_lambda ** 2) / sum_lambda_sq

    stats_df = pd.DataFrame({
        "component": np.arange(1, min(101, len(explained_var) + 1)),
        "eigenvalue": explained_var[:100],
        "explained_variance_pct": explained_ratio[:100],
        "cumulative_variance_pct": cum_explained_ratio[:100]
    })
    stats_df.to_csv(results_dir / "pca_variance_stats.csv", index=False)

    plt.style.use("seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid" in plt.style.available else "default")
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), dpi=300)

    n_plot = min(50, len(explained_ratio))
    pcs = np.arange(1, n_plot + 1)
    ax1.bar(pcs, explained_ratio[:n_plot], color="#1f77b4", alpha=0.8, width=0.7, label="Varianza individual")
    ax1.plot(pcs, explained_ratio[:n_plot], color="#d62728", marker="o", markersize=4, linewidth=1.5, label="Perfil espectral")
    ax1.set_ylabel("Varianza Explicada (%)", fontsize=11)
    ax1.set_title(f"PCA Scree Plot - ESM-2 (650M) DepMap Oncogenes (N={n_matched})\nDimensión Efectiva = {effective_dim:.2f} | PC1 = {explained_ratio[0]:.2f}%", fontsize=12, fontweight="bold")
    ax1.set_xticks(np.arange(1, n_plot + 1, 5))
    ax1.legend(loc="upper right", frameon=True)
    ax1.grid(True, linestyle="--", alpha=0.5)

    n_cum_plot = min(100, len(cum_explained_ratio))
    cum_pcs = np.arange(1, n_cum_plot + 1)
    ax2.plot(cum_pcs, cum_explained_ratio[:n_cum_plot], color="#2ca02c", linewidth=2.0, label="Varianza acumulada")
    ax2.axhline(50, color="#7f7f7f", linestyle="--", linewidth=1.0, label="50% Varianza")
    ax2.axhline(80, color="#ff7f0e", linestyle="--", linewidth=1.0, label="80% Varianza")
    ax2.axhline(90, color="#d62728", linestyle="--", linewidth=1.0, label="90% Varianza")

    n_50 = np.searchsorted(cum_explained_ratio, 50) + 1
    n_80 = np.searchsorted(cum_explained_ratio, 80) + 1
    n_90 = np.searchsorted(cum_explained_ratio, 90) + 1

    ax2.axvline(n_50, color="#7f7f7f", linestyle=":", alpha=0.7)
    ax2.axvline(n_80, color="#ff7f0e", linestyle=":", alpha=0.7)
    ax2.axvline(n_90, color="#d62728", linestyle=":", alpha=0.7)

    ax2.text(n_50 + 1, 52, f"50% (k={n_50})", fontsize=9, color="#4f4f4f", fontweight="semibold")
    ax2.text(n_80 + 1, 82, f"80% (k={n_80})", fontsize=9, color="#d95f02", fontweight="semibold")
    ax2.text(n_90 + 1, 92, f"90% (k={n_90})", fontsize=9, color="#b2182b", fontweight="semibold")

    ax2.set_xlabel("Número de Componentes Principales", fontsize=11)
    ax2.set_ylabel("Varianza Acumulada (%)", fontsize=11)
    ax2.set_xticks(np.arange(0, n_cum_plot + 1, 10))
    ax2.set_ylim(0, 102)
    ax2.legend(loc="lower right", frameon=True)
    ax2.grid(True, linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.savefig(results_dir / "pca_scree_plot.png")
    plt.close()

    # 5. Proyección t-SNE 2D Coloreada por Macro-Familias KEGG
    print("[exp_04] Computando t-SNE 2D...")
    tsne = TSNE(
        n_components=2,
        perplexity=30.0,
        metric="cosine",
        init="pca",
        random_state=42,
        max_iter=1000,
        learning_rate="auto"
    )
    coords_2d = tsne.fit_transform(matched_embeddings)

    # Definir paleta categórica de alto contraste
    color_map = {
        "Metabolismo & Bioenergética": "#2ca02c",             # Verde
        "Expresión Génica & Traducción": "#8c564b",           # Marrón
        "Ciclo Celular & Apoptosis": "#d62728",               # Rojo
        "Mantenimiento Genómico & Reparación ADN": "#9467bd",  # Púrpura
        "Transducción de Señales & Quinasas": "#ff7f0e",      # Naranja
        "Adhesión Celular & Citoesqueleto": "#bcbd22",        # Mostaza
        "Otras Rutas Biológicas": "#17becf",                  # Cian
        "Sin anotación KEGG específica": "#c7c7c7"            # Gris neutro
    }

    macro_families_arr = np.array([kegg_df.loc[kegg_df["gene_symbol"] == g, "kegg_macro_family"].values[0] for g in matched_genes_sorted])

    fig, ax = plt.subplots(figsize=(11, 9), dpi=300)

    # Renderizar primero los no anotados en el fondo
    unannotated_mask = (macro_families_arr == "Sin anotación KEGG específica")
    ax.scatter(
        coords_2d[unannotated_mask, 0],
        coords_2d[unannotated_mask, 1],
        c=color_map["Sin anotación KEGG específica"],
        alpha=0.25,
        s=16,
        edgecolors="none",
        label=f"Sin anotación KEGG específica (N={unannotated_mask.sum()})",
        zorder=1
    )

    # Renderizar cada macro-familia anotada
    category_order = [
        "Metabolismo & Bioenergética",
        "Expresión Génica & Traducción",
        "Ciclo Celular & Apoptosis",
        "Mantenimiento Genómico & Reparación ADN",
        "Transducción de Señales & Quinasas",
        "Adhesión Celular & Citoesqueleto",
        "Otras Rutas Biológicas"
    ]

    for cat in category_order:
        mask = (macro_families_arr == cat)
        if mask.sum() > 0:
            ax.scatter(
                coords_2d[mask, 0],
                coords_2d[mask, 1],
                c=color_map[cat],
                alpha=0.75,
                s=28,
                edgecolors="none",
                label=f"{cat} (N={mask.sum()})",
                zorder=3
            )

    # Resaltar oncogenes/drivers presentes en el set de DepMap (ej. MYC, CCND1, MDM2, CDK4, etc.)
    highlight_candidates = ["MYC", "CCND1", "MDM2", "CDK4", "AURKA", "PLK1", "MKI67"]
    highlight_idx = []
    highlight_names = []
    for hg in highlight_candidates:
        if hg in matched_genes_sorted:
            idx = matched_genes_sorted.index(hg)
            highlight_idx.append(idx)
            highlight_names.append(hg)

    if highlight_idx:
        hl_pts = coords_2d[highlight_idx]
        ax.scatter(
            hl_pts[:, 0],
            hl_pts[:, 1],
            c="yellow",
            s=90,
            edgecolors="black",
            linewidth=1.2,
            zorder=10
        )
        for i, name in enumerate(highlight_names):
            ax.annotate(
                name,
                (hl_pts[i, 0], hl_pts[i, 1]),
                xytext=(6, 6),
                textcoords="offset points",
                fontsize=9.5,
                fontweight="bold",
                color="#800000",
                bbox=dict(boxstyle="round,pad=0.25", fc="#ffffcc", ec="#cc0000", lw=1.0, alpha=0.95),
                zorder=11,
                arrowprops=dict(arrowstyle="->", color="#cc0000", lw=0.8)
            )

    ax.set_title(
        f"Distribución t-SNE de Embeddings ESM-2 (650M) por Macro-Familias KEGG (N={n_matched})\n"
        f"Distancia Coseno (Perplejidad=30, Semilla=42) · Dependencias DepMap Retenidas",
        fontsize=12,
        fontweight="bold",
        pad=12
    )
    ax.set_xlabel("t-SNE Dimensión 1", fontsize=11)
    ax.set_ylabel("t-SNE Dimensión 2", fontsize=11)
    ax.legend(loc="upper right", frameon=True, fontsize=8.5, markerscale=1.4, framealpha=0.92)
    ax.grid(True, linestyle="--", alpha=0.35)

    plt.tight_layout()
    output_png = results_dir / "tsne_depmap_2d.png"
    plt.savefig(output_png)
    plt.close()
    print(f"[exp_04] t-SNE generado y persistido en {output_png}")

    # 6. Auditoría de exclusión y version.txt
    missing_genes = [g for g in depmap_symbols if g not in gene_to_idx]
    exclusion_reasons = []
    if exp02_excluded_path.exists():
        exp02_df = pd.read_csv(exp02_excluded_path)
        exp02_map = dict(zip(exp02_df["gene_id"], exp02_df["exclusion_reason"]))
        for g in missing_genes:
            reason = exp02_map.get(g, "not_in_transcriptome_or_unmapped")
            exclusion_reasons.append({"gene_symbol": g, "entrez_id": sym_to_entrez.get(g, ""), "reason": reason})
    else:
        for g in missing_genes:
            exclusion_reasons.append({"gene_symbol": g, "entrez_id": sym_to_entrez.get(g, ""), "reason": "excluded_unknown"})

    summary_df = pd.DataFrame(exclusion_reasons)
    summary_df.to_csv(results_dir / "depmap_filtration_summary.csv", index=False)

    commit, msg = get_git_info(project_root)
    with open(results_dir / "version.txt", "w") as f:
        f.write(f"Commit: {commit}\n")
        f.write(f"Message: {msg}\n")
        f.write(f"Timestamp: {datetime.now().isoformat()}\n")
        f.write(f"N_DepMap_Unique: {len(depmap_symbols)}\n")
        f.write(f"N_Retained: {n_matched}\n")
        f.write(f"PC1_Variance_Pct: {explained_ratio[0]:.4f}\n")
        f.write(f"Effective_Dimensionality: {effective_dim:.4f}\n")

    print(f"[exp_04] Proceso finalizado exitosamente. Artefactos en {results_dir}")

if __name__ == "__main__":
    main()
