#!/usr/bin/env python3
"""
exp_03: Reducción Dimensional t-SNE 2D y Visualización de Embeddings ESM-2 (650M).
Lee goals/01/data/esm2_650m_embeddings.h5 y genera reports/figures/tsne_embeddings.png.
"""

import time
from pathlib import Path
import numpy as np
import h5py
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from scipy.stats import gaussian_kde


def main():
    t0 = time.time()
    script_dir = Path(__file__).resolve().parent
    exp_dir = script_dir.parent
    goal_dir = exp_dir.parent.parent

    h5_path = goal_dir / "data/esm2_650m_embeddings.h5"
    results_dir = exp_dir / "results"
    figures_dir = exp_dir / "reports/figures"
    results_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    print(f"[{time.strftime('%H:%M:%S')}] Iniciando análisis t-SNE para exp_03")
    print(f"[{time.strftime('%H:%M:%S')}] Leyendo embeddings desde {h5_path}...")

    with h5py.File(h5_path, "r") as f:
        embeddings = f["embeddings"][:]
        raw_genes = f["gene_id"][:]
        # Decodificar strings si vienen en bytes
        if isinstance(raw_genes[0], bytes):
            genes = [g.decode("utf-8") for g in raw_genes]
        else:
            genes = [str(g) for g in raw_genes]

    n_samples, n_features = embeddings.shape
    print(f"[{time.strftime('%H:%M:%S')}] Matriz cargada: {n_samples:,} secuencias x {n_features} dimensiones")

    # 1. Pre-reducción con PCA a 50 componentes para aceleración y preservación global
    print(f"[{time.strftime('%H:%M:%S')}] Computando PCA (50 componentes)...")
    pca = PCA(n_components=50, random_state=42)
    x_pca = pca.fit_transform(embeddings)
    var_exp = float(np.sum(pca.explained_variance_ratio_) * 100.0)
    print(f"[{time.strftime('%H:%M:%S')}] Varianza explicada por 50 PCs: {var_exp:.2f}%")

    # 2. Reducción t-SNE 2D
    print(f"[{time.strftime('%H:%M:%S')}] Computando t-SNE (perplexity=30, random_state=42)...")
    tsne = TSNE(
        n_components=2,
        perplexity=30,
        random_state=42,
        init="pca",
        learning_rate="auto",
        max_iter=1000,
        n_jobs=-1
    )
    coords_2d = tsne.fit_transform(x_pca)
    print(f"[{time.strftime('%H:%M:%S')}] t-SNE completado en {time.time() - t0:.1f} s")

    # 3. Guardar coordenadas tabuladas
    coords_df = pd.DataFrame({
        "gene_id": genes,
        "tsne_dim1": coords_2d[:, 0],
        "tsne_dim2": coords_2d[:, 1]
    })
    coords_csv = results_dir / "tsne_coordinates.csv"
    coords_df.to_csv(coords_csv, index=False)
    print(f"[{time.strftime('%H:%M:%S')}] Coordenadas guardadas en {coords_csv}")

    # 4. Estimar densidad para coloración de alta fidelidad
    print(f"[{time.strftime('%H:%M:%S')}] Calculando densidad kernel para la visualización...")
    # Muestreo representativo para cálculo rápido de densidad KDE si N es grande
    sample_idx = np.random.RandomState(42).choice(n_samples, size=min(4000, n_samples), replace=False)
    kde = gaussian_kde(coords_2d[sample_idx].T)
    density = kde(coords_2d.T)

    # 5. Generar Figura t-SNE 2D
    print(f"[{time.strftime('%H:%M:%S')}] Renderizando figura t-SNE...")
    fig, ax = plt.subplots(figsize=(10, 8), dpi=300)

    # Ordenar puntos por densidad para que los picos queden al frente
    sort_idx = np.argsort(density)
    x = coords_2d[sort_idx, 0]
    y = coords_2d[sort_idx, 1]
    c = density[sort_idx]

    sc = ax.scatter(
        x, y,
        c=c,
        cmap="viridis",
        s=4.0,
        alpha=0.6,
        edgecolors="none",
        rasterized=True
    )

    cbar = plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label("Densidad Relativa en Espacio Latente", fontsize=11, labelpad=8)
    cbar.ax.tick_params(labelsize=9)

    ax.set_title(
        "Distribución Latente de Embeddings Proteicos (ESM-2 650M)\n"
        "Proyección t-SNE 2D · Residue Mean Pooling ($D=1280$)",
        fontsize=13,
        fontweight="bold",
        pad=14
    )
    ax.set_xlabel("t-SNE Dimensión 1", fontsize=11)
    ax.set_ylabel("t-SNE Dimensión 2", fontsize=11)
    ax.grid(True, linestyle="--", alpha=0.3)

    # Cuadro informativo con metadatos del experimento
    info_text = (
        f"Proteoma Curado: $N = {n_samples:,}$\n"
        f"Modelo: ESM-2 650M (Capa 33)\n"
        f"Dimensión original: {n_features}\n"
        f"Varianza PCA (50 PCs): {var_exp:.1f}%\n"
        f"Perplejidad: 30 · Semilla: 42"
    )
    ax.text(
        0.03, 0.03, info_text,
        transform=ax.transAxes,
        fontsize=9.5,
        verticalalignment="bottom",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="white", alpha=0.85, edgecolor="#cccccc")
    )

    plt.tight_layout()
    output_png = figures_dir / "tsne_embeddings.png"
    plt.savefig(output_png, dpi=300)
    plt.close()
    print(f"[{time.strftime('%H:%M:%S')}] Figura guardada en {output_png}")
    print(f"[{time.strftime('%H:%M:%S')}] Proceso de análisis t-SNE concluido.")


if __name__ == "__main__":
    main()
