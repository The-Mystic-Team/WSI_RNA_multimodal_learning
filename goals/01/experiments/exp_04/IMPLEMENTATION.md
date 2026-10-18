# Implementation Plan: exp_04 — Filtrado DepMap, Proyección t-SNE con Macro-Familias KEGG y Scree Plot PCA

**Experiment ID**: `exp_04`  
**Goal**: `01`  
**Companion Design**: [`DESIGN.md`](./DESIGN.md)  
**Date**: 2026-10-17  
**Status**: In Progress  

---

## 1. Overview & Operational Scope

Implementar y ejecutar el script `scripts/filter_and_analyze_depmap.py` sobre los **1,944 genes retenidos** para:
1. **Curaduría y Persistencia de Embeddings ($N = 1,944$):**
   - Confirmar el artefacto `goals/01/data/esm2_650m_depmap_embeddings.h5` con los $1,944$ genes filtrados y su hash criptográfico SHA-256.
2. **Anotación Funcional y Ontológica (KEGG Macro-Familias - Opción A):**
   - Consultar la jerarquía de KEGG BRITE (`br08901`) y el mapeo de genes humanos (`rest.kegg.jp/link/pathway/hsa`).
   - Mapear cada uno de los $1,944$ genes a su macro-familia funcional primaria mediante jerarquía de precedencia:
     - *Metabolismo & Bioenergética*
     - *Expresión Génica & Traducción*
     - *Ciclo Celular & Apoptosis*
     - *Mantenimiento Genómico & Reparación ADN*
     - *Transducción de Señales & Quinasas*
     - *Adhesión Celular & Citoesqueleto*
     - *Otras Rutas Biológicas*
     - *Sin anotación KEGG específica*
   - Guardar la tabla formal de anotaciones en `results/kegg_pathway_assignments.csv`.
3. **Análisis Espectral (PCA & Scree Plot):**
   - Varianza explicada y acumulada de componentes principales y cálculo de la dimensión efectiva ($\text{ED}$).
   - Guardar `results/pca_scree_plot.png` y `results/pca_variance_stats.csv`.
4. **Proyección Topológica t-SNE 2D:**
   - Proyectar los $1,944$ embeddings con métrica coseno (`perplexity=30`, `random_state=42`).
   - Colorear los puntos de acuerdo a las macro-familias KEGG.
   - Anotar de forma destacada los oncogenes/drivers presentes en el panel (ej. *MYC*, *CCND1*, *MDM2*, *CDK4*).
   - Generar la figura `results/tsne_depmap_2d.png`.
5. **Auditoría de Filtración:**
   - Documentar la causa de exclusión de los genes no incluidos en `results/depmap_filtration_summary.csv`.

---

## 2. Directory & Artifact Structure

```
goals/01/
├── data/
│   ├── esm2_650m_embeddings.h5               ← Entrada: Embeddings completos (15,940 genes)
│   └── esm2_650m_depmap_embeddings.h5        ← ARTEFACTO PRIMARIO GENERADO (1,944 genes, D=1280)
└── experiments/exp_04/
    ├── DESIGN.md                             ← Diseño experimental ratificado
    ├── IMPLEMENTATION.md                     ← Este documento
    ├── scripts/
    │   └── filter_and_analyze_depmap.py      ← Script autocontenido
    ├── results/
    │   ├── depmap_filtration_summary.csv     ← Auditoría completa de retención y exclusiones
    │   ├── kegg_pathway_assignments.csv      ← Asignaciones ontológicas de macro-familias KEGG
    │   ├── pca_variance_stats.csv            ← Razones de varianza explicada y acumulada
    │   ├── pca_scree_plot.png                ← Gráfica Scree Plot (dominancia vs. dimensionalidad)
    │   ├── tsne_depmap_2d.png                ← Gráfica t-SNE 2D coloreada por KEGG Macro-Familias
    │   ├── esm2_650m_depmap_embeddings.sha256← Checksum SHA-256 del artefacto HDF5
    │   └── version.txt                       ← Hash git y timestamp de ejecución
    └── reports/
        └── summary.md                        ← Reporte científico consolidado post-ejecución
```
        └── summary.md                        ← Reporte científico consolidado post-ejecución
```

---

## 3. Algorithm & Code Architecture (`scripts/filter_and_analyze_depmap.py`)

El script se diseñará de forma modular, con código legible y sin boilerplate innecesario, estructurado en las siguientes fases:

### Fase 1: Ingesta y Cruce
```python
# Rutas relativas a la raíz del proyecto
DEP_MAP_PATH = "data/tcga-brca/preprocessed/genes/DepMap_oncogenes.csv"
INPUT_H5_PATH = "goals/01/data/esm2_650m_embeddings.h5"
OUTPUT_H5_PATH = "goals/01/data/esm2_650m_depmap_embeddings.h5"
EXP_DIR = "goals/01/experiments/exp_04"
```
1. Cargar `DepMap_oncogenes.csv`, aislar símbolos HGNC limpios y eliminar duplicados.
2. Abrir `goals/01/data/esm2_650m_embeddings.h5`, decodificar `gene_id` a lista de strings.
3. Crear mapeo de índices: `gene_to_idx = {g: i for i, g in enumerate(h5_genes)}`.
4. Identificar `matched_genes` preservando orden alfabético determinístico y extraer los índices correspondientes.

### Fase 2: Escritura HDF5 y Checksum
1. Leer los vectores de la matriz `embeddings` indexados por `matched_indices`.
2. Escribir `goals/01/data/esm2_650m_depmap_embeddings.h5`:
   - Dataset `gene_id`: array 1D de strings UTF-8.
   - Dataset `embeddings`: matriz $[1944, 1280]$ en `float32` con compresión `gzip` ligera (`compression_opts=4`).
3. Calcular hash SHA-256 del archivo y registrar en `results/esm2_650m_depmap_embeddings.sha256`.

### Fase 3: Descomposición Espectral PCA y Scree Plot
1. Centrar los embeddings restando la media por dimensión: $\tilde{\mathbf{X}} = \mathbf{X} - \bar{\mathbf{X}}$.
2. Ajustar `sklearn.decomposition.PCA(n_components=min(100, N))` sobre $\tilde{\mathbf{X}}$.
3. Computar autovalores y ratios de varianza explicada.
4. Calcular Dimensión Efectiva (Participation Ratio):
   $$\text{ED} = \frac{(\sum \lambda_i)^2}{\sum \lambda_i^2}$$
5. Generar figura `results/pca_scree_plot.png` con Matplotlib/Seaborn:
   - Panel superior: Barras de varianza individual por PC ($PC_1 \dots PC_{50}$) con línea de codo.
   - Panel inferior: Curva de varianza acumulada con líneas guía en 50%, 80% y 95%.
   - Anotación clara del valor de $\text{ED}$ y varianza de $PC_1$.

### Fase 4: Proyección t-SNE 2D
1. Instanciar `sklearn.manifold.TSNE`:
   - `n_components=2`, `perplexity=30`, `metric='cosine'`, `random_state=42`, `init='pca'`.
2. Proyectar a $\mathbf{Y} \in \mathbb{R}^{1944 \times 2}$.
3. Generar `results/tsne_depmap_2d.png`:
   - Dispersión con densidad (*hexbin* o *scatter* con alpha).
   - Resaltar y etiquetar oncogenes canónicos de cáncer de mama (*TP53*, *PIK3CA*, *ERBB2*, *ESR1*, *GATA3*, *PTEN*, *CDH1*, *BRCA1*).

### Fase 5: Auditoría y Resumen
1. Cruzar los genes faltantes ($N=319$) con `goals/01/experiments/exp_02/results/excluded_genes_summary.csv`.
2. Guardar `results/depmap_filtration_summary.csv`.
3. Guardar `results/pca_variance_stats.csv`.
4. Registrar `results/version.txt` con el hash git actual y timestamp.

---

## 4. Execution Protocol & Telemetry

- **Entorno:** `histo-DL` (`/home/jmalagont/miniconda3/envs/histo-DL/bin/python`).
- **Tiempo estimado:** ~15 a 25 segundos (ejecución síncrona en primer plano según directrices para cómputos $< 60$ s).
- **Comando de Ejecución:**
  ```bash
  conda activate histo-DL && python goals/01/experiments/exp_04/scripts/filter_and_analyze_depmap.py
  ```
- **Telemetría Condensada (Formato Esperado):**
  Una línea condensada por hito:
  ```
  [exp_04] Ingesta: 2263 genes DepMap parseados | 1944 coincidentes (85.90%) | 319 excluidos
  [exp_04] HDF5 generado: goals/01/data/esm2_650m_depmap_embeddings.h5 [1944, 1280] | SHA256 calculado
  [exp_04] PCA completado: PC1 = X.X% | ED = XX.X | Codo en PC_k
  [exp_04] t-SNE 2D completado en X.Xs | Figuras y tablas persistidas en results/
  ```

---

## 5. Verification & Acceptance Criteria

| Criterio | Umbral / Especificación | Validación |
| :--- | :--- | :--- |
| **Recuperación DepMap** | $\ge 85.0\%$ ($N \ge 1,920$) | `results/depmap_filtration_summary.csv` |
| **Dimensiones HDF5** | Exactamente `[1944, 1280]` | `h5py.File['embeddings'].shape` |
| **Varianza PC1** | $< 25.0\%$ | `results/pca_variance_stats.csv` |
| **Dimensión Efectiva** | $\text{ED} \ge 15.0$ | `results/pca_variance_stats.csv` |
| **Calidad Gráfica** | Figuras legibles $\ge 300$ DPI | `pca_scree_plot.png` y `tsne_depmap_2d.png` |

---

## 6. Next Steps
Tras la aprobación de este plan de implementación por parte del usuario, se procederá a:
1. Crear el script `goals/01/experiments/exp_04/scripts/filter_and_analyze_depmap.py`.
2. Actualizar el estado a `running` en `goals/01/experiments/index.md`.
3. Ejecutar el script bajo el entorno `histo-DL`.
4. Evaluar los resultados generados y consolidar el reporte en `goals/01/experiments/exp_04/reports/summary.md` y `.logbook.md`.
