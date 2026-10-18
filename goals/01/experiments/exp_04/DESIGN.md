# Experiment Design: Filtrado DepMap, Proyección t-SNE y Scree Plot PCA
**Experiment**: `goals/01/experiments/exp_04/`  
**Project**: `WSI_RNA_multimodal_learning` (Goal 01)  
**Date**: 2026-10-17  
**Author**: Co-Investigator (AI under active role) & Principal Investigator  
**Status**: Draft

---

## 1. Hypothesis
El subconjunto de embeddings ESM-2 (650M) filtrado a partir del panel de dependencias de DepMap ($N = 1,944$) en `goals/01/data/esm2_650m_depmap_embeddings.h5` exhibe una organización topológica no aleatoria en el espacio latente $D=1280$, donde las proteínas se agrupan en cúmulos consistentes con las macro-familias funcionales y metabólicas de KEGG (Aproximación A), permitiendo una diferenciación clara entre el metabolismo central, la expresión génica/traducción, y los puntos de control del ciclo celular.

---

## 2. Experimental Setup
- **Dataset de Entrada (Embeddings):** `goals/01/data/esm2_650m_embeddings.h5` ($15,940 \times 1280$, `float32`, generado en `exp_03`).
- **Lista de Oncogenes/Dependencias:** `../../data/tcga-brca/preprocessed/genes/DepMap_oncogenes.csv` ($2,263$ genes únicos en formato `GENE_SYMBOL (ENTREZ_ID)`).
- **Alcance de Genes:** Filtrado estricto a los $N = 1,944$ genes retenidos tras el cribado determinístico, omitiendo genes excluidos en fases previas.
- **Anotación Funcional y Ontológica (KEGG - Opción A):**
  - Mapeo determinístico de cada uno de los $1,944$ genes a través de la jerarquía oficial de KEGG BRITE (`br08901`) y la API REST de KEGG (`rest.kegg.jp`):
    1. *Metabolismo & Bioenergética* (Metabolism)
    2. *Expresión Génica & Traducción* (Translation, Transcription, Protein processing)
    3. *Ciclo Celular & Apoptosis* (Cell growth and death)
    4. *Mantenimiento Genómico & Reparación ADN* (Replication and repair)
    5. *Transducción de Señales & Quinasas* (Signal transduction)
    6. *Adhesión Celular & Citoesqueleto* (Cellular community, Motility)
    7. *Otras Rutas Biológicas* (Otras vías ontológicas anotadas)
    8. *Sin anotación KEGG específica* (Fondo neutro gris)
- **Análisis Espectral y Reducción de Dimensionalidad:**
  - **PCA (Scree Plot):** Centrado de matriz $\mathbf{X} \in \mathbb{R}^{1944 \times 1280}$, autovalores y cálculo de dimensión efectiva ($\text{ED}$).
  - **t-SNE 2D:** Reducción con métrica coseno ($\text{perplexity}=30$, $\text{seed}=42$), coloreada por las macro-familias de KEGG con anotación destacada de los oncogenes y dependencias driver presentes (ej. *MYC*, *CCND1*, *MDM2*, *CDK4*).
- **Hardware y Entorno:** Servidor local, entorno conda `histo-DL`.

---

## 3. File Layout for This Experiment
```
goals/01/experiments/exp_04/
├── DESIGN.md                                 ← Este documento (diseño experimental actualizado)
├── IMPLEMENTATION.md                         ← Plan de implementación
├── scripts/
│   └── filter_and_analyze_depmap.py          ← Script autocontenido de filtrado, KEGG, PCA y t-SNE
├── results/
│   ├── depmap_filtration_summary.csv         ← Métricas de recuperación e inventario de exclusiones
│   ├── kegg_pathway_assignments.csv          ← Mapeo formal de los 1,944 genes a macro-familias KEGG
│   ├── pca_variance_stats.csv                ← Varianza explicada y acumulada de componentes
│   ├── pca_scree_plot.png                    ← Scree Plot con varianza individual y acumulada
│   ├── tsne_depmap_2d.png                    ← Proyección topológica t-SNE 2D coloreada por KEGG
│   └── esm2_650m_depmap_embeddings.sha256    ← Hash SHA-256 del artefacto HDF5
└── reports/
    └── summary.md                            ← Reporte formal de consolidación
```

**Artefacto persistente a nivel de Goal:**
- `goals/01/data/esm2_650m_depmap_embeddings.h5` (datasets: `gene_id` $[1944]$, `embeddings` $[1944, 1280]$ en `float32`).

---

## 4. Baselines & Referencias
| Referencia | Origen | Dimensionalidad | Propósito |
| :--- | :--- | :--- | :--- |
| **Embeddings Completos (`exp_03`)** | `goals/01/data/esm2_650m_embeddings.h5` | $15,940 \times 1280$ | Espacio global no filtrado (baseline de dispersión t-SNE) |
| **Panel DepMap Crudo** | `data/tcga-brca/preprocessed/genes/DepMap_oncogenes.csv` | 2,263 genes | Universo de consulta para medir recall y pérdida por longitud |

---

## 5. Proposed Conditions
- **Condición Única (Criterio Determinístico DepMap):**
  - Intersección de símbolos únicos de `DepMap_oncogenes.csv` con los identificadores presentes en `esm2_650m_embeddings.h5`.
  - Normalización estándar (centrado por media) previa a PCA.
  - Distancia coseno para cálculo de afinidad t-SNE.

---

## 6. Evaluation Protocol
- **Métricas Primarias:**
  1. **Tasa de Recuperación DepMap:** $\frac{N_{\text{recuperados}}}{N_{\text{DepMap\_unicos}}} \ge 85.0\%$.
  2. **Varianza Explicada por PC1:** $\le 25.0\%$ (ausencia de dominancia monodimensional colapsada).
  3. **Dimensión Efectiva ($\text{ED}$):** $\ge 15.0$ (riqueza anisotrópica multivariada).
- **Métricas Secundarias:**
  - Varianza explicada acumulada por los primeros 10 y 50 componentes ($PC_{1-10}, PC_{1-50}$).
  - Cohesión topológica observada visualmente en t-SNE 2D.
- **Presupuesto Computacional:**
  - Tiempo estimado de cómputo: $< 30$ segundos (in-memory con `scikit-learn` y `h5py`).

---

## 7. Expected Results & Decision Rules
- **Si la hipótesis se cumple ($N \approx 1,944$, $PC_1 < 25\%$, $\text{ED} \ge 15$):**
  - Se valida el artefacto `goals/01/data/esm2_650m_depmap_embeddings.h5` como insumo primario para las fases de correlación multimodal (Goal 02 / arquitectura WSI-RNA).
  - El Scree Plot confirmará que los embeddings retienen información diversa a lo largo de múltiples dimensiones y no se rigen por un gradiente trivial.
- **Si la hipótesis falla (ej. $PC_1 > 50\%$ o colapso dimensional):**
  - Investigar si el mean-pooling de ESM-2 sobreexpone un sesgo composicional (ej. longitud de secuencia o abundancia de aminoácidos hidrofóbicos) y formular una normalización por capas o estandarización por z-score en un experimento de ablación.

---

## 8. Risks & Mitigations
- **Riesgo 1: Desalineación de Símbolos Genéticos (Sinónimos/Alias):**
  - *Mitigación:* Parseo robusto mediante expresiones regulares para aislar el símbolo primario previo a los paréntesis con el Entrez ID. Documentar los 319 genes no recuperados en `results/depmap_filtration_summary.csv` cruzándolos contra las exclusiones de `exp_02` ($L > 1022$ AA).
- **Riesgo 2: Sensibilidad de t-SNE a Semilla o Perplejidad:**
  - *Mitigación:* Fijar estrictamente `random_state=42`, utilizar `metric='cosine'` para respetar la geometría de embeddings de transformers biológicos, y contrastar con el Scree Plot analítico de PCA que es estrictamente invariante a hiperparámetros estocásticos.

---

## 9. Reproducibility Checklist
- [ ] Semillas fijadas (`random_state=42`) para t-SNE y PCA.
- [ ] Versión del entorno conda documentada (`histo-DL`).
- [ ] Artefacto de embeddings filtrados persistido en `goals/01/data/esm2_650m_depmap_embeddings.h5`.
- [ ] Tablas de auditoría persistidas en `results/`.
- [ ] Commit / Hash de versión registrado en `results/version.txt`.

---

## 10. Next Steps
1. Aprobar este plan de diseño experimental (`DESIGN.md`).
2. Elaborar y aprobar el plan de implementación (`IMPLEMENTATION.md`) especificando la estructura de `filter_and_analyze_depmap.py`.
3. Ejecutar el script y generar los artefactos en `goals/01/data/` y `results/`.
