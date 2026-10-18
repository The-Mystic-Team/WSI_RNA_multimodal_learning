# Experiment Summary: exp_04 — Filtrado DepMap, Proyección t-SNE con Macro-Familias KEGG y Scree Plot PCA

**Experiment ID**: `exp_04`  
**Goal**: `01`  
**Companion Design**: [`DESIGN.md`](../DESIGN.md)  
**Companion Implementation**: [`IMPLEMENTATION.md`](../IMPLEMENTATION.md)  
**Date**: 2026-10-17  
**Status**: Complete  
**Verdict**: Success (Hypothesis Confirmed — Organización Topológica por Macro-Familias KEGG y Diagnóstico de Dominancia Espectral)  

---

## 1. Executive Summary

El experimento `exp_04` filtró el espacio completo de representaciones latentes de proteínas de `exp_03` ([`goals/01/data/esm2_650m_embeddings.h5`](file:///home/jmalagont/project/WSI_RNA_multimodal_learning/goals/01/data/esm2_650m_embeddings.h5), $N=15,940$) utilizando el catálogo de dependencias esenciales de DepMap ([`data/tcga-brca/preprocessed/genes/DepMap_oncogenes.csv`](file:///home/jmalagont/project/WSI_RNA_multimodal_learning/data/tcga-brca/preprocessed/genes/DepMap_oncogenes.csv)), reteniendo con éxito **$1,944$ genes** ($85.90\%$ de recuperación del panel).

En respuesta a la directriz de la investigación, se implementó la **Aproximación A (Macro-Familias Biológicas)** empleando la ontología oficial de **KEGG** (*Kyoto Encyclopedia of Genes and Genomes*, jerarquía BRITE `br08901`) para anotar y colorear la distribución de los $1,944$ genes en el espacio latente t-SNE 2D, omitiendo deliberadamente los genes excluidos en etapas previas para preservar la integridad del panel cribado.

### Conclusiones Principales:
1. **Anotación Funcional y Metabólica KEGG (Opción A):**
   - **Metabolismo & Bioenergética:** **225 genes** (Glicólisis, Ciclo de Krebs, Fosforilación Oxidativa, Metabolismo de lípidos, nucleótidos y aminoácidos).
   - **Expresión Génica & Traducción:** **510 genes** (Subunidades ribosomales, biogénesis de ribosomas, espliceosoma y factores de iniciación).
   - **Ciclo Celular & Apoptosis:** **109 genes** (Puntos de control mitóticos, reguladores de ciclo, caspasas).
   - **Mantenimiento Genómico & Reparación ADN:** **85 genes** (Recombinación homóloga, reparación por escisión, helicasas).
   - **Transducción de Señales & Quinasas:** **73 genes** (Cascadas de señalización MAPK, mTOR, quinasas reguladoras).
   - **Adhesión Celular & Citoesqueleto:** **20 genes** (Uniones celulares, adhesión focal, dinámica de actina).
   - **Otras Rutas Biológicas:** **145 genes** (Transportadores y respuestas a estímulos).
   - **Sin anotación KEGG específica:** **777 genes** (Representados en gris neutro en el fondo).
2. **Proyección t-SNE 2D Estructurada:**
   - La figura [`results/tsne_depmap_2d.png`](file:///home/jmalagont/project/WSI_RNA_multimodal_learning/goals/01/experiments/exp_04/results/tsne_depmap_2d.png) evidencia una estructuración topológica clara: las enzimas metabólicas y los complejos ribosomales/espliceosómicos forman cúmulos densos diferenciados, mientras que los reguladores del ciclo celular se ubican en ramas especializadas.
   - Se destacan y anotan con callouts específicos los oncogenes y dependencias driver presentes en el panel (ej. ***MYC***, ***CCND1***, ***MDM2***, ***CDK4***, ***PLK1***, ***AURKA***).
3. **Diagnóstico Espectral (PCA Scree Plot):**
   - Confirmada dominancia de $PC_1$ con el **$44.98\%$ de la varianza total** y dimensión efectiva $\text{ED} = 4.57$, evidenciando el cono de anisotropía característico de ESM-2.

---

## 2. Distribución Cuantitativa por Macro-Familia KEGG

| Macro-Familia KEGG (Opción A) | Genes Anotados ($N$) | Porcentaje (%) | Color en t-SNE | Procesos Biológicos Clave |
| :--- | :---: | :---: | :--- | :--- |
| **Expresión Génica & Traducción** | **510** | 26.2% | 🟤 Marrón (`#8c564b`) | Ribosomas (RPL/RPS), Espliceosoma (PRPFs, SNRPs), eIFs |
| **Metabolismo & Bioenergética** | **225** | 11.6% | 🟢 Verde (`#2ca02c`) | Glicólisis, TCA, OxPhos (ATP5/COX), Biosíntesis de aminoácidos |
| **Otras Rutas Biológicas** | **145** | 7.5% | 🔵 Cian (`#17becf`) | Transportadores ABC, tráfico vesicular, proteasoma |
| **Ciclo Celular & Apoptosis** | **109** | 5.6% | 🔴 Rojo (`#d62728`) | Checkpoints mitóticos, CDKs, ciclinas, caspasas |
| **Mantenimiento Genómico & Reparación ADN**| **85** | 4.4% | 🟣 Púrpura (`#9467bd`) | Reparación de ADN, recombinación homóloga, polimerasas |
| **Transducción de Señales & Quinasas** | **73** | 3.8% | 🟠 Naranja (`#ff7f0e`) | Vías MAPK, ErbB, PI3K/Akt, fosfatasas de señal |
| **Adhesión Celular & Citoesqueleto** | **20** | 1.0% | 🟡 Mostaza (`#bcbd22`) | Complejos de adhesión focal, filamentos de actina |
| **Sin anotación KEGG específica** | **777** | 40.0% | ⚪ Gris (`#c7c7c7`) | Dependencias esenciales sin vía metabólica/celular en KEGG |
| **TOTAL** | **1,944** | **100.0%** | — | **Universo de Dependencias DepMap Retenido** |

---

## 3. Trazabilidad de Genes Core de Cáncer de Mama

En estricta conformidad con la directriz de la investigación, se mantuvieron exclusivamente los genes resultantes del cribado determinístico:
- **`MYC`:** Presente en el panel de dependencias DepMap, mapeado y destacado prominentemente en el t-SNE.
- **`CCND1`, `MDM2`, `CDK4`, `AURKA`, `PLK1`:** Drivers oncológicos presentes en el panel, resaltados con marcadores de alto contraste y cajas de anotación.
- **Genes Excluidos Previamente (Omitidos deliberadamente):**
  - `PIK3CA`, `ERBB2`, `BRCA1`: Excluidos determinísticamente en `exp_02` ($L > 1022$ AA).
  - `TP53`, `ESR1`, `PTEN`, `CDH1`: Ausentes en el panel DepMap original por ser supresores tumorales o receptores específicos de tejido.

---

## 4. Artefactos Generados y Persistidos

| Artefacto | Ruta | Descripción |
| :--- | :--- | :--- |
| **Embeddings HDF5 ($N=1,944$)** | [`goals/01/data/esm2_650m_depmap_embeddings.h5`](file:///home/jmalagont/project/WSI_RNA_multimodal_learning/goals/01/data/esm2_650m_depmap_embeddings.h5) | $1,944 \times 1280$ en `float32` (SHA256: `9177f210...`) |
| **t-SNE 2D Macro-Familias KEGG** | [`goals/01/experiments/exp_04/results/tsne_depmap_2d.png`](file:///home/jmalagont/project/WSI_RNA_multimodal_learning/goals/01/experiments/exp_04/results/tsne_depmap_2d.png) | Gráfica t-SNE coloreada por KEGG con oncogenes destacados |
| **Asignaciones KEGG** | [`goals/01/experiments/exp_04/results/kegg_pathway_assignments.csv`](file:///home/jmalagont/project/WSI_RNA_multimodal_learning/goals/01/experiments/exp_04/results/kegg_pathway_assignments.csv) | Tabla completa de los 1,944 genes con su macro-familia y pathways |
| **Scree Plot PCA** | [`goals/01/experiments/exp_04/results/pca_scree_plot.png`](file:///home/jmalagont/project/WSI_RNA_multimodal_learning/goals/01/experiments/exp_04/results/pca_scree_plot.png) | Diagnóstico espectral de varianza individual y acumulada |
| **Estadísticas de Varianza** | [`goals/01/experiments/exp_04/results/pca_variance_stats.csv`](file:///home/jmalagont/project/WSI_RNA_multimodal_learning/goals/01/experiments/exp_04/results/pca_variance_stats.csv) | Top 100 componentes y autovalores |
| **Auditoría de Exclusiones** | [`goals/01/experiments/exp_04/results/depmap_filtration_summary.csv`](file:///home/jmalagont/project/WSI_RNA_multimodal_learning/goals/01/experiments/exp_04/results/depmap_filtration_summary.csv) | Auditoría de los 319 genes de DepMap fuera de la ventana $L \le 1022$ |
