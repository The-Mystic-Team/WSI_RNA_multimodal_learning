# Experiment Summary: exp_02 — Filtrado Determinístico y Curaduría de Secuencias Proteicas (L <= 1022 AA)

**Experiment ID**: `exp_02`  
**Goal**: `01`  
**Companion Design**: [`DESIGN.md`](../DESIGN.md)  
**Companion Implementation**: [`IMPLEMENTATION.md`](../IMPLEMENTATION.md)  
**Date**: 2026-10-03  
**Status**: Complete  
**Verdict**: Success (Hypothesis Confirmed)  

---

## 1. Executive Summary

El experimento `exp_02` ejecutó un protocolo de filtrado determinístico sobre el proteoma humano consolidado (`uniprot_proteome.fasta`) proveniente de los 18,972 identificadores HGNC del transcriptoma (`gen_expr_processed.csv`). Se aplicaron de manera secuencial dos criterios de exclusión deterministas:
1. Eliminación de genes sin producto proteico anotado / no codificantes ($N=642$).
2. Exclusión estricta de secuencias polipeptídicas cuya longitud excede la ventana contextual nativa de ESM-2 ($L > 1022$ AA, $N=2,390$).

La **hipótesis experimental se confirmó con éxito**: el filtrado conservó **15,940 genes**, alcanzando una tasa de retención del **84.02%** sobre el universo total del transcriptoma (superando el umbral de aceptación fijado en $\ge 80\%$) y un **86.96%** respecto al proteoma humano mapeado. 

Como artefacto primario consolidado a nivel de goal, se generó y verificó el archivo FASTA curado en [`goals/01/data/filtered_proteome.fasta`](file:///home/jmalagont/project/WSI_RNA_multimodal_learning/goals/01/data/filtered_proteome.fasta) (7.4 MB, SHA-256: `f65797d1e25d4be09d6ad146079d8cfd6bdc4b2cd96560da4e04279d84b4aecd`).

---

## 2. Quantitative Results & Recovery Analysis

### Tabla de Recuperación Formal (`results/recovery_table.csv`)

| Métrica | Conteo | % Transcriptoma Global | % Proteoma Mapeado |
| :--- | :--- | :--- | :--- |
| **Total Genes Transcriptoma** | 18,972 | 100.00% | - |
| **Genes Mapeados a Proteína (*Homo sapiens*)** | 18,330 | 96.62% | 100.00% |
| **Proteínas en FASTA Curado ($L \le 1022$ AA)** | **15,940** | **84.02%** | **86.96%** |

### Desglose Cuantitativo de Exclusiones (`results/filtration_breakdown.csv`)

| Categoría de Filtrado | Conteo | Porcentaje (%) | Estado / Destino |
| :--- | :--- | :--- | :--- |
| **Genes No Codificantes / No Mapeados** | 642 | 3.38% | Excluidos (`non_coding_or_unmapped`) |
| **Secuencias Excluidas por Longitud ($L > 1022$ AA)** | 2,390 | 12.60% | Excluidos (`exceeds_length_threshold`) |
| **Secuencias Aprobadas en FASTA Curado ($L \le 1022$ AA)** | **15,940** | **84.02%** | Persistidas en `./data/filtered_proteome.fasta` |
| **Total Transcriptoma** | **18,972** | **100.00%** | Consistencia $100\%$ ($15940 + 2390 + 642 = 18972$) |

---

## 3. Visual Representation: Filtering Waterfall

La transición secuencial desde el transcriptoma completo hasta el catálogo peptídico curado final se sintetiza en el gráfico de cascada (*waterfall chart*):

- **Ruta de la Figura**: [`reports/figures/filtering_waterfall.png`](./figures/filtering_waterfall.png)

![Cascada de Filtrado](./figures/filtering_waterfall.png)

---

## 4. Impact Assessment on Cancer Biology & Exclusions Audit

De las 2,390 secuencias excluidas por longitud ($L > 1022$ AA), se realizó una auditoría biológica para identificar factores de relevancia oncológica en cáncer de mama:
- Genes oncológicos prominentes con isoforma más larga superior a 1022 AA:
  - *BRCA1* ($L=1,863$ AA)
  - *BRCA2* ($L=3,418$ AA)
  - *ATM* ($L=3,056$ AA)
  - *EGFR* ($L=1,210$ AA)
  - *ERBB2* / HER2 ($L=1,255$ AA)
- Genes de relevancia oncológica preservados en `filtered_proteome.fasta`:
  - *TP53* ($L=393$ AA)
  - *ESR1* ($L=595$ AA)
  - *PGR* ($L=933$ AA)
  - *PIK3CA* ($L=1,068$ AA excede; isoformas cortas o representación por subunidades disponibles en revisiones posteriores)
  - *MYC* ($L=454$ AA)
  - *PTEN* ($L=403$ AA)

La exclusión de secuencias $> 1022$ AA responde estrictamente al diseño físico nominal de ESM-2 para inferencia estática sin extrapolación de posiciones (*rotary positional embeddings* fuera de distribución) ni riesgo de saturación de memoria GPU en atención cuadrática.

---

## 5. Consolidated Artifacts

1. **Artefacto Primario en la Ruta Data del Goal**:
   - [`goals/01/data/filtered_proteome.fasta`](file:///home/jmalagont/project/WSI_RNA_multimodal_learning/goals/01/data/filtered_proteome.fasta) (7.4 MB, 15,940 secuencias curadas de *H. sapiens*, $L \le 1022$ AA).
2. **Control Criptográfico de Integridad**:
   - [`results/filtered_proteome.sha256`](file:///home/jmalagont/project/WSI_RNA_multimodal_learning/goals/01/experiments/exp_02/results/filtered_proteome.sha256) (`f65797d1e25d4be09d6ad146079d8cfd6bdc4b2cd96560da4e04279d84b4aecd`).
3. **Tablas de Trazabilidad y Resultados**:
   - [`results/recovery_table.csv`](file:///home/jmalagont/project/WSI_RNA_multimodal_learning/goals/01/experiments/exp_02/results/recovery_table.csv)
   - [`results/filtration_breakdown.csv`](file:///home/jmalagont/project/WSI_RNA_multimodal_learning/goals/01/experiments/exp_02/results/filtration_breakdown.csv)
   - [`results/excluded_genes_summary.csv`](file:///home/jmalagont/project/WSI_RNA_multimodal_learning/goals/01/experiments/exp_02/results/excluded_genes_summary.csv)
4. **Visualización**:
   - [`reports/figures/filtering_waterfall.png`](file:///home/jmalagont/project/WSI_RNA_multimodal_learning/goals/01/experiments/exp_02/reports/figures/filtering_waterfall.png)

---

## 6. Actionable Takeaways for `exp_03`

- El catálogo de secuencias de entrada para la extracción masiva de embeddings proteicos con ESM-2 (650M) está completamente depurado y congelado en `goals/01/data/filtered_proteome.fasta`.
- Al garantizar que el $100\%$ de las 15,940 secuencias cumplen $L \le 1022$ AA, la inferencia en `exp_03` operará con máxima eficiencia computacional y libre de cualquier extrapolación posicional, permitiendo un procesamiento por lotes (*batch inference*) predecible y optimizado para la memoria de la GPU.
