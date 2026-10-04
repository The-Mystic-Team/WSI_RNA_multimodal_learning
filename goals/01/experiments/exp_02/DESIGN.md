# Experiment Design: exp_02 — Filtrado Determinístico y Curaduría de Secuencias Proteicas (L <= 1022 AA)

**Experiment**: `goals/01/experiments/exp_02/`  
**Project**: `WSI_RNA_multimodal_learning`  
**Goal**: `01`  
**Date**: 2026-10-03  
**Status**: Draft  

---

## 1. Hypothesis
El filtrado determinístico del proteoma humano recuperado (`uniprot_proteome.fasta`), restringiendo estrictamente a secuencias polipeptídicas de genes codificantes con longitud menor o igual a 1022 aminoácidos ($L \le 1022$ AA) para coincidir con la ventana contextual nativa libre de extrapolación posicional de ESM-2 (650M), retendrá al menos el $80\%$ del universo de genes del transcriptoma humano ($N=18,972$), generando como artefacto primario consolidado un archivo FASTA curado directamente en la ruta de datos del goal (`./data/filtered_proteome.fasta`), apto para inferencia estática sin riesgo de desbordamiento de memoria GPU (*CUDA OOM*) ni degradación de atención cuadrática.

---

## 2. Experimental Setup

- **Archivos de Entrada**:
  - Catálogo proteico completo obtenido en `exp_01`: `goals/01/data/uniprot_proteome.fasta` (referenciado como `./data/uniprot_proteome.fasta` desde el goal o `../../data/uniprot_proteome.fasta` desde el script; 18,330 secuencias polipeptídicas de isoforma más larga para *Homo sapiens*).
  - Matriz transcriptómica de referencia: `../../data/tcga-brca/preprocessed/expression/gen_expr_processed.csv` (18,972 genes HGNC del transcriptoma de TCGA-BRCA).
  - Reporte de mapeo previo: `goals/01/experiments/exp_01/results/mapping_stats.csv`.

- **Reglas Heurísticas y Filtros Determinísticos**:
  1. **Filtro de Codificación Proteica**:
     - Se excluyen todos los genes de la matriz transcriptómica que no posean una secuencia de proteína canónica anotada en UniProtKB (ej. ARN no codificantes lncRNA, microARNs, pseudogenes sin traducción demostrada).
     - Estado de exclusión asignado: `non_coding_or_unmapped`.
  2. **Filtro de Longitud Posicional Nominal ($L \le 1022$ AA)**:
     - Se descartan de forma estricta todas las secuencias cuya longitud sea mayor a 1022 aminoácidos ($L > 1022$ AA).
     - Justificación: el modelo `esm2_t33_650M_UR50D` posee un contexto nativo de 1024 tokens. Descontando los tokens sintácticos delimitadores (`<cls>` y `<eos>`), la capacidad máxima de procesamiento estandarizada es de 1022 residuos. Exceder este umbral demandaría extrapolación de *rotary position embeddings* fuera de distribución o generaría *out-of-memory* (CUDA OOM) por la complejidad de atención $O(L^2)$.
     - Estado de exclusión asignado: `exceeds_length_threshold`.
  3. **Verificación de Calidad y Validez Sintáctica Peptídica**:
     - Confirmación de que las secuencias contengan exclusivamente los 20 aminoácidos estándar de la IUPAC (A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y). En caso de existir residuos ambiguos (`X`, `B`, `Z`, `U`, `O`) o asteriscos de parada (`*`), se registrará la incidencia para tratamiento determinista.

- **Generación y Destino del Artefacto Primario (Ruta Data del Goal)**:
  - **Ruta de Salida**: `goals/01/data/filtered_proteome.fasta` (accesible localmente como `./data/filtered_proteome.fasta` desde la raíz del goal o `../../data/filtered_proteome.fasta` desde el directorio de ejecución de scripts de `exp_02`).
  - **Formato del Artefacto**: Archivo FASTA estándar con las secuencias curadas que superaron ambos filtros determinísticos.
  - **Formato de Encabezado**: `>GENE_SYMBOL|UNIPROT_ACCESSION|LENGTH_AA|Homo_sapiens` para preservar trazabilidad biunívoca con las columnas de expresión génica.
  - **Control de Integridad**: Generación automática de hash SHA-256 (`results/filtered_proteome.sha256`) tras la escritura del archivo para garantizar persistencia e inmutabilidad.

- **Entorno de Ejecución**: Conda `histo-DL` (Python 3.10+, pandas, biopython).

---

## 3. File Layout for This Experiment

```
goals/01/
├── data/
│   ├── uniprot_proteome.fasta         ← Entrada (exp_01): proteoma humano completo (18,330 seqs)
│   └── filtered_proteome.fasta        ← ARTEFACTO GENERADO (exp_02): FASTA filtrado (L <= 1022 AA)
└── experiments/
    └── exp_02/
        ├── DESIGN.md                  ← Este diseño experimental
        ├── IMPLEMENTATION.md          ← Plan de implementación (a aprobar antes de editar código)
        ├── scripts/                   ← Scripts de filtrado determinístico y verificación
        │   └── filter_proteome.py     ← Genera data/filtered_proteome.fasta y métricas
        ├── results/                   ← Artefactos tabulares de resultados
        │   ├── recovery_table.csv     ← Tabla formal de recuperación global y relativa
        │   ├── filtration_breakdown.csv ← Desglose cuantitativo de exclusiones
        │   ├── excluded_genes_summary.csv ← Listado detallado de exclusiones por gen
        │   └── filtered_proteome.sha256 ← Hash de verificación del FASTA generado en data/
        └── reports/
            ├── figures/
            │   └── filtering_waterfall.png ← Gráfico de cascada / retención secuencial
            └── summary.md             ← Reporte científico consolidado
```

---

## 4. Experimental Conditions

| Condición | Criterio de Inclusión | Criterio de Exclusión | Salida en `goals/01/data/` |
| :--- | :--- | :--- | :--- |
| `unfiltered_baseline` (exp_01) | Todos los genes mapeados en UniProtKB (*Homo sapiens*) | Genes sin producto proteico anotado | `./data/uniprot_proteome.fasta` ($N=18,330$) |
| `curated_esm2_nominal` (exp_02) | Genes con producto proteico Y $1 \le L \le 1022$ AA | 1. No codificantes ($N=642$)<br>2. Longitud $L > 1022$ AA | `./data/filtered_proteome.fasta` ($N \approx 15,939$) |

---

## 5. Evaluation Protocol & Deliverables

### A. Generación y Validación del Artefacto Primario (`goals/01/data/filtered_proteome.fasta`)
1. **Generación Determinista**:
   - Lectura de `goals/01/data/uniprot_proteome.fasta`.
   - Aplicación del filtro determinista $L \le 1022$ AA.
   - Escritura atómica en `goals/01/data/filtered_proteome.fasta`.
2. **Validación Sintáctica y Bioinformática**:
   - Conteo de registros en el FASTA generado debe ser exactamente igual a $N_{\text{passed}}$.
   - Ninguna secuencia en `filtered_proteome.fasta` puede tener $L > 1022$ AA ($\max(L) \le 1022$).
   - Todos los identificadores `GENE_SYMBOL` deben coincidir con columnas válidas del transcriptoma original.

### B. Artefactos de Resultados y Métricas Cuantitativas

1. **Tabla de Recuperación (`results/recovery_table.csv`)**:
   Deberá computar y tabular formalmente:
   - $N_{\text{total}}$: Total de genes en la matriz de transcriptoma ($18,972$).
   - $N_{\text{mapped}}$: Genes con producto proteico recuperado en UniProt ($18,330$).
   - $N_{\text{passed}}$: Proteínas que superan el filtro de longitud ($L \le 1022$ AA) incorporadas en `filtered_proteome.fasta`.
   - Tasa de Recuperación sobre el Transcriptoma Global:
     $$R_{\text{transcriptome}} = \frac{N_{\text{passed}}}{N_{\text{total}}} \times 100$$
   - Tasa de Retención sobre el Proteoma Mapeado:
     $$R_{\text{proteome}} = \frac{N_{\text{passed}}}{N_{\text{mapped}}} \times 100$$

2. **Tabla de Desglose de Exclusiones (`results/filtration_breakdown.csv`)**:
   - Total excluidos por no codificantes / sin mapeo ($N_{\text{non\_coding}}$).
   - Total excluidos por longitud excesiva ($N_{L > 1022}$).
   - Porcentajes respectivos sobre el total del transcriptoma ($N=18,972$).

3. **Listado Detallado de Exclusiones (`results/excluded_genes_summary.csv`)**:
   - Columnas: `gene_id`, `uniprot_acc`, `length_aa`, `exclusion_reason`.

4. **Visualización de Cascada de Filtrado (`reports/figures/filtering_waterfall.png`)**:
   - Diagrama de barras / cascada ilustrando la retención secuencial: Transcriptoma Inicial ($18,972$) $\to$ Proteoma Mapeado ($18,330$) $\to$ FASTA Curado en `data/` ($N_{\text{passed}}$).

### C. Criterios de Aceptación:
- **Tasa de Recuperación**: $R_{\text{transcriptome}} \ge 80\%$ (conforme a los resultados de `exp_01`, el valor esperado es exactamente **$84.01\%$**, equivalente a **$15,939$ genes**).
- **Consistencia Numérica**: $N_{\text{passed}} + N_{L > 1022} = N_{\text{mapped}}$ y $N_{\text{mapped}} + N_{\text{non\_coding}} = N_{\text{total}}$.
- **Presencia en `data/`**: Comprobación física de la existencia de `goals/01/data/filtered_proteome.fasta` y su hash SHA-256.

---

## 6. Risks & Mitigations

1. **Exclusión de Genes Oncológicos de Gran Tamaño**:
   - *Riesgo*: Genes oncológicos clave con secuencias canónicas $> 1022$ AA (ej. *BRCA1* [1,863 AA], *BRCA2* [3,418 AA], *ATM* [3,056 AA], *EGFR* [1,210 AA]) no formarán parte de `filtered_proteome.fasta`.
   - *Mitigación*: Documentar explícitamente en `reports/summary.md` y `results/excluded_genes_summary.csv` la lista de estos genes, justificando que la exclusión responde a la ventana física no extrapolada del modelo de lenguaje proteico.
2. **Inconsistencia de Rutas o Nombres de Archivo**:
   - *Riesgo*: Escribir el FASTA curado dentro del subdirectorio del experimento en lugar de la carpeta compartida del goal.
   - *Mitigación*: Fijar en el diseño e implementación la ruta estricta `goals/01/data/filtered_proteome.fasta` (o `../../data/filtered_proteome.fasta` desde el script).
3. **Residuos Ambigüos o Caracteres Inválidos**:
   - *Riesgo*: Presencia de selenocisteína (`U`) o códigos ambiguos (`X`) en las secuencias que originen problemas al tokenizar en `exp_03`.
   - *Mitigación*: Verificación exhaustiva en el script de filtrado para emitir alertas si se detectan caracteres no estándar.

---

## 7. Reproducibility Checklist

- [ ] Secuencia de entrada verificada en `goals/01/data/uniprot_proteome.fasta` (SHA-256: `314d3bcb...`).
- [ ] Generación directa del archivo `goals/01/data/filtered_proteome.fasta`.
- [ ] Hash SHA-256 del artefacto generado guardado en `results/filtered_proteome.sha256`.
- [ ] Cota de longitud fijada estrictamente en $L \le 1022$ aminoácidos sin parámetros dinámicos no controlados.
- [ ] Validación de consistencia con las 18,972 columnas del transcriptoma.
- [ ] Entorno de ejecución registrado bajo Conda `histo-DL`.

---

## 8. Next Steps

1. Revisión y ratificación de este diseño corregido (`DESIGN.md`).
2. Elaboración del plan de implementación en `experiments/exp_02/IMPLEMENTATION.md`.
3. Implementación del script generador de `filtered_proteome.fasta` en `./data/`, tabla de recuperación y figuras de reporte.
