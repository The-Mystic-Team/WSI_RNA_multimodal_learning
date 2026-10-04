# Experiment Design: exp_01 — Descarga del Proteoma Humano vía UniProt API

**Experiment**: `goals/01/experiments/exp_01/`  
**Project**: `WSI_RNA_multimodal_learning`  
**Goal**: `01`  
**Date**: 2026-10-03  
**Status**: Draft  

---

## 1. Hypothesis
La consulta estructurada por lotes (*batch mapping*) a la API de UniProtKB/IdMapping a partir de los 18,972 identificadores HGNC del transcriptoma (`gen_expr_processed.csv`), restringida estrictamente a organismo *Homo sapiens* (Taxonomy ID: 9606) y seleccionando la **isoforma proteica más larga ($L_{\max}$)** por gen, recuperará de forma determinista las secuencias polipeptídicas de $\ge 90\%$ de los genes codificantes humanos anotados, produciendo un archivo FASTA desacoplado en `./data/uniprot_proteome.fasta` junto con un reporte exploratorio completo de cobertura y distribución de longitudes.

---

## 2. Experimental Setup
- **Fuente de Entrada**: Matriz de expresión génica preprocesada en `../../data/tcga-brca/preprocessed/expression/gen_expr_processed.csv`.
  - Universo total: 18,972 columnas correspondientes a símbolos oficiales de genes HGNC.
- **Servicio Remoto y Restricción Taxonómica**: 
  - UniProtKB REST API (`https://rest.uniprot.org/`).
  - **Filtro estricto de organismo**: ÚNICAMENTE *Homo sapiens* (`organism_id:9606`). Prohibida la incorporación de secuencias u homólogos de otros taxones.
- **Regla de Selección de Isoforma**:
  - Para cada gen que presente múltiples entradas o variantes de empalme alternativo (*alternative splicing isoforms*) en UniProtKB para *Homo sapiens*, se ordenan por longitud en aminoácidos y se extrae **únicamente la isoforma más larga ($L_{\max}$)** con secuencia completa.
- **Destino del Artefacto Primario**:
  - `goals/01/data/uniprot_proteome.fasta` (accesible como `./data/uniprot_proteome.fasta` desde el goal o `../../data/uniprot_proteome.fasta` desde el script del experimento).
- **Entorno de Ejecución**: Conda `histo-DL` (Python 3.10+, requests / urllib3 / matplotlib / pandas).

---

## 3. File Layout for This Experiment
```
goals/01/experiments/exp_01/
├── DESIGN.md                  ← Este documento de diseño científico
├── IMPLEMENTATION.md          ← Plan de implementación (a generar tras aprobación de este diseño)
├── scripts/                   ← Script(s) de extracción, selección de Lmax y análisis exploratorio
├── results/                   ← Métricas, logs de descarga y tablas cuantitativas
│   ├── recovery_table.csv     ← Tabla de porcentaje de recuperación y conteos
│   ├── mapping_stats.csv      ← Detalle por gen (HGNC, UniProt ID, Lmax, status)
│   └── unmapped_genes.txt     ← Lista de genes sin secuencia proteica
└── reports/
    ├── figures/
    │   └── sequence_length_histogram.png ← Histograma de distribución de tamaños (AA)
    └── summary.md             ← Resumen científico consolidado
```

---

## 4. Experimental Conditions
Se evaluará la recuperación de secuencias bajo una estrategia única de extracción en lotes con filtrado taxonómico e isoforma máxima:

| Condición | Organismo | Estrategia de Selección | Tamaño de Lote (Chunk) | Manejo de Reintentos |
| :--- | :--- | :--- | :--- | :--- |
| `batch_uniprot_lmax_human` | *Homo sapiens* (9606) estrictamente | Isoforma más larga ($L_{\max}$) por gen | 1,000 genes / petición | Backoff exponencial (3 reintentos por lote) |

---

## 5. Evaluation Protocol & Exploratory Results

### A. Resultados Exploratorios Obligatorios:
1. **Tabla de Porcentaje de Recuperación (`results/recovery_table.csv` & `reports/summary.md`)**:
   Deberá tabular formalmente:
   - Total de genes consultados en el transcriptoma ($N_{\text{total}} = 18,972$).
   - Total de genes mapeados exitosamente a isoforma humana ($N_{\text{mapped}}$).
   - Porcentaje de recuperación global:
     $$R_{\text{global}} = \frac{N_{\text{mapped}}}{N_{\text{total}}} \times 100$$
   - Total de genes no codificantes / no encontrados en UniProt ($N_{\text{unmapped}}$).
   - Porcentaje de pérdida:
     $$P_{\text{unmapped}} = \frac{N_{\text{unmapped}}}{N_{\text{total}}} \times 100$$
2. **Histograma de Distribución de Tamaños de Secuencia (`reports/figures/sequence_length_histogram.png`)**:
   - Histograma de frecuencias de las longitudes de los péptidos recuperados ($L_{\max}$ en número de aminoácidos).
   - Inclusión de líneas de referencia vertical: mediana, media y la cota teórica de 1022 AA (umbral de pre-entrenamiento de ESM-2 que se aplicará en `exp_02`).
   - Métricas de dispersión reportadas en texto: mínimo, $Q_1$, mediana, $Q_3$, percentil 95 y máximo ($L_{\max}^{\text{proteome}}$).

### B. Métricas Primarias y Criterios de Aceptación:
- **Tasa de Recuperación**: $R_{\text{coding}} \ge 90\%$ respecto a los genes codificantes humanos esperados (~16,500).
- **Pureza Taxonómica**: 100% de las entradas deben provenir de *Homo sapiens* (Taxonomy ID 9606).
- **Integridad Sintáctica del FASTA**:
  - Encabezado uniforme: `>GENE_SYMBOL|UNIPROT_ACCESSION|LENGTH_AA|Homo_sapiens`.
  - Secuencias no vacías, sin codones de parada ni residuos no estándar.

---

## 6. Risks & Mitigations
1. **Contaminación por Homólogos de Otras Especies**:
   - *Mitigación*: Filtrado explícito a nivel de query (`AND organism_id:9606`) y verificación secundaria post-descarga en los metadatos de UniProt.
2. **Límites de Petición de la API (HTTP 429)**:
   - *Mitigación*: Agrupación en lotes (*chunks* de 1,000 IDs) con pausas de 500 ms y backoff exponencial ante respuestas HTTP 429/503.
3. **Múltiples Isoformas por Gen**:
   - *Mitigación*: Algoritmo determinista de selección que agrupa entradas por símbolo HGNC, descarta fragmentos y toma $\arg\max_i (\text{len}(\text{seq}_i))$.
4. **Desconexiones de Red en Procesos Largos**:
   - *Mitigación*: Persistencia incremental (*checkpointing* local de lotes) para permitir reanudación transparente.

---

## 7. Reproducibility Checklist
- [ ] Versión fija del endpoint de la API registrado (`UniProt Release 2024_x`).
- [ ] Taxonomía fijada estrictamente a 9606 (*Homo sapiens*).
- [ ] Criterio de desempate determinista para secuencias de igual longitud (orden alfabético por UniProt Accession).
- [ ] Hash SHA-256 de `uniprot_proteome.fasta` registrado en `results/`.
- [ ] Script parametrizado ejecutable en sesión desacoplada de `tmux` (`01_exp_01`).

---

## 8. Next Steps
1. Revisión y aprobación del presente diseño ajustado (`DESIGN.md`).
2. Elaboración del plan de implementación en `experiments/exp_01/IMPLEMENTATION.md`.
3. Implementación de los scripts de consulta, extracción de $L_{\max}$ humana, generación de la tabla de recuperación e histograma de longitudes.
