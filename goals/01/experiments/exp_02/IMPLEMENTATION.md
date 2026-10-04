# Implementation Plan: exp_02 — Filtrado Determinístico y Curaduría de Secuencias Proteicas (L <= 1022 AA)

**Experiment ID**: `exp_02`  
**Goal**: `01`  
**Companion Design**: [`DESIGN.md`](./DESIGN.md)  
**Date**: 2026-10-03  
**Status**: Approved (Ready for Build)  

---

## 1. Overview & Operational Scope

Implementar un pipeline determinístico y reproducible en Python para:
1. Cargar el catálogo completo de proteínas humanas generado en `exp_01` (`goals/01/data/uniprot_proteome.fasta`).
2. Cruzar con el universo transcriptómico de referencia (`../../data/tcga-brca/preprocessed/expression/gen_expr_processed.csv`) y los metadatos de mapeo (`../exp_01/results/mapping_stats.csv`).
3. Aplicar las dos reglas de filtrado determinístico:
   - Exclusión de genes no codificantes / sin producto proteico.
   - Exclusión de secuencias con longitud mayor a 1022 aminoácidos ($L > 1022$ AA).
4. **Generar y persistir el artefacto primario** directamente en la ruta de datos del goal:
   - `goals/01/data/filtered_proteome.fasta` (accesible localmente como `./data/filtered_proteome.fasta` desde el goal o `../../data/filtered_proteome.fasta` desde el script).
5. Calcular el hash de integridad SHA-256 de `filtered_proteome.fasta` y guardarlo en `results/filtered_proteome.sha256`.
6. Generar las tablas y gráficos de resultados cuantitativos:
   - `results/recovery_table.csv`: Métricas de retención sobre el transcriptoma y el proteoma.
   - `results/filtration_breakdown.csv`: Desglose cuantificado de exclusiones (no codificantes vs. $L > 1022$ AA).
   - `results/excluded_genes_summary.csv`: Lista completa de genes descartados con motivo y longitud.
   - `reports/figures/filtering_waterfall.png`: Diagrama de cascada / retención secuencial de genes.

---

## 2. Directory & Artifact Structure

```
goals/01/
├── data/
│   ├── uniprot_proteome.fasta         ← Entrada: Proteoma completo de H. sapiens (18,330 seqs)
│   └── filtered_proteome.fasta        ← ARTEFACTO PRIMARIO GENERADO (exp_02) (L <= 1022 AA)
└── experiments/exp_02/
    ├── DESIGN.md                      ← Diseño experimental ratificado
    ├── IMPLEMENTATION.md              ← Este documento
    ├── scripts/
    │   └── filter_proteome.py         ← Script autocontenido de filtrado y generación de artefactos
    ├── results/
    │   ├── recovery_table.csv         ← Métricas formales de recuperación
    │   ├── filtration_breakdown.csv   ← Desglose de motivos de descarte
    │   ├── excluded_genes_summary.csv ← Detalle exhaustivo de genes excluidos
    │   └── filtered_proteome.sha256   ← Hash SHA-256 de data/filtered_proteome.fasta
    └── reports/
        ├── figures/
        │   └── filtering_waterfall.png ← Gráfico de cascada de retención
        └── summary.md                 ← Resumen científico post-ejecución
```

---

## 3. Algorithm & Component Design (`scripts/filter_proteome.py`)

### Fase A: Ingesta y Validación de Fuentes
- Leer los identificadores del transcriptoma desde la cabecera de `../../data/tcga-brca/preprocessed/expression/gen_expr_processed.csv` (18,972 genes HGNC).
- Cargar `../exp_01/results/mapping_stats.csv` para identificar directamente el subconjunto de genes no mapeados/no codificantes ($N=642$).
- Parsear `../../data/uniprot_proteome.fasta` registrando para cada entrada: `gene_symbol`, `uniprot_acc`, `length_aa`, `sequence`.

### Fase B: Filtrado Determinístico y Generación de `filtered_proteome.fasta`
1. **Partición de Secuencias**:
   - Para cada entrada en `uniprot_proteome.fasta`:
     - Si $\text{length} \le 1022$:
       - Clasificación: `passed`.
       - Formatear encabezado: `>{gene_symbol}|{uniprot_acc}|{length_aa}|Homo_sapiens`.
       - Añadir a la lista de retención para `filtered_proteome.fasta`.
     - Si $\text{length} > 1022$:
       - Clasificación: `exceeds_length_threshold`.
       - Registrar en la lista de exclusiones con su longitud exacta.
2. **Escritura del Artefacto Primario en la Ruta Data del Goal**:
   - Abrir en modo escritura `../../data/filtered_proteome.fasta`.
   - Escribir las secuencias retenidas en formato FASTA estándar (líneas de secuencia envueltas a 60 u 80 caracteres para legibilidad y estándar bioinformático).
   - Calcular inmediatamente el hash SHA-256 del archivo generado y persistirlo en `results/filtered_proteome.sha256`.

### Fase C: Integración de Metadatos de Exclusión
- Consolidar todos los genes excluidos en un DataFrame unificado (`results/excluded_genes_summary.csv`):
  - Genes sin mapeo / no codificantes (`non_coding_or_unmapped`): 642 registros.
  - Genes con proteína $> 1022$ AA (`exceeds_length_threshold`): ~2,391 registros.

### Fase D: Computación de Tablas de Resultados
1. **Tabla de Recuperación (`results/recovery_table.csv`)**:
   | Metrica | Conteo | Porcentaje_Transcriptoma | Porcentaje_Proteoma |
   | :--- | :--- | :--- | :--- |
   | Total Genes Transcriptoma | 18972 | 100.00% | - |
   | Genes Mapeados a Proteína (Homo sapiens) | 18330 | 96.62% | 100.00% |
   | Proteínas que Superan Filtro (L <= 1022 AA) | 15939 | 84.01% | 86.96% |

2. **Desglose de Filtrado (`results/filtration_breakdown.csv`)**:
   | Categoria | Conteo | Porcentaje |
   | :--- | :--- | :--- |
   | Genes No Codificantes / No Mapeados | 642 | 3.38% |
   | Excluidos por Longitud (> 1022 AA) | 2391 | 12.60% |
   | Aprobados en FASTA Curado (L <= 1022 AA) | 15939 | 84.01% |
   | Total Transcriptoma | 18972 | 100.00% |

### Fase E: Generación del Diagrama de Cascada (`reports/figures/filtering_waterfall.png`)
- Generar mediante `matplotlib` un gráfico de cascada (*waterfall chart*) de alta resolución (300 DPI):
  - Barra 1 (Azul): Transcriptoma TCGA-BRCA ($18,972$).
  - Barra 2 (Rojo): Excluidos No Codificantes ($-642$).
  - Barra 3 (Gris): Proteoma Humano Mapeado ($18,330$).
  - Barra 4 (Naranja): Excluidos Longitud $> 1022$ AA ($-2,391$).
  - Barra 5 (Verde): FASTA Curado Final ($15,939$, $84.01\%$).
  - Anotaciones numéricas y porcentuales sobre cada barra.

---

## 4. Execution Directives & Verification

### Entorno y Dependencias
- Entorno Conda: `histo-DL`.
- Dependencias: `python`, `pandas`, `matplotlib`, `hashlib`.

### Régimen de Ejecución
- Debido a que el procesamiento de 18,330 registros en memoria y generación de gráficos toma $< 10$ segundos, la validación se puede ejecutar en primer plano (< 60s) conforme a las directivas de `AGENTS.md`.

### Criterios de Aceptación Post-Ejecución
1. El archivo `goals/01/data/filtered_proteome.fasta` existe físicamente y no está vacío.
2. El conteo de encabezados en `goals/01/data/filtered_proteome.fasta` coincide exactamente con el valor $N_{\text{passed}}$ de `results/recovery_table.csv`.
3. Ninguna secuencia dentro de `filtered_proteome.fasta` excede los 1022 aminoácidos.
4. El archivo `results/filtered_proteome.sha256` contiene el hash criptográfico exacto del FASTA generado.
5. Los archivos de resultados y figuras se ubican en sus respectivas carpetas del experimento.

---

## 5. Gatekeeping & Next Steps
- Solicitar aprobación explícita del usuario para crear `goals/01/experiments/exp_02/scripts/filter_proteome.py` y ejecutar el script.
