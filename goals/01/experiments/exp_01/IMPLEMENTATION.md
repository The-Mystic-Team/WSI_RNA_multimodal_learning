# Implementation Plan: exp_01 — Descarga del Proteoma Humano vía UniProt API

**Experiment ID**: `exp_01`  
**Goal**: `01`  
**Companion Design**: [`DESIGN.md`](./DESIGN.md)  
**Date**: 2026-10-03  
**Status**: Approved (Ready for Build)  

---

## 1. Overview & Operational Scope
Implementar un pipeline robusto y reproducible en Python para:
1. Extraer los 18,972 identificadores génicos oficiales (HGNC) desde la cabecera de `../../data/tcga-brca/preprocessed/expression/gen_expr_processed.csv`.
2. Consultar la API de UniProtKB restringida estrictamente a *Homo sapiens* (Taxonomy ID: 9606) para recuperar todas las secuencias peptídicas asociadas a los genes.
3. Para cada gen humano con múltiples entradas/isoformas, seleccionar deterministamente la **isoforma proteica más larga ($L_{\max}$)**.
4. Generar el artefacto primario en `./data/uniprot_proteome.fasta` (accesible a nivel de goal en `goals/01/data/uniprot_proteome.fasta`).
5. Producir los análisis exploratorios obligatorios:
   - Tabla de recuperación en `results/recovery_table.csv`.
   - Histograma de distribución de longitudes en `reports/figures/sequence_length_histogram.png`.
   - Listado de genes no codificantes / no mapeados en `results/unmapped_genes.txt`.

---

## 2. Directory & Artifact Structure

```
goals/01/
├── data/
│   └── uniprot_proteome.fasta     ← Artefacto primario de secuencias más largas
└── experiments/exp_01/
    ├── DESIGN.md
    ├── IMPLEMENTATION.md          ← Este documento
    ├── scripts/
    │   └── download_uniprot_proteome.py  ← Script autocontenido de extracción y reporte
    ├── results/
    │   ├── recovery_table.csv     ← Métricas globales de cobertura
    │   ├── mapping_stats.csv      ← Detalle por gen (HGNC, Accession, Lmax, Status)
    │   └── unmapped_genes.txt     ← Lista de genes sin secuencia humana
    └── reports/
        ├── figures/
        │   └── sequence_length_histogram.png ← Distribución de tamaños de secuencia
        └── summary.md             ← Reporte final generado post-ejecución
```

---

## 3. Algorithm & Component Design

### Fase A: Ingesta de Identificadores
- Leer únicamente la primera línea (cabecera) de `../../data/tcga-brca/preprocessed/expression/gen_expr_processed.csv` para evitar cargar la matriz densa en RAM innecesariamente.
- Descartar la columna `patient_id`.
- Obtener la lista ordenada y única de 18,972 símbolos génicos.

### Fase B: Recuperación de Secuencias desde UniProtKB (Estrategia Híbrida Eficiente)
Para garantizar máxima cobertura taxonómica y evitar saturación por 18,972 peticiones HTTP individuales:
1. **Descarga en Lote del Proteoma Canónico + Isoformas de *Homo sapiens* (UP000005640)**:
   - Stream masivo desde UniProt REST (`https://rest.uniprot.org/uniprotkb/stream?format=fasta&includeIsoform=true&query=(proteome:UP000005640)`).
   - Extraer metadatos de encabezado: accession, identificador de gen (`GN=...`), longitud y secuencia.
   - Indexar en un diccionario en memoria: `gene_symbol -> list of (accession, length, sequence)`.
2. **Fallback Dirigido para Genes No Resueltos**:
   - Para cualquier gen de la lista de 18,972 que no se encuentre en el proteoma de referencia principal (debido a alias o cambios de nomenclatura HGNC), ejecutar consultas en lotes (*batch query*) mediante el endpoint de búsqueda:
     `https://rest.uniprot.org/uniprotkb/search?query=(organism_id:9606)+AND+(gene_exact:<GENE>)+AND+reviewed:true&format=fasta`
   - Aplicar backoff exponencial (3 intentos) y pausas de 500 ms ante respuestas con código 429/503.

### Fase C: Selección Determinista de la Isoforma Más Larga ($L_{\max}$)
Para cada gen $g \in \{1, \dots, 18972\}$:
- Si existen secuencias humanas asociadas:
  - Ordenar por longitud descendente $\text{len}(\text{seq})$ y por orden alfabético de accession como desempate:
    $$(acc^*, seq^*) = \arg\max_{(acc, seq)} (\text{len}(seq), -acc)$$
  - Registrar en `mapping_stats`: `status="mapped"`, `uniprot_id=acc*`, `length=len(seq*)`.
  - Escribir en `uniprot_proteome.fasta` con encabezado:
    `>{gene_symbol}|{acc*}|{length}|Homo_sapiens`
- Si no existen secuencias:
  - Registrar en `mapping_stats`: `status="unmapped"`, `uniprot_id="None"`, `length=0`.
  - Añadir a `unmapped_genes.txt`.

### Fase D: Análisis Exploratorio y Métricas
1. **Cálculo de Cobertura (`recovery_table.csv`)**:
   - $N_{\text{total}} = 18,972$.
   - $N_{\text{mapped}} =$ conteo de genes con secuencia proteica humana.
   - $R_{\text{global}} = \frac{N_{\text{mapped}}}{N_{\text{total}}} \times 100$.
   - $N_{\text{unmapped}} =$ conteo de genes sin proteína.
   - $P_{\text{unmapped}} = \frac{N_{\text{unmapped}}}{N_{\text{total}}} \times 100$.
2. **Generación del Histograma (`sequence_length_histogram.png`)**:
   - Histograma de frecuencias de $L_{\max}$ (ancho de bins = 50 AA, rango 0 a $\max(L)$).
   - Línea vertical discontinua roja: mediana de la distribución.
   - Línea vertical punteada naranja: media.
   - Línea vertical sólida verde: umbral de 1022 AA (ventana nativa pre-entrenada de ESM-2 para `exp_02`).
   - Caja de anotación con estadísticos: $Q_1$, Mediana, Media, $Q_3$, P95, Máx, y porcentaje de secuencias $\le 1022$ AA.

---

## 4. Execution Directives & Telemetry

### Entorno y Dependencias
- Entorno Conda: `histo-DL`.
- Paquetes: `requests`, `urllib3`, `pandas`, `matplotlib`, `numpy`.

### Modo de Ejecución (Tmux Persistente - Fire-and-Forget)
Dado que la descarga y procesamiento del proteoma supera los 60 segundos, se ejecutará dentro de una sesión persistente de `tmux` bajo el identificador determinista `01_exp_01`:

```bash
tmux new-session -d -s 01_exp_01 \
  "source /home/jmalagont/miniconda3/etc/profile.d/conda.sh && conda activate histo-DL && \
   python -u goals/01/experiments/exp_01/scripts/download_uniprot_proteome.py; bash"
```

### Formato de Telemetría (Una línea por hito)
El script emitirá prints condensados en tiempo real:
```text
[Init] Cargados 18,972 identificadores HGNC desde gen_expr_processed.csv
[Stream] Descargando proteoma humano UP000005640 con isoformas...
[Progress] 42,000 entradas de UniProt procesadas e indexadas.
[Fallback] Consultando 312 genes sin coincidencia en proteoma principal...
[Selection] Selección de Lmax finalizada: 16,840 genes mapeados | 2,132 no codificantes.
[Artifacts] FASTA escrito en ../../data/uniprot_proteome.fasta (SHA256: ...)
[Plot] Histograma guardado en reports/figures/sequence_length_histogram.png
[Done] Exp 01 completado exitosamente. Tasa de recuperación: 88.76% global.
```

---

## 5. Verification & Acceptance Criteria
1. El archivo `goals/01/data/uniprot_proteome.fasta` existe, tiene tamaño no nulo (>20 MB), y cada entrada contiene una secuencia válida sin caracteres anómalos.
2. `results/recovery_table.csv` reporta cuantitativamente el porcentaje de recuperación.
3. `reports/figures/sequence_length_histogram.png` visualiza claramente la distribución con la línea de corte en 1022 AA.
4. `results/unmapped_genes.txt` lista todos los genes que no arrojaron secuencia proteica para su auditoría.
