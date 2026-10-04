# Experiment Design: exp_03 — Extracción de Embeddings con ESM-2 (650M) y Reducción Dimensional t-SNE

**Experiment**: `goals/01/experiments/exp_03/`  
**Project**: `WSI_RNA_multimodal_learning`  
**Goal**: `01`  
**Date**: 2026-10-03  
**Status**: Approved  

---

## 1. Hypothesis
La extracción de representaciones densas mediante el modelo de lenguaje proteico preentrenado ESM-2 (`esm2_t33_650M_UR50D`, capa 33, $D=1280$) aplicando *mean pooling* estrictamente sobre los tokens de aminoácidos para las 15,940 secuencias polipeptídicas de *Homo sapiens* curadas en `exp_02` ($L \le 1022$ AA) generará una matriz de embeddings proteicos de dimensión $15,940 \times 1280$ libre de valores espurios (cero NaNs, cero Infs) y con variabilidad latente no degenerada (sin colapso dimensional), la cual quedará empaquetada como artefacto primario en formato HDF5 (`goals/01/data/esm2_650m_embeddings.h5`) indexada exclusivamente por identificadores génicos, y cuya proyección bidimensional vía t-SNE evidenciará una estructura topológica continua y agrupaciones asociadas a diversidad biológica en el espacio de representación.

---

## 2. Experimental Setup

- **Archivos de Entrada**:
  - Catálogo curado de secuencias polipeptídicas humanas: `goals/01/data/filtered_proteome.fasta` ($N=15,940$ secuencias con $L \le 1022$ AA, SHA256: `f65797d1e25d4be09d6ad146079d8cfd6bdc4b2cd96560da4e04279d84b4aecd`).
  - Identificadores esperados: Formato de encabezado `>GENE_SYMBOL|UNIPROT_ACCESSION|LENGTH_AA|Homo_sapiens`.

- **Modelo y Capa de Extracción**:
  - **Arquitectura**: ESM-2 650M (`facebook/esm2_t33_650M_UR50D`).
  - **Parámetros**: 650 millones, 33 capas de Transformer con atención rotatoria (RoPE).
  - **Dimensión de Embedding ($D$)**: 1280 características por residuo.
  - **Capa Objetivo**: Capa final de representaciones ocultas (capa 33 / *last hidden state*).
  - **Modo de Operación**: Inferencia pura (`model.eval()`, `torch.inference_mode()` / `torch.no_grad()`).

- **Estrategia de Pooling (Mean Pooling sobre Aminoácidos)**:
  - Dada una secuencia polipeptídica de longitud $L$, el tokenizador inserta delimitadores de inicio `<cls>` (índice 0) y fin `<eos>` (índice $L+1$).
  - La matriz de estados ocultos para la capa 33 es $H \in \mathbb{R}^{(L+2) \times 1280}$ (o con relleno de padding en caso de lotes).
  - El vector de representación del gen $e_{\text{gene}} \in \mathbb{R}^{1280}$ se calcula formalmente como el promedio de los residuos biológicos genuinos, excluyendo `<cls>`, `<eos>` y tokens de padding:
    $$e_{\text{gene}} = \frac{1}{L} \sum_{i=1}^{L} H_{i+1}$$
  - De esta forma, la representación agregada captura exclusivamente el contenido informacional de la cadena de aminoácidos sin sesgos de los tokens sintácticos del modelo de lenguaje.

- **Diseño del Artefacto Primario (`goals/01/data/esm2_650m_embeddings.h5`)**:
  - **Ruta**: `goals/01/data/esm2_650m_embeddings.h5` (accesible localmente como `./data/esm2_650m_embeddings.h5` desde la raíz del goal).
  - **Contenido del HDF5**: Contiene **únicamente** la identificación de cada gen y sus respectivos embeddings:
    - Dataset `gene_id` (o `genes`): Array 1D de cadenas de texto (formato string HDF5 / ASCII) con los 15,940 símbolos de genes HGNC ordenados de forma determinista.
    - Dataset `embeddings`: Matriz 2D de dimensiones $(15940, 1280)$ en precisión `float32`.
    - Atributos globales del archivo (metadatos de trazabilidad): `model_name="esm2_t33_650M_UR50D"`, `layer=33`, `pooling="amino_acid_mean_pooling"`, `input_fasta_sha256="f65797..."`, `created_date="2026-10-03"`.
  - **Integridad y Persistencia**: Generación de hash SHA-256 (`results/esm2_650m_embeddings.sha256`) tras la escritura completa del archivo HDF5.

- **Resultados Adicionales: Reducción Dimensional t-SNE 2D**:
  - **Algoritmo**: t-SNE (*t-Distributed Stochastic Neighbor Embedding*) aplicado sobre la matriz completa de $15,940 \times 1280$ (o tras inicialización PCA estándar con 50 componentes para aceleración y preservación global de varianza).
  - **Parámetros**: `n_components=2`, `perplexity=30`, `metric="cosine"` o `"euclidean"`, `random_state=42`.
  - **Visualización Entregable**: `reports/figures/tsne_embeddings.png` (gráfico de dispersión 2D de alta resolución con mapas de densidad para visualizar la distribución topológica del espacio latente proteico).

- **Hardware y Entorno de Ejecución**:
  - Entorno Conda: `histo-DL` (PyTorch 2.2.2+cu118, transformers 4.51.3, h5py 3.14.0, scikit-learn 1.9.0).
  - Acelerador: Inferencia en **una sola GPU** (NVIDIA GeForce GTX TITAN X, 12,288 MiB VRAM, Maxwell CC 5.2).
  - Precisión Numérica: `torch.float32` (precisión estándar óptima nativa para hardware Maxwell).
  - Estrategia de Batching: **Tamaño de lote fijo uniforme $B = 4$ secuencias por lote** (Opción A). Dado que $L \le 1022$ AA está estrictamente garantizado desde `exp_02`, un lote de 4 secuencias demanda un consumo pico de activaciones de ~3.0 GB, totalizando ~6.2 GB de VRAM con el modelo cargado, garantizando un margen libre de casi 6 GB y **cero riesgo de CUDA OOM**.

---

## 3. File Layout for This Experiment

```
goals/01/
├── data/
│   ├── filtered_proteome.fasta         ← Entrada (exp_02): 15,940 proteínas curadas (L <= 1022 AA)
│   └── esm2_650m_embeddings.h5         ← ARTEFACTO PRIMARIO (exp_03): HDF5 con gene_id y embeddings (15940 x 1280)
└── experiments/
    └── exp_03/
        ├── DESIGN.md                   ← Este diseño experimental
        ├── IMPLEMENTATION.md           ← Plan de implementación (a aprobar antes de editar código)
        ├── scripts/                    ← Scripts de inferencia y análisis dimensional
        │   ├── extract_embeddings.py   ← Inferencia ESM-2, mean-pooling y exportación HDF5
        │   └── analyze_tsne.py         ← Reducción t-SNE 2D y generación de figuras
        ├── results/                    ← Métricas y artefactos de validación
        │   ├── embedding_validation.csv← Normas, varianza, verificación NaN/Inf y rangos
        │   ├── extraction_summary.csv  ← Estadísticas de ejecución (tiempo, throughput tokens/s)
        │   └── esm2_650m_embeddings.sha256 ← Hash de verificación del HDF5 generado en data/
        └── reports/
            ├── figures/
            │   └── tsne_embeddings.png ← Proyección t-SNE 2D del espacio latente
            └── summary.md              ← Reporte científico consolidado
```

---

## 4. Baselines and Proposed Conditions

| Condición | Modelo / Método | Pooling | Salida Primaria en `goals/01/data/` | Validación Esperada |
| :--- | :--- | :--- | :--- | :--- |
| `esm2_650m_mean_pool` (Propuesta) | `esm2_t33_650M_UR50D` (Capa 33, FP32, $B=4$, 1x GPU) | Residue Mean Pooling ($1 \dots L$) | `./data/esm2_650m_embeddings.h5` | Matriz $(15940, 1280)$, 0 NaNs/Infs, varianza $> 0$, t-SNE estructurado |
| *Baseline de Control Sintáctico (Sanity Check)* | Matriz sintética o embeddings aleatorios | N/A | En memoria (verificación de contraste) | Verificar que las representaciones ESM-2 tienen similitud coseno no trivial frente a ruido blanco |

---

## 5. Evaluation Protocol & Deliverables

### A. Validación Numérica y Estructural del Artefacto HDF5
1. **Completitud y Cardinalidad**:
   - Total de registros: Exactamente $15,940$ filas en `gene_id` y en `embeddings`.
   - Correspondencia biunívoca con los identificadores presentes en `filtered_proteome.fasta`.
2. **Integridad Numérica Estricta**:
   - Cero valores `NaN` ($\sum \text{isnan} = 0$).
   - Cero valores infinitos ($\sum \text{isinf} = 0$).
   - Cero vectores nulos ($\|e_{\text{gene}}\|_2 > 0$ para todo gen).
3. **Distribución de Características Latentes**:
   - Distribución de normas $L_2$ de los vectores de embedding: media, desviación estándar, percentiles 5 y 95.
   - Varianza por dimensión: confirmación de que ninguna dimensión presenta varianza cero (ausencia de colapso representacional).
   - Rango efectivo y análisis espectral básico para verificar que la información está distribuida a lo largo de las 1280 dimensiones.

### B. Proyección y Visualización t-SNE 2D (`reports/figures/tsne_embeddings.png`)
1. **Cálculo t-SNE**:
   - Ejecución determinista sobre los 15,940 vectores normalizados o estandarizados con semilla fija (`seed=42`).
   - Generación de coordenadas bidimensionales $(z_1, z_2)$.
2. **Generación de Figura**:
   - Gráfico de dispersión 2D con visualización de densidad y gradiente de dispersión.
   - Evaluación visual de continuidad: verificación de formación de cúmulos y ausencia de artefactos de fragmentación anómala.

### C. Métricas de Rendimiento Computacional
- Tiempo total de inferencia y throughput promedio en secuencias/segundo y tokens/segundo.
- Consumo pico de VRAM registrado en la GPU GTX TITAN X durante la pasada hacia adelante.

---

## 6. Expected Results & Decision Rules

- **Criterio de Aceptación (Éxito)**:
  - Generación íntegra de `goals/01/data/esm2_650m_embeddings.h5` conteniendo exactamente 15,940 genes y matriz de forma $(15940, 1280)$.
  - 100% de valores numéricos válidos (0 NaNs, 0 Infs).
  - Hash SHA-256 generado y verificado.
  - Proyección t-SNE 2D guardada en `reports/figures/tsne_embeddings.png`.
- **Regla de Parada y Contingencia (Fallo)**:
  - Si ocurre un desbordamiento de memoria GPU (CUDA OOM), el script abortará de forma segura y se recalibrará el presupuesto de tokens por lote (`max_tokens_per_batch`).
  - Si se detecta algún NaN en los embeddings, se inspeccionará la secuencia causante (ej. caracteres desconocidos no filtrados) antes de consolidar el archivo.

---

## 7. Risks & Mitigations

| Riesgo | Probabilidad | Impacto | Mitigación |
| :--- | :--- | :--- | :--- |
| **CUDA OOM en secuencias largas ($L \approx 1000$ AA)** | Baja | Alto | Tamaño de lote fijo conservador $B = 4$ secuencias, con inferencia en `torch.inference_mode()`, aprovechando que la longitud máxima está acotada estrictamente a $L \le 1022$ AA. Consumo pico estimado $\approx 6.2$ GB sobre 11.9 GB de VRAM. |
| **Sesgo por inclusión de tokens `<cls>` / `<eos>` en el pooling** | Baja | Medio | Formulación explícita de máscara de atención que selecciona únicamente las posiciones de residuos biológicos genuinos ($1 \le i \le L$). |
| **Inconsistencia de orden de genes entre HDF5 y FASTA** | Baja | Crítico | Indexación secuencial preservando estrictamente el orden de iteración del FASTA original, registrando tanto `gene_id` como metadatos explícitos en el HDF5. |
| **Tiempo de inferencia prolongado (>60 min)** | Media | Bajo | Ejecución en segundo plano vía `tmux` persistente (siguiendo las reglas de `AGENTS.md`), permitiendo ejecución desacoplada sin riesgo de interrupción. |

---

## 8. Reproducibility Checklist

- [ ] Semillas fijadas (`seed=42`) para t-SNE y particionado de lotes.
- [ ] Versión del modelo registrada: `facebook/esm2_t33_650M_UR50D` (commit/revision hash de Hugging Face).
- [ ] Hash SHA-256 del archivo de entrada verificado (`filtered_proteome.fasta`).
- [ ] Hash SHA-256 del artefacto HDF5 generado registrado en `results/esm2_650m_embeddings.sha256`.
- [ ] Entorno congelado en Conda `histo-DL`.
- [ ] Código modularizado y parametrizado bajo `goals/01/experiments/exp_03/scripts/`.

---

## 9. Next Steps
1. Revisión y aprobación formal del usuario de este plan experimental (`DESIGN.md`).
2. Tras la aprobación, formulación del plan de implementación (`IMPLEMENTATION.md`) especificando los scripts de inferencia (`extract_embeddings.py`) y reducción t-SNE (`analyze_tsne.py`).
3. Ejecución de la inferencia en sesión persistente `tmux` de acuerdo con el protocolo de `AGENTS.md`.
