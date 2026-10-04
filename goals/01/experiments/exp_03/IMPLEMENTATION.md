# Implementation Plan: exp_03 — Extracción de Embeddings con ESM-2 (650M) y Reducción Dimensional t-SNE

**Experiment ID**: `exp_03`  
**Goal**: `01`  
**Companion Design**: [`DESIGN.md`](./DESIGN.md)  
**Date**: 2026-10-03  
**Status**: Approved (Ready for Execution)  

---

## 1. Overview & Operational Scope

Implementar un pipeline de inferencia profunda y reducción dimensional en Python para:
1. Cargar las 15,940 secuencias de proteínas humanas curadas en `exp_02` desde `goals/01/data/filtered_proteome.fasta`.
2. Ejecutar inferencia por lotes estáticos ($B=4$) en una sola GPU (NVIDIA GTX TITAN X, 12 GB VRAM) empleando el modelo `facebook/esm2_t33_650M_UR50D` en modo evaluación (`torch.inference_mode()`, `float32`).
3. Aplicar **mean pooling** estrictamente sobre los residuos de aminoácidos ($1 \dots L$), excluyendo los tokens sintácticos `<cls>`, `<eos>` y cualquier relleno de *padding*.
4. **Generar y persistir el artefacto primario** en la ruta de datos del goal:
   - `goals/01/data/esm2_650m_embeddings.h5` (accesible como `../../data/esm2_650m_embeddings.h5` desde el script o `./data/esm2_650m_embeddings.h5` desde el goal), conteniendo **únicamente** `gene_id` (15,940 símbolos HGNC) y `embeddings` (matriz $[15940, 1280]$ en `float32`).
5. Calcular y registrar el hash de integridad SHA-256 en `results/esm2_650m_embeddings.sha256`.
6. Validar numéricamente la ausencia de anomalías (0 NaNs, 0 Infs, varianza $>0$, normas $L_2$) en `results/embedding_validation.csv` y registrar telemetría de cómputo en `results/extraction_summary.csv`.
7. Ejecutar reducción no lineal t-SNE en 2 dimensiones sobre la matriz latente completa y generar la figura de alta resolución en `reports/figures/tsne_embeddings.png`.

---

## 2. Directory & Artifact Structure

```
goals/01/
├── data/
│   ├── filtered_proteome.fasta         ← Entrada (exp_02): 15,940 proteínas curadas (L <= 1022 AA)
│   └── esm2_650m_embeddings.h5         ← ARTEFACTO PRIMARIO GENERADO (exp_03) (15940 x 1280, float32)
└── experiments/exp_03/
    ├── DESIGN.md                       ← Diseño experimental ratificado
    ├── IMPLEMENTATION.md               ← Este documento
    ├── scripts/
    │   ├── extract_embeddings.py       ← Inferencia ESM-2 650M, pooling y empaquetado HDF5
    │   └── analyze_tsne.py             ← Reducción t-SNE 2D y generación de visualizaciones
    ├── results/
    │   ├── esm2_650m_embeddings.sha256 ← Hash de verificación del HDF5 generado en data/
    │   ├── embedding_validation.csv    ← Verificación de NaNs, normas L2, varianza y rangos
    │   └── extraction_summary.csv      ← Tiempo total, throughput (seqs/s, tokens/s), memoria VRAM
    └── reports/
        ├── figures/
        │   └── tsne_embeddings.png     ← Gráfico 2D t-SNE de distribución de embeddings
        └── summary.md                  ← Reporte científico consolidado post-ejecución
```

---

## 3. Algorithm & Component Design

### Componente 1: Inferencia y Empaquetado HDF5 (`scripts/extract_embeddings.py`)

#### Fase 1: Ingesta del FASTA y Preprocesamiento
- Lectura secuencial de `../../data/filtered_proteome.fasta`.
- Parseo de los encabezados `>GENE_SYMBOL|UNIPROT_ACCESSION|LENGTH_AA|Homo_sapiens` extrayendo el `GENE_SYMBOL` y la cadena polipeptídica.
- Verificación de consistencia: Total de secuencias = $15,940$, longitud máxima comprobada $\le 1022$ AA.

#### Fase 2: Configuración del Modelo y Dispositivo
- Selección explícita de dispositivo: `cuda:0` (NVIDIA GeForce GTX TITAN X).
- Carga de tokenizador y modelo desde Hugging Face / Transformers:
  - `AutoTokenizer.from_pretrained("facebook/esm2_t33_650M_UR50D")`
  - `AutoModel.from_pretrained("facebook/esm2_t33_650M_UR50D", torch_dtype=torch.float32)`
- Activación de `model.eval()` e inferencia bajo `torch.inference_mode()`.

#### Fase 3: Inferencia por Lotes ($B=4$) y Mean Pooling de Residuos
- Iteración en lotes fijos de 4 proteínas ($B=4$).
- Tokenización con padding dinámico por lote hasta la secuencia más larga del lote (`padding=True`, `return_tensors="pt"`).
- Para cada secuencia $j$ en el lote, cuya longitud biológica real es $L_j$:
  - Los estados ocultos de la última capa (capa 33) tienen forma `(B, seq_len_padded, 1280)`.
  - El token `<cls>` ocupa el índice 0; los residuos aminoacídicos ocupan los índices $1$ a $L_j$; el token `<eos>` ocupa el índice $L_j + 1$.
  - El embedding a nivel de gen $e_j \in \mathbb{R}^{1280}$ se calcula formalmente como:
    $$e_j = \frac{1}{L_j} \sum_{i=1}^{L_j} H_{j, i, :}$$
- Conversión a NumPy (`float32`).
- Telemetría en tiempo real: log condensado de progreso en una sola línea cada 500 secuencias procesadas (indicando índice, porcentaje completado, throughput de tokens/s y tiempo transcurrido).

#### Fase 4: Persistencia en HDF5 (`../../data/esm2_650m_embeddings.h5`)
- Creación del archivo HDF5 utilizando `h5py`:
  - Dataset `gene_id`: Array 1D de cadenas de texto (tipo `h5py.string_dtype(encoding='utf-8')`) de longitud 15,940 con los símbolos de genes en orden exacto.
  - Dataset `embeddings`: Matriz 2D de dimensiones $(15940, 1280)$ en precisión `float32`.
  - Atributos del archivo: `model_name="facebook/esm2_t33_650M_UR50D"`, `layer=33`, `pooling="residue_mean_pooling"`, `source_fasta="filtered_proteome.fasta"`.
- Cierre seguro del descriptor de archivo HDF5.

#### Fase 5: Validación Numérica y Cómputo de Hashes
- Lectura de verificación sobre el HDF5 recién creado:
  - Comprobación de dimensiones exactas: $(15940, 1280)$.
  - Conteo de `NaN` e `Inf`: debe ser estrictamente 0.
  - Cálculo de estadísticas de normas $L_2$ (mínimo, mediana, media, percentiles 5 y 95, máximo).
  - Cálculo de varianza por dimensión (confirmar que las 1280 dimensiones tienen varianza $> 0$).
  - Guardar resumen en `results/embedding_validation.csv`.
- Cálculo del hash SHA-256 de `../../data/esm2_650m_embeddings.h5` y persistencia en `results/esm2_650m_embeddings.sha256`.
- Guardar resumen de telemetría de inferencia (tiempo total, secuencias/s, consumo pico de VRAM) en `results/extraction_summary.csv`.

---

### Componente 2: Reducción Dimensional t-SNE 2D (`scripts/analyze_tsne.py`)

#### Fase 1: Carga de Embeddings
- Apertura en modo lectura de `../../data/esm2_650m_embeddings.h5`.
- Carga de la matriz de $(15940, 1280)$ y la lista de genes.

#### Fase 2: Reducción Dimensional Determinista
- Estandarización / normalización de características (opcional para t-SNE, o reducción previa con PCA a 50 componentes para optimizar tiempo y preservar estructura global).
- Ajuste y transformación con `TSNE`:
  - `n_components=2`
  - `perplexity=30`
  - `metric="cosine"` (adecuada para espacios de alta dimensión normalizados)
  - `random_state=42`
  - `init="pca"`
- Coordenadas bidimensionales resultantes: $Z \in \mathbb{R}^{15940 \times 2}$.

#### Fase 3: Renderizado y Exportación de Visualización
- Generación de figura con `matplotlib`:
  - Gráfico de dispersión 2D (`reports/figures/tsne_embeddings.png`) con mapa de densidad (hexbin o scatter con `alpha=0.35` y color por densidad kernel).
  - Anotación de metadatos: Modelo (ESM-2 650M), $N=15,940$, pooling (residue mean), perplejidad = 30.
  - Exportación en alta resolución (300 DPI).

---

## 4. Execution & Compute Protocol (Conforme a AGENTS.md)

Dado que la inferencia de 15,940 secuencias en una sola GPU GTX TITAN X tomará más de 60 segundos (~20–40 minutos):
1. **Ambiente Conda**: `histo-DL`.
2. **Sesión Persistente tmux**:
   ```bash
   tmux new-session -d -s 01_exp_03 \
     "source /home/jmalagont/miniconda3/etc/profile.d/conda.sh && conda activate histo-DL && \
      cd goals/01/experiments/exp_03 && \
      python -u scripts/extract_embeddings.py && \
      python -u scripts/analyze_tsne.py; bash"
   ```
3. **Monitoreo desacoplado (Fire-and-forget)**:
   - El agente no esperará pasivamente el proceso.
   - Para inspeccionar el progreso en tiempo real:
     ```bash
     tmux attach -t 01_exp_03
     ```
   - Desconectar de la sesión tmux: `Ctrl+B`, luego `D`.

---

## 5. Verification & Acceptance Criteria

| Artefacto / Requisito | Criterio de Aceptación |
| :--- | :--- |
| `goals/01/data/esm2_650m_embeddings.h5` | Existe en `data/` del goal; tamaño esperado ~80 MB; contiene datasets `gene_id` ($N=15,940$) y `embeddings` ($15,940 \times 1280$, `float32`). |
| Integridad Numérica | Cero `NaN`, cero `Inf`, varianza $> 0$ en todas las 1280 dimensiones (`results/embedding_validation.csv`). |
| Control Criptográfico | Hash SHA-256 verificado y registrado en `results/esm2_650m_embeddings.sha256`. |
| Figura t-SNE | Imagen generada en `reports/figures/tsne_embeddings.png` (300 DPI, dispersión 2D estructurada sin artefactos numéricos). |
| Telemetría | Resumen completo en `results/extraction_summary.csv`. |

---

## 6. Next Step
Una vez que el usuario apruebe este plan de implementación, se procederá a:
1. Crear los scripts `extract_embeddings.py` y `analyze_tsne.py` en `goals/01/experiments/exp_03/scripts/`.
2. Solicitar aprobación explícita del usuario para iniciar la ejecución en `tmux`.
