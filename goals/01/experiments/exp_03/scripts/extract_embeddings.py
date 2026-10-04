#!/usr/bin/env python3
"""
exp_03: Inferencia ESM-2 650M, Mean Pooling de Residuos y Exportación HDF5.
Procesa las 15,940 secuencias de filtered_proteome.fasta en lotes de tamaño 4 sobre una GPU.
"""

import time
import hashlib
from pathlib import Path
import numpy as np
import h5py
import torch
import pandas as pd
from transformers import AutoTokenizer, AutoModel


def parse_fasta(fasta_path):
    records = []
    current_header = None
    current_seq = []
    with open(fasta_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header:
                    seq_str = "".join(current_seq)
                    gene_symbol = current_header.split("|")[0]
                    records.append((gene_symbol, seq_str))
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            gene_symbol = current_header.split("|")[0]
            records.append((gene_symbol, "".join(current_seq)))
    return records


def main():
    start_total_time = time.time()
    script_dir = Path(__file__).resolve().parent
    exp_dir = script_dir.parent
    goal_dir = exp_dir.parent.parent

    input_fasta = goal_dir / "data/filtered_proteome.fasta"
    output_h5 = goal_dir / "data/esm2_650m_embeddings.h5"
    results_dir = exp_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    output_h5.parent.mkdir(parents=True, exist_ok=True)

    print(f"[{time.strftime('%H:%M:%S')}] Iniciando exp_03: Extracción de embeddings ESM-2 650M")
    print(f"[{time.strftime('%H:%M:%S')}] Archivo de entrada: {input_fasta}")

    records = parse_fasta(input_fasta)
    n_total = len(records)
    print(f"[{time.strftime('%H:%M:%S')}] Cargadas {n_total:,} secuencias de proteínas curadas")

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    print(f"[{time.strftime('%H:%M:%S')}] Dispositivo de cómputo: {device} ({torch.cuda.get_device_name(0)})")

    model_name = "facebook/esm2_t33_650M_UR50D"
    print(f"[{time.strftime('%H:%M:%S')}] Cargando tokenizador y modelo: {model_name}")
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModel.from_pretrained(model_name, torch_dtype=torch.float32)
    model.eval()
    model.to(device)
    print(f"[{time.strftime('%H:%M:%S')}] Modelo cargado en GPU. VRAM inicial: {torch.cuda.memory_allocated(0)/(1024**2):.1f} MB")

    batch_size = 4
    embedding_dim = 1280
    all_embeddings = np.zeros((n_total, embedding_dim), dtype=np.float32)
    all_gene_ids = []

    total_tokens = 0
    t0 = time.time()

    print(f"[{time.strftime('%H:%M:%S')}] Procesando en lotes fijos de {batch_size} con mean pooling de residuos...")

    for i in range(0, n_total, batch_size):
        batch_records = records[i:i + batch_size]
        batch_genes = [r[0] for r in batch_records]
        batch_seqs = [r[1] for r in batch_records]

        # Tokenización dinámica por lote con padding al elemento más largo del lote
        inputs = tokenizer(batch_seqs, return_tensors="pt", padding=True)
        inputs = {k: v.to(device) for k, v in inputs.items()}

        with torch.inference_mode():
            outputs = model(**inputs)
            hidden_states = outputs.last_hidden_state  # forma: (B, L_padded + 2, 1280)

        # Mean pooling exclusivamente sobre los residuos biológicos genuinos
        for j, seq in enumerate(batch_seqs):
            seq_len = len(seq)
            total_tokens += seq_len
            # Residuo 1 hasta L (índice 1 a seq_len + 1, ignorando <cls>=0 y <eos>=seq_len+1)
            aa_rep = hidden_states[j, 1:seq_len + 1, :]
            mean_vector = aa_rep.mean(dim=0).cpu().numpy()
            idx = i + j
            all_embeddings[idx] = mean_vector
            all_gene_ids.append(batch_genes[j])

        # Telemetría condensada cada 500 secuencias o al finalizar
        processed = min(i + batch_size, n_total)
        if processed % 500 == 0 or processed == n_total:
            elapsed = time.time() - t0
            throughput_seq = processed / elapsed
            throughput_tok = total_tokens / elapsed
            vram_mb = torch.cuda.max_memory_allocated(0) / (1024**2)
            pct = (processed / n_total) * 100.0
            print(
                f"[{time.strftime('%H:%M:%S')}] Avance: {processed:5d}/{n_total} ({pct:5.1f}%) | "
                f"Throughput: {throughput_seq:4.1f} seq/s ({throughput_tok:6.0f} tok/s) | "
                f"VRAM pico: {vram_mb:5.0f} MB"
            )

    inference_duration = time.time() - start_total_time
    print(f"[{time.strftime('%H:%M:%S')}] Inferencia completada en {inference_duration:.1f} s ({inference_duration/60:.2f} min)")

    # 1. Guardar artefacto HDF5 primario en goals/01/data/
    print(f"[{time.strftime('%H:%M:%S')}] Exportando artefacto HDF5 a {output_h5}...")
    with h5py.File(output_h5, "w") as f:
        dt = h5py.string_dtype(encoding="utf-8")
        f.create_dataset("gene_id", data=np.array(all_gene_ids, dtype=object), dtype=dt)
        f.create_dataset("embeddings", data=all_embeddings, dtype=np.float32)
        f.attrs["model_name"] = model_name
        f.attrs["layer"] = 33
        f.attrs["pooling"] = "residue_mean_pooling"
        f.attrs["total_genes"] = n_total
        f.attrs["embedding_dim"] = embedding_dim
        f.attrs["created_at"] = time.strftime("%Y-%m-%d %H:%M:%S")

    # 2. Computar hash SHA-256 de esm2_650m_embeddings.h5
    sha256_hash = hashlib.sha256()
    with open(output_h5, "rb") as f:
        for byte_block in iter(lambda: f.read(65536), b""):
            sha256_hash.update(byte_block)
    h5_sha256 = sha256_hash.hexdigest()
    sha256_file = results_dir / "esm2_650m_embeddings.sha256"
    with open(sha256_file, "w", encoding="utf-8") as f:
        f.write(f"{h5_sha256}  esm2_650m_embeddings.h5\n")
    print(f"[{time.strftime('%H:%M:%S')}] SHA-256 generado: {h5_sha256}")

    # 3. Validación numérica estricta
    print(f"[{time.strftime('%H:%M:%S')}] Validando integridad numérica...")
    nan_count = int(np.isnan(all_embeddings).sum())
    inf_count = int(np.isinf(all_embeddings).sum())
    norms = np.linalg.norm(all_embeddings, axis=1)
    variances = np.var(all_embeddings, axis=0)
    zero_var_dims = int((variances == 0).sum())

    val_df = pd.DataFrame([{
        "total_genes": n_total,
        "embedding_dim": embedding_dim,
        "nan_count": nan_count,
        "inf_count": inf_count,
        "zero_variance_dims": zero_var_dims,
        "mean_norm_l2": float(np.mean(norms)),
        "std_norm_l2": float(np.std(norms)),
        "min_norm_l2": float(np.min(norms)),
        "max_norm_l2": float(np.max(norms)),
        "p05_norm_l2": float(np.percentile(norms, 5)),
        "p95_norm_l2": float(np.percentile(norms, 95)),
        "mean_dimension_variance": float(np.mean(variances)),
        "min_dimension_variance": float(np.min(variances)),
    }])
    val_file = results_dir / "embedding_validation.csv"
    val_df.to_csv(val_file, index=False)
    print(f"[{time.strftime('%H:%M:%S')}] Validación guardada en {val_file}")

    # 4. Resumen de extracción
    summary_df = pd.DataFrame([{
        "model": model_name,
        "batch_size": batch_size,
        "total_sequences": n_total,
        "total_tokens": total_tokens,
        "inference_seconds": round(inference_duration, 2),
        "sequences_per_second": round(n_total / inference_duration, 2),
        "tokens_per_second": round(total_tokens / inference_duration, 2),
        "peak_vram_mb": round(torch.cuda.max_memory_allocated(0) / (1024**2), 1),
        "sha256": h5_sha256
    }])
    summary_file = results_dir / "extraction_summary.csv"
    summary_df.to_csv(summary_file, index=False)
    print(f"[{time.strftime('%H:%M:%S')}] Telemetría guardada en {summary_file}")
    print(f"[{time.strftime('%H:%M:%S')}] Extracción completada exitosamente.")


if __name__ == "__main__":
    main()
