"""Malignant-cell pseudobulk for the six GSE276139 CSF samples.

The gene lists and the cell rule are read from results/locked. This script does not
change either list.
"""

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread

ROOT = Path(r"g:\Cancer Metastasis")
TENX = ROOT / "data" / "geo" / "GSE276139" / "10x"
OUT = ROOT / "results" / "compartment"
OUT.mkdir(parents=True, exist_ok=True)

EPITHELIAL = ["KRT8", "KRT18", "KRT19", "EPCAM", "NAPSA", "NKX2-1"]
IMMUNE = ["PTPRC", "CD3D", "CD3E", "CD68", "AIF1", "C1QA"]
TCELL = ["CD3D", "CD3E"]
MACROPHAGE = ["CD68", "AIF1", "C1QA"]

SAMPLES = [
    ("CSF_CSFLUAD1", "GSM8491849_CSFLUAD1"),
    ("CSF_CSFLUAD2", "GSM8491850_CSFLUAD2"),
    ("CSF_CSFLUAD3", "GSM8491851_CSFLUAD3"),
    ("CSF_CWB", "GSM9505951_CWB"),
    ("CSF_WLC", "GSM9505952_WLC"),
    ("CSF_ZRQ", "GSM9505953_ZRQ"),
]


def read_symbols(path: Path):
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def features_for(prefix: str):
    path = TENX / f"{prefix}-features.tsv.gz"
    symbols = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            symbols.append(parts[1].upper() if len(parts) > 1 else parts[0].upper())
    return symbols


def row_index(symbols):
    """Sum duplicate symbols. Returns gene -> integer row positions."""
    index = {}
    for i, symbol in enumerate(symbols):
        index.setdefault(symbol, []).append(i)
    return index


def sum_rows(matrix, index, genes):
    total = np.zeros(matrix.shape[1], dtype=np.float64)
    for gene in genes:
        rows = index.get(gene)
        if not rows:
            continue
        for row in rows:
            total += np.asarray(matrix[row, :].todense()).ravel()
    return total


def main():
    lepto = read_symbols(ROOT / "results" / "locked" / "leptomeningeal_genes_human.txt")
    paren = read_symbols(ROOT / "results" / "locked" / "parenchyma_genes_human.txt")
    genes = list(dict.fromkeys(lepto + paren))
    rows = []
    qc = []
    for sample_id, prefix in SAMPLES:
        print(sample_id, flush=True)
        symbols = features_for(prefix)
        index = row_index(symbols)
        matrix = mmread(TENX / f"{prefix}-matrix.mtx.gz").tocsr()
        if matrix.shape[0] != len(symbols):
            raise RuntimeError(f"{sample_id} matrix rows {matrix.shape[0]} != features {len(symbols)}")
        epithelial = sum_rows(matrix, index, EPITHELIAL)
        immune = sum_rows(matrix, index, IMMUNE)
        tcell = sum_rows(matrix, index, TCELL)
        macrophage = sum_rows(matrix, index, MACROPHAGE)
        malignant = (epithelial >= 2) & (epithelial > immune) & (tcell == 0) & (macrophage == 0)
        n_mal = int(malignant.sum())
        library = np.asarray(matrix[:, malignant].sum(axis=0)).ravel() if n_mal else np.array([])
        total_umi = float(library.sum()) if n_mal else 0.0
        record = {
            "sample": sample_id,
            "compartment": "CSF",
            "n_cells": int(matrix.shape[1]),
            "n_malignant": n_mal,
            "malignant_umi": total_umi,
        }
        qc.append(record)
        print(f"  cells {matrix.shape[1]} malignant {n_mal} umi {total_umi:.0f}", flush=True)
        if n_mal == 0 or total_umi == 0:
            continue
        gene_sums = {}
        for gene in genes:
            positions = index.get(gene, [])
            if not positions:
                gene_sums[gene] = np.nan
                continue
            total = 0.0
            for position in positions:
                total += float(matrix[position, malignant].sum())
            gene_sums[gene] = total
        cpm = {gene: (np.nan if not np.isfinite(value) else np.log2(value / total_umi * 1e6 + 1)) for gene, value in gene_sums.items()}
        cpm["sample"] = sample_id
        cpm["compartment"] = "CSF"
        rows.append(cpm)
        del matrix
    pd.DataFrame(qc).to_csv(OUT / "csf_malignant_qc.csv", index=False)
    pd.DataFrame(rows).to_csv(OUT / "csf_pseudobulk_log2cpm.csv", index=False)
    print("CSF_SAMPLES_WRITTEN", len(rows))


if __name__ == "__main__":
    main()
