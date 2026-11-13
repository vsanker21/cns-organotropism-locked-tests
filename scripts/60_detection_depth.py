"""Median genes detected per malignant cell, for a post hoc covariate.

A gene is detected when a Gene Expression feature for that symbol has a raw
count above zero. Duplicate symbols are counted once. The malignant-cell
rule is the one already used for the compartment scores. This file does not
change a locked gene list.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

ROOT = Path(__file__).resolve().parents[1]
TENX = ROOT / "data" / "geo" / "GSE276139" / "10x"
OUT = ROOT / "results" / "posthoc"
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


def features_for(prefix: str) -> tuple[list[str], list[str]]:
    symbols, kinds = [], []
    with gzip.open(TENX / f"{prefix}-features.tsv.gz", "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            symbols.append(parts[1].upper() if len(parts) > 1 else parts[0].upper())
            kinds.append(parts[2] if len(parts) > 2 else "Gene Expression")
    return symbols, kinds


def row_index(symbols: list[str]) -> dict[str, list[int]]:
    index: dict[str, list[int]] = {}
    for i, symbol in enumerate(symbols):
        index.setdefault(symbol, []).append(i)
    return index


def sum_rows(matrix: sparse.spmatrix, index: dict[str, list[int]], genes: list[str]) -> np.ndarray:
    total = np.zeros(matrix.shape[1], dtype=np.float64)
    for gene in genes:
        for row in index.get(gene, []):
            total += np.asarray(matrix[row, :].todense()).ravel()
    return total


def genes_detected(matrix: sparse.csr_matrix, symbols: list[str], gene_rows: np.ndarray) -> np.ndarray:
    """Unique symbols with a positive count, one value per cell."""
    sub = matrix[gene_rows, :].tocsr()
    if sub.data.size and np.any(sub.data <= 0):
        sub = sub.copy()
        sub.data = np.where(sub.data > 0, 1, 0)
        sub.eliminate_zeros()
    used = np.array([symbols[int(row)] for row in gene_rows], dtype=object)
    rows, cols = sub.nonzero()
    if rows.size == 0:
        return np.zeros(sub.shape[1], dtype=np.float64)
    codes, _ = pd.factorize(used, sort=False)
    symbol_code = codes[rows]
    order = np.lexsort((symbol_code, cols))
    symbol_code = symbol_code[order]
    cols = cols[order]
    keep = np.empty(cols.size, dtype=bool)
    keep[0] = True
    keep[1:] = (cols[1:] != cols[:-1]) | (symbol_code[1:] != symbol_code[:-1])
    return np.bincount(cols[keep], minlength=sub.shape[1]).astype(np.float64)


def parenchyma_median_genes() -> None:
    """GSE131907 mBrain. One gene row at a time, so the full frame is not copied."""
    annotation = pd.read_csv(
        ROOT / "data" / "geo" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz",
        sep="\t",
    )
    brain = annotation.loc[annotation["Sample_Origin"] == "mBrain", ["Index", "Sample"]]
    sample_of = dict(zip(brain["Index"], brain["Sample"]))
    matrix_path = ROOT / "data" / "geo" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    with gzip.open(matrix_path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    cells = header[1:]
    positions = np.array([i for i, name in enumerate(cells) if name in sample_of], dtype=int)
    if len(positions) != len(sample_of):
        raise RuntimeError(f"matched {len(positions)} mBrain cells, annotation has {len(sample_of)}")
    samples = np.array([sample_of[cells[i]] for i in positions])
    n_cells = len(positions)
    detected = np.zeros(n_cells, dtype=np.int32)
    epithelial = np.zeros(n_cells, dtype=np.float64)
    immune = np.zeros(n_cells, dtype=np.float64)
    tcell = np.zeros(n_cells, dtype=np.float64)
    macrophage = np.zeros(n_cells, dtype=np.float64)
    first_positive: dict[str, np.ndarray] = {}
    genes_read = 0
    with gzip.open(matrix_path, "rt", encoding="utf-8", errors="replace") as handle:
        handle.readline()
        for line in handle:
            symbol, values = line.split("\t", 1)
            symbol = symbol.upper()
            counts = np.fromstring(values, sep="\t", dtype=np.int32)
            if counts.size != len(cells):
                raise RuntimeError(f"{symbol} has {counts.size} cells, header has {len(cells)}")
            picked = counts[positions]
            positive = picked > 0
            previous = first_positive.get(symbol)
            if previous is None:
                detected += positive
                first_positive[symbol] = positive
            else:
                detected += positive & ~previous
                previous |= positive
            if symbol in EPITHELIAL:
                epithelial += picked
            if symbol in IMMUNE:
                immune += picked
            if symbol in TCELL:
                tcell += picked
            if symbol in MACROPHAGE:
                macrophage += picked
            genes_read += 1
            if genes_read % 4000 == 0:
                print(f"  genes {genes_read}", flush=True)
            del counts, picked, positive
    print(f"  genes {genes_read}", flush=True)
    del first_positive
    malignant = (epithelial >= 2) & (epithelial > immune) & (tcell == 0) & (macrophage == 0)
    expected = pd.read_csv(ROOT / "results" / "compartment" / "parenchyma_malignant_qc.csv")
    rows = []
    for sample_id in sorted(set(samples)):
        use = malignant & (samples == sample_id)
        n_mal = int(use.sum())
        stored_name = f"BrM_{sample_id}"
        stored = int(expected.loc[expected["sample"] == stored_name, "n_malignant"].iloc[0])
        if n_mal != stored:
            raise RuntimeError(f"{stored_name} malignant {n_mal} != stored {stored}")
        values = detected[use]
        rows.append(
            {
                "sample": stored_name,
                "compartment": "Parenchyma",
                "n_malignant": n_mal,
                "median_genes": float(np.median(values)),
                "min_genes": float(values.min()) if n_mal else np.nan,
                "max_genes": float(values.max()) if n_mal else np.nan,
            }
        )
        print(f"{stored_name} malignant {n_mal} median_genes {rows[-1]['median_genes']:.1f}", flush=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "parenchyma_median_genes.csv", index=False)
    print("PARENCHYMA_DEPTH_WRITTEN", len(frame), flush=True)


def main() -> None:
    expected = pd.read_csv(ROOT / "results" / "compartment" / "csf_malignant_qc.csv").set_index("sample")
    rows = []
    for sample_id, prefix in SAMPLES:
        print(sample_id, flush=True)
        symbols, kinds = features_for(prefix)
        index = row_index(symbols)
        matrix = mmread(TENX / f"{prefix}-matrix.mtx.gz").tocsr()
        if matrix.shape[0] != len(symbols):
            raise RuntimeError(f"{sample_id} rows {matrix.shape[0]} != features {len(symbols)}")
        epithelial = sum_rows(matrix, index, EPITHELIAL)
        immune = sum_rows(matrix, index, IMMUNE)
        tcell = sum_rows(matrix, index, TCELL)
        macrophage = sum_rows(matrix, index, MACROPHAGE)
        malignant = (epithelial >= 2) & (epithelial > immune) & (tcell == 0) & (macrophage == 0)
        n_mal = int(malignant.sum())
        if n_mal != int(expected.loc[sample_id, "n_malignant"]):
            raise RuntimeError(f"{sample_id} malignant {n_mal} != stored {int(expected.loc[sample_id, 'n_malignant'])}")
        gene_rows = np.array([i for i, kind in enumerate(kinds) if kind == "Gene Expression"], dtype=int)
        detected = genes_detected(matrix, symbols, gene_rows)[malignant]
        rows.append(
            {
                "sample": sample_id,
                "compartment": "CSF",
                "n_malignant": n_mal,
                "median_genes": float(np.median(detected)),
                "min_genes": float(detected.min()),
                "max_genes": float(detected.max()),
            }
        )
        print(
            f"  malignant {n_mal} median_genes {rows[-1]['median_genes']:.1f}",
            flush=True,
        )
        del matrix
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "csf_median_genes.csv", index=False)
    print("CSF_DEPTH_WRITTEN", len(frame))


if __name__ == "__main__":
    main()
