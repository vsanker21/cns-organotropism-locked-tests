"""Proliferation-gene pseudobulk for the six CSF samples, same malignant-cell rule."""

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread

ROOT = Path(r"g:\Cancer Metastasis")
TENX = ROOT / "data" / "geo" / "GSE276139" / "10x"
OUT = ROOT / "results" / "compartment"

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


def main():
    genes = [
        line.strip()
        for line in (ROOT / "results" / "locked" / "proliferation_genes.txt").read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    rows = []
    for sample_id, prefix in SAMPLES:
        print(sample_id, flush=True)
        symbols = []
        with gzip.open(TENX / f"{prefix}-features.tsv.gz", "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                symbols.append(parts[1].upper() if len(parts) > 1 else parts[0].upper())
        index = {}
        for i, symbol in enumerate(symbols):
            index.setdefault(symbol, []).append(i)
        matrix = mmread(TENX / f"{prefix}-matrix.mtx.gz").tocsr()

        def total(wanted):
            acc = np.zeros(matrix.shape[1], dtype=np.float64)
            for gene in wanted:
                for row in index.get(gene, []):
                    acc += np.asarray(matrix[row, :].todense()).ravel()
            return acc

        malignant = (total(EPITHELIAL) >= 2) & (total(EPITHELIAL) > total(IMMUNE)) & (total(TCELL) == 0) & (total(MACROPHAGE) == 0)
        # Recompute the four sums once. The calls above scan repeatedly. Keep the mask.
        # Recompute cleanly below if the mask is empty.
        library = float(np.asarray(matrix[:, malignant].sum()))
        record = {"sample": sample_id}
        for gene in genes:
            positions = index.get(gene, [])
            if not positions or library == 0:
                record[gene] = np.nan
                continue
            count = sum(float(matrix[position, malignant].sum()) for position in positions)
            record[gene] = np.log2(count / library * 1e6 + 1)
        rows.append(record)
        del matrix
    pd.DataFrame(rows).to_csv(OUT / "csf_proliferation_log2cpm.csv", index=False)
    print("WROTE", len(rows))


if __name__ == "__main__":
    main()
