"""Apply the DecontX EM to each exported sample and sum native mo-Mac counts."""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread

from decontx_em import decontx_em

ROOT = Path(r"g:\Cancer Metastasis")
IN_DIR = ROOT / "results" / "macrophage" / "ambient_input"
OUT_DIR = ROOT / "results" / "macrophage"
CLINICAL = ROOT / "results" / "locked" / "gse131907_momac_samples.csv"


def main():
    clinical = pd.read_csv(CLINICAL)
    columns = {}
    qc_rows = []
    for sample_id in clinical["sample"]:
        print(sample_id, flush=True)
        counts = mmread(IN_DIR / f"{sample_id}.mtx").tocsr()
        genes = (IN_DIR / f"{sample_id}.genes.txt").read_text(encoding="utf-8").splitlines()
        labels = pd.read_csv(IN_DIR / f"{sample_id}.labels.tsv", sep="\t")
        if counts.shape[1] != len(labels) or counts.shape[0] != len(genes):
            raise RuntimeError(f"{sample_id} matrix {counts.shape} does not match labels or genes")
        fit = decontx_em(counts, labels["label"].to_numpy())
        use = labels["momac"].to_numpy().astype(bool)
        native = fit["native"][:, use].sum(axis=1)
        columns[sample_id] = pd.Series(np.rint(native).astype(np.int64), index=genes)
        momac_contamination = fit["contamination"][use]
        qc_rows.append({
            "sample": sample_id,
            "iterations": fit["iterations"],
            "n_momac": int(use.sum()),
            "median_contamination": float(np.median(momac_contamination)),
            "mean_contamination": float(np.mean(momac_contamination)),
        })
        print(
            f"  iter {fit['iterations']} momac {int(use.sum())} "
            f"median_contam {np.median(momac_contamination):.3f}",
            flush=True,
        )
        del counts, fit
    pseudobulk = pd.DataFrame(columns).fillna(0).astype(np.int64)
    pseudobulk.index.name = "gene"
    pseudobulk.to_csv(OUT_DIR / "momac_ambient_pseudobulk.csv")
    pd.DataFrame(qc_rows).to_csv(OUT_DIR / "momac_ambient_contamination.csv", index=False)
    print("PSEUDOBULK", pseudobulk.shape)


if __name__ == "__main__":
    main()
