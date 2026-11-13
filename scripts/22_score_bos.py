"""Score the frozen three-gene Bos set.

The genes were locked before these cohorts were scored. A 15-gene minimum
had marked the complete set invalid. A set is scored here when every locked
gene is present. GSE12276 is not an independent test of these genes.
"""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(r"g:\Cancer Metastasis\scripts")))
import importlib

score_mod = importlib.import_module("04_score_cohorts")

GENES = score_mod.read_genes("bos_genes.txt")
OUT = score_mod.OUT


def score_complete(expr, genes):
    present = [gene for gene in genes if gene in expr.index]
    missing = [gene for gene in genes if gene not in expr.index]
    if missing or not present:
        return None, {"n_present": len(present), "n_genes": len(genes), "missing": missing, "valid": False}
    sub = expr.loc[present]
    sd = sub.std(axis=1, ddof=0)
    kept = sub.loc[sd > 0]
    z = kept.sub(kept.mean(axis=1), axis=0).div(kept.std(axis=1, ddof=0), axis=0)
    return z.mean(axis=0), {
        "n_present": int(kept.shape[0]),
        "n_genes": len(genes),
        "missing": missing,
        "valid": int(kept.shape[0]) == len(genes),
    }


def main():
    report = {
        "genes": GENES,
        "rule": "mean of gene-wise z-scores when every locked gene is present",
        "inferential_cohort": "GSE76714",
        "circular_cohort": "GSE12276",
    }
    expr, meta = score_mod.prepare_expression("GSE76714")
    meta = meta.set_index("gsm")
    scores, coverage = score_complete(expr, GENES)
    groups = meta.loc[scores.index, "chars"].map(lambda chars: score_mod.char_map(chars).get("patient group", ""))
    y = groups.map(lambda value: 1 if "with brain metastasis" in value.lower() else 0 if "without brain metastasis" in value.lower() else np.nan)
    result = score_mod.test_binary(y, scores, "GSE76714:bos")
    result["coverage"] = coverage
    report["GSE76714"] = result

    expr, meta = score_mod.prepare_expression("GSE12276")
    meta = meta.set_index("gsm")
    scores, coverage = score_complete(expr, GENES)
    relapse = meta.loc[scores.index, "chars"].map(lambda chars: score_mod.char_map(chars).get("site of relapse (brain or other)", ""))
    y = relapse.map(lambda value: 1 if "brain" in value.lower() else 0 if value.lower() == "other" else np.nan)
    result = score_mod.test_binary(y, scores, "GSE12276:bos")
    result["coverage"] = coverage
    result["role"] = "pipeline check only; these genes are tied to this series"
    report["GSE12276"] = result
    (OUT / "bos_three_gene.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
