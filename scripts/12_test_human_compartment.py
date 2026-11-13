"""Sample-level test of the frozen compartment programs.

Primary test: two-sided Mann-Whitney U on sample scores.
A sample is eligible only when it has at least 20 malignant cells.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(r"g:\Cancer Metastasis")
OUT = ROOT / "results" / "compartment"
LOCKED = ROOT / "results" / "locked"
RNG = np.random.default_rng(20261112)


def read_symbols(path: Path):
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def load_expression():
    csf = pd.read_csv(OUT / "csf_pseudobulk_log2cpm.csv")
    par = pd.read_csv(OUT / "parenchyma_pseudobulk_log2cpm.csv")
    qc = pd.concat(
        [
            pd.read_csv(OUT / "csf_malignant_qc.csv"),
            pd.read_csv(OUT / "parenchyma_malignant_qc.csv"),
        ],
        ignore_index=True,
    )
    expr = pd.concat([csf, par], ignore_index=True, sort=False)
    expr = expr.merge(qc[["sample", "n_malignant", "malignant_umi"]], on="sample", how="left")
    eligible = expr[expr["n_malignant"] >= 20].copy()
    return eligible, qc


def score_program(eligible: pd.DataFrame, genes):
    present = [gene for gene in genes if gene in eligible.columns and eligible[gene].notna().all()]
    coverage = len(present) / len(genes) if genes else 0
    if coverage < 0.70:
        return None, present, coverage
    values = eligible[present].astype(float)
    z = (values - values.mean(axis=0)) / values.std(axis=0, ddof=0)
    z = z.replace([np.inf, -np.inf], np.nan).dropna(axis=1, how="any")
    if z.shape[1] / len(genes) < 0.70:
        return None, list(z.columns), z.shape[1] / len(genes)
    scores = z.mean(axis=1)
    return scores, list(z.columns), z.shape[1] / len(genes)


def test_scores(scores: pd.Series, labels: pd.Series, higher_group: str):
    csf = scores[labels == "CSF"].to_numpy()
    par = scores[labels == "Parenchyma"].to_numpy()
    stat, p = mannwhitneyu(csf, par, alternative="two-sided")
    # Permutation of the difference in means. The observed sign is retained.
    observed = float(csf.mean() - par.mean())
    pooled = np.concatenate([csf, par])
    n_csf = len(csf)
    null = np.empty(10000)
    for i in range(10000):
        draw = RNG.permutation(pooled)
        null[i] = draw[:n_csf].mean() - draw[n_csf:].mean()
    perm_p = float((np.sum(np.abs(null) >= abs(observed)) + 1) / (len(null) + 1))
    # Descriptive AUC for CSF versus parenchyma, with the higher_group coded as the positive class.
    y = (labels == higher_group).to_numpy()
    order = np.argsort(scores.to_numpy())
    ranks = np.empty(len(scores), dtype=float)
    ranks[order] = np.arange(1, len(scores) + 1)
    # Average ranks for ties.
    # Recompute with scipy-free tie handling via pandas rank.
    rank = scores.rank(method="average").to_numpy()
    n_pos = int(y.sum())
    n_neg = int(len(y) - n_pos)
    auc = float((rank[y == 1].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))
    return {
        "n_csf": int(len(csf)),
        "n_parenchyma": int(len(par)),
        "mean_csf": float(csf.mean()),
        "mean_parenchyma": float(par.mean()),
        "median_csf": float(np.median(csf)),
        "median_parenchyma": float(np.median(par)),
        "mannwhitney_u": float(stat),
        "mannwhitney_p": float(p),
        "mean_difference_csf_minus_parenchyma": observed,
        "permutation_p": perm_p,
        "auc_for_predicted_higher_group": auc,
        "predicted_higher_group": higher_group,
        "direction_matches_prediction": bool(
            (higher_group == "CSF" and observed > 0) or (higher_group == "Parenchyma" and observed < 0)
        ),
    }


def attach_control(eligible: pd.DataFrame, genes):
    extra = pd.concat(
        [
            pd.read_csv(OUT / "csf_proliferation_log2cpm.csv"),
            pd.read_csv(OUT / "parenchyma_proliferation_log2cpm.csv"),
        ],
        ignore_index=True,
        sort=False,
    )
    keep = ["sample"] + [gene for gene in genes if gene in extra.columns]
    eligible = eligible.drop(columns=[gene for gene in genes if gene in eligible.columns])
    return eligible.merge(extra[keep], on="sample", how="left")


def gene_level(eligible: pd.DataFrame, genes):
    present = [gene for gene in genes if gene in eligible.columns and eligible[gene].notna().all()]
    csf = eligible.loc[eligible["compartment"] == "CSF", present]
    par = eligible.loc[eligible["compartment"] == "Parenchyma", present]
    delta = csf.median(axis=0) - par.median(axis=0)
    return {
        "n_genes": len(present),
        "fraction_higher_in_csf": float((delta > 0).mean()) if len(present) else None,
        "median_log2cpm_difference_csf_minus_parenchyma": float(delta.median()) if len(present) else None,
    }


def sensitivity_min_cells(scores: pd.Series, labels: pd.Series, n_malignant: pd.Series, minimum: int, higher_group: str):
    keep = n_malignant >= minimum
    if int(keep.sum()) < 8 or int((labels[keep] == "CSF").sum()) < 4 or int((labels[keep] == "Parenchyma").sum()) < 4:
        return {"status": "too_few_samples", "minimum_malignant_cells": minimum}
    result = test_scores(scores[keep], labels[keep], higher_group)
    result["minimum_malignant_cells"] = minimum
    result["n_excluded"] = int((~keep).sum())
    return result


def main():
    eligible, qc = load_expression()
    prolif_genes = read_symbols(LOCKED / "proliferation_genes.txt")
    eligible = attach_control(eligible, prolif_genes)
    labels = eligible["compartment"]
    programs = {
        "leptomeningeal": (read_symbols(LOCKED / "leptomeningeal_genes_human.txt"), "CSF"),
        "parenchyma": (read_symbols(LOCKED / "parenchyma_genes_human.txt"), "Parenchyma"),
        "proliferation_control": (read_symbols(LOCKED / "proliferation_genes.txt"), "none"),
    }
    results = {"eligible_samples": eligible["sample"].tolist(), "qc": qc.to_dict(orient="records")}
    score_table = eligible[["sample", "compartment", "n_malignant"]].copy()
    for name, (genes, higher) in programs.items():
        scores, used, coverage = score_program(eligible, genes)
        entry = {"n_locked": len(genes), "n_used": len(used), "coverage": coverage}
        if scores is None:
            entry["status"] = "not_scored_coverage_below_0.70"
            results[name] = entry
            continue
        if higher == "none":
            higher = "CSF"
            entry["note"] = "No compartment prediction. AUC treats CSF as the positive class only so the number is comparable."
        entry.update(test_scores(scores, labels, higher))
        entry["gene_level"] = gene_level(eligible, genes)
        if higher != "none":
            entry["sensitivity_at_least_50_malignant_cells"] = sensitivity_min_cells(
                scores, labels, eligible["n_malignant"], 50, higher
            )
        entry["status"] = "tested"
        results[name] = entry
        score_table[name] = scores.to_numpy()
        print(name, json.dumps({k: entry[k] for k in entry if k != "note"}))
    score_table.to_csv(OUT / "human_sample_scores.csv", index=False)
    (OUT / "human_compartment_test.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    print("ELIGIBLE", len(eligible))


if __name__ == "__main__":
    main()
