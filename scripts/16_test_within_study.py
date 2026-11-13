"""Sample-level tests inside GSE131907. Gene lists stay frozen.

Contrasts, fixed in results/locked/within_study_rule.txt:
PE versus solid metastases (brain and lymph node together), and brain versus lymph node.
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


def score_program(frame: pd.DataFrame, genes):
    present = [gene for gene in genes if gene in frame.columns and frame[gene].notna().all()]
    coverage = len(present) / len(genes) if genes else 0.0
    if coverage < 0.70:
        return None, present, coverage
    values = frame[present].astype(float)
    z = (values - values.mean(axis=0)) / values.std(axis=0, ddof=0)
    z = z.replace([np.inf, -np.inf], np.nan).dropna(axis=1, how="any")
    if z.shape[1] / len(genes) < 0.70:
        return None, list(z.columns), z.shape[1] / len(genes)
    return z.mean(axis=1), list(z.columns), z.shape[1] / len(genes)


def compare(scores: pd.Series, labels: pd.Series, group_a: str, group_b: str):
    a = scores[labels == group_a].to_numpy()
    b = scores[labels == group_b].to_numpy()
    stat, p = mannwhitneyu(a, b, alternative="two-sided")
    observed = float(a.mean() - b.mean())
    pooled = np.concatenate([a, b])
    null = np.empty(10000)
    for i in range(10000):
        draw = RNG.permutation(pooled)
        null[i] = draw[: len(a)].mean() - draw[len(a) :].mean()
    perm_p = float((np.sum(np.abs(null) >= abs(observed)) + 1) / (len(null) + 1))
    return {
        "group_a": group_a,
        "group_b": group_b,
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "mean_a": float(a.mean()),
        "mean_b": float(b.mean()),
        "median_a": float(np.median(a)),
        "median_b": float(np.median(b)),
        "mean_difference_a_minus_b": observed,
        "mannwhitney_u": float(stat),
        "mannwhitney_p": float(p),
        "permutation_p": perm_p,
    }


def bh(p_values):
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    n = len(p)
    q = np.empty(n)
    running = 1.0
    for i in range(n - 1, -1, -1):
        running = min(running, ranked[i] * n / (i + 1))
        q[i] = running
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


def main():
    qc = pd.read_csv(OUT / "gse131907_within_qc.csv")
    expr = pd.read_csv(OUT / "gse131907_within_log2cpm.csv")
    expr = expr.merge(qc[["sample", "n_malignant", "malignant_umi"]], on="sample", how="left")
    eligible = expr[expr["n_malignant"] >= 20].copy()
    eligible["site"] = eligible["origin"].map({"PE": "Effusion", "mBrain": "Brain", "mLN": "LymphNode"})
    eligible["solid"] = np.where(eligible["origin"] == "PE", "Effusion", "Solid")

    programs = {
        "leptomeningeal": read_symbols(LOCKED / "leptomeningeal_genes_human.txt"),
        "parenchyma": read_symbols(LOCKED / "parenchyma_genes_human.txt"),
        "g2m": read_symbols(LOCKED / "proliferation_genes.txt"),
    }
    score_table = eligible[["sample", "origin", "site", "n_malignant"]].copy()
    results = {
        "eligible_by_origin": eligible.groupby("origin").size().astype(int).to_dict(),
        "qc": qc.to_dict(orient="records"),
        "tests": [],
    }
    stored = {}
    for name, genes in programs.items():
        scores, used, coverage = score_program(eligible, genes)
        stored[name] = scores
        if scores is None:
            results["tests"].append({"program": name, "status": "not_scored", "coverage": coverage, "n_locked": len(genes)})
            continue
        score_table[name] = scores.to_numpy()
        for contrast, labels, a, b in (
            ("effusion_vs_solid", eligible["solid"], "Effusion", "Solid"),
            ("brain_vs_lymph_node", eligible["site"], "Brain", "LymphNode"),
        ):
            if contrast == "brain_vs_lymph_node":
                mask = eligible["origin"].isin(["mBrain", "mLN"])
                entry = compare(scores[mask], labels[mask], a, b)
            else:
                entry = compare(scores, labels, a, b)
            entry.update({"program": name, "contrast": contrast, "n_locked": len(genes), "n_used": len(used), "coverage": coverage, "status": "tested"})
            results["tests"].append(entry)

    tested = [row for row in results["tests"] if row.get("status") == "tested"]
    q_values = bh([row["mannwhitney_p"] for row in tested])
    for row, q in zip(tested, q_values):
        row["fdr_among_six"] = float(q)

    by_key = {(row["program"], row["contrast"]): row for row in tested}

    def fluid_supported(program):
        row = by_key.get((program, "effusion_vs_solid"))
        return bool(row and row["mean_difference_a_minus_b"] > 0 and row["mannwhitney_p"] < 0.05)

    def separates_brain(program):
        row = by_key.get((program, "brain_vs_lymph_node"))
        return bool(row and row["mannwhitney_p"] < 0.05)

    def lepto_beats(contrast):
        lepto = by_key.get(("leptomeningeal", contrast))
        g2m = by_key.get(("g2m", contrast))
        if lepto is None or g2m is None:
            return False
        opposite = np.sign(lepto["mean_difference_a_minus_b"]) == -np.sign(g2m["mean_difference_a_minus_b"])
        return bool(
            opposite
            and lepto["mannwhitney_p"] < 0.05
            and abs(lepto["mean_difference_a_minus_b"]) > abs(g2m["mean_difference_a_minus_b"])
            and np.sign(g2m["mean_difference_a_minus_b"]) != 0
        )

    results["reading"] = {
        "fluid_versus_solid_supported": fluid_supported("parenchyma") and fluid_supported("g2m"),
        "parenchyma_higher_in_effusion_and_significant": fluid_supported("parenchyma"),
        "g2m_higher_in_effusion_and_significant": fluid_supported("g2m"),
        "parenchyma_brain_specific": separates_brain("parenchyma"),
        "g2m_separates_brain_from_lymph_node": separates_brain("g2m"),
        "leptomeningeal_beats_g2m_on_effusion_contrast": lepto_beats("effusion_vs_solid"),
        "leptomeningeal_beats_g2m_on_brain_contrast": lepto_beats("brain_vs_lymph_node"),
    }
    score_table.to_csv(OUT / "gse131907_within_scores.csv", index=False)
    (OUT / "gse131907_within_study_test.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    print("ELIGIBLE", results["eligible_by_origin"])
    for row in tested:
        print(
            row["program"],
            row["contrast"],
            "diff",
            round(row["mean_difference_a_minus_b"], 3),
            "p",
            f"{row['mannwhitney_p']:.4g}",
            "q",
            f"{row['fdr_among_six']:.4g}",
            "means",
            round(row["mean_a"], 3),
            round(row["mean_b"], 3),
            "n",
            row["n_a"],
            row["n_b"],
        )
    print("READING", json.dumps(results["reading"]))


if __name__ == "__main__":
    main()
