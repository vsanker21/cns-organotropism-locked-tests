"""Confirm the frozen mo-Mac genes in GSE123902 metastases.

Primary tumors are not used. No gene is added from this cohort.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(r"g:\Cancer Metastasis")
DATA = ROOT / "data" / "geo" / "GSE123902"
OUT = ROOT / "results" / "macrophage"
LOCKED = ROOT / "results" / "locked"

MAC = ["CD68", "CSF1R", "AIF1", "C1QA"]
OTHER = ["CD3D", "CD3E", "MS4A1", "CD79A", "EPCAM", "KRT8", "KRT18", "NKX2-1"]
GLIA_GATE = ["TMEM119", "P2RY12", "CX3CR1", "SALL1"]
GLIA_SENSITIVITY = [
    "GLDN", "SLC1A3", "GFAP", "PLP1", "CLDN11", "APOD",
    "MBP", "MOG", "OLIG1", "OLIG2", "AQP4", "S100B",
]
SAMPLES = [
    ("LX255B", "Brain", "GSM3516668_MSK_LX255B_METASTASIS_dense.csv.gz"),
    ("LX681", "Brain", "GSM3516671_MSK_LX681_METASTASIS_dense.csv.gz"),
    ("LX701", "Brain", "GSM3516678_MSK_LX701_METASTASIS_dense.csv.gz"),
    ("LX666", "Bone", "GSM3516664_MSK_LX666_METASTASIS_dense.csv.gz"),
    ("LX699", "Adrenal", "GSM3516677_MSK_LX699_METASTASIS_dense.csv.gz"),
]


def sum_genes(frame, genes):
    total = np.zeros(len(frame), dtype=float)
    for gene in genes:
        if gene in frame.columns:
            total += frame[gene].to_numpy(dtype=float)
    return total


def pseudobulk(path):
    frame = pd.read_csv(path, index_col=0)
    frame.columns = [str(column).upper() for column in frame.columns]
    mac = sum_genes(frame, MAC)
    other = sum_genes(frame, OTHER)
    glia = sum_genes(frame, GLIA_GATE)
    keep = (mac >= 2) & (mac > other) & (glia == 0)
    n_mac = int(keep.sum())
    if n_mac == 0:
        return None, {"n_cells": int(frame.shape[0]), "n_macrophage": 0, "umi": 0.0}
    cells = frame.loc[keep]
    library = float(cells.to_numpy(dtype=float).sum())
    totals = cells.sum(axis=0)
    return totals / library * 1e6, {
        "n_cells": int(frame.shape[0]),
        "n_macrophage": n_mac,
        "umi": library,
    }


def score_programs(log2cpm, programs):
    scores = {}
    used = {}
    for name, genes in programs.items():
        present = [gene for gene in genes if gene in log2cpm.columns and log2cpm[gene].notna().all()]
        values = log2cpm[present].astype(float)
        sd = values.std(axis=0, ddof=0)
        values = values.loc[:, sd > 0]
        z = (values - values.mean(axis=0)) / values.std(axis=0, ddof=0)
        scores[name] = z.mean(axis=1)
        used[name] = {"n_locked": len(genes), "n_used": int(z.shape[1])}
    return scores, used


def completely_separated(scores, sites, higher_group):
    brain = scores[sites == "Brain"].to_numpy()
    other = scores[sites != "Brain"].to_numpy()
    if higher_group == "Brain":
        return bool(len(brain) and len(other) and np.min(brain) > np.max(other))
    return bool(len(brain) and len(other) and np.max(brain) < np.min(other))


def test_one(scores, sites):
    brain = scores[sites == "Brain"].to_numpy()
    other = scores[sites != "Brain"].to_numpy()
    stat, p = mannwhitneyu(brain, other, alternative="two-sided")
    return {
        "n_brain": int(len(brain)),
        "n_other": int(len(other)),
        "mean_brain": float(brain.mean()),
        "mean_other": float(other.mean()),
        "min_brain": float(brain.min()),
        "max_brain": float(brain.max()),
        "min_other": float(other.min()),
        "max_other": float(other.max()),
        "mannwhitney_p": float(p),
        "mannwhitney_u": float(stat),
    }


def main():
    frozen = pd.read_csv(LOCKED / "momac_ambient_genes.csv")
    programs = {
        "brain_high": frozen.loc[frozen["direction"] == "brain_high", "gene"].str.upper().tolist(),
        "lymph_node_high": frozen.loc[frozen["direction"] == "lymph_node_high", "gene"].str.upper().tolist(),
    }
    reduced = {
        name: [gene for gene in genes if gene not in GLIA_SENSITIVITY]
        for name, genes in programs.items()
    }
    rows = []
    qc = []
    for sample_id, site, filename in SAMPLES:
        print(sample_id, site, flush=True)
        cpm, info = pseudobulk(DATA / filename)
        info.update({"sample": sample_id, "site": site})
        qc.append(info)
        print(f"  cells {info['n_cells']} macrophages {info['n_macrophage']}", flush=True)
        if cpm is None:
            continue
        record = {"sample": sample_id, "site": site}
        record.update(np.log2(cpm + 1).to_dict())
        rows.append(record)
    qc_df = pd.DataFrame(qc)
    qc_df.to_csv(OUT / "gse123902_macrophage_qc.csv", index=False)
    eligible_qc = qc_df[qc_df["n_macrophage"] >= 20]
    n_brain = int((eligible_qc["site"] == "Brain").sum())
    n_other = int((eligible_qc["site"] != "Brain").sum())
    result = {
        "qc": qc,
        "eligible_brain": n_brain,
        "eligible_other": n_other,
        "confirmed": False,
        "reason": "",
    }
    if n_brain < 2 or n_other < 2:
        result["reason"] = "fewer_than_two_eligible_samples_on_one_side"
        (OUT / "gse123902_confirmation.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
        print("STOP", result["reason"])
        return
    expr = pd.DataFrame(rows)
    expr = expr[expr["sample"].isin(eligible_qc["sample"])].set_index("sample")
    sites = expr["site"]
    expr = expr.drop(columns=["site"])
    primary_scores, primary_used = score_programs(expr, programs)
    reduced_scores, reduced_used = score_programs(expr, reduced)
    score_table = pd.DataFrame({"site": sites})
    comparisons = {}
    for label, score_map, usage in (
        ("primary", primary_scores, primary_used),
        ("without_glial_markers", reduced_scores, reduced_used),
    ):
        comparisons[label] = {}
        for name, higher in (("brain_high", "Brain"), ("lymph_node_high", "Other")):
            scores = score_map[name]
            score_table[f"{label}_{name}"] = scores
            entry = test_one(scores, sites)
            entry.update(usage[name])
            entry["completely_separated_in_predicted_direction"] = completely_separated(scores, sites, higher)
            comparisons[label][name] = entry
            print(label, name, json.dumps(entry))
    primary_ok = all(
        comparisons["primary"][name]["completely_separated_in_predicted_direction"]
        for name in ("brain_high", "lymph_node_high")
    )
    reduced_ok = all(
        comparisons["without_glial_markers"][name]["completely_separated_in_predicted_direction"]
        for name in ("brain_high", "lymph_node_high")
    )
    if primary_ok and reduced_ok:
        result["confirmed"] = True
        result["reason"] = "both_programs_separate_with_and_without_glial_markers"
    elif primary_ok and not reduced_ok:
        result["confirmed"] = False
        result["reason"] = "separation_depends_on_glial_markers"
    else:
        result["confirmed"] = False
        result["reason"] = "samples_overlap_or_direction_does_not_match"
    result["comparisons"] = comparisons
    score_table.to_csv(OUT / "gse123902_confirmation_scores.csv")
    (OUT / "gse123902_confirmation.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print("CONFIRMED", result["confirmed"], result["reason"])


if __name__ == "__main__":
    main()
