"""Orthogonal tests that do not change the locked Step A list.

Mouse models and the bone-tropic derivative were not used to choose genes.
The Smid brain-relapse set is a published comparator. On GSE12276 it is
circular, because that series is where the set was derived. On GSE76714 it
is an external check.
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

from importlib.machinery import SourceFileLoader

ROOT = Path(r"g:\Cancer Metastasis")
score_mod = SourceFileLoader(
    "score_cohorts",
    str(ROOT / "scripts" / "04_score_cohorts.py"),
).load_module()

OUT = ROOT / "results" / "scores"
LOCK = ROOT / "results" / "locked"


def classify_mouse(sample: str) -> str:
    if "BrM2" in sample:
        return "Parenchyma"
    if "_Par_" in sample or sample.endswith("_Par"):
        return "Parental"
    if "_IC_" in sample:
        return "Leptomeningeal"
    raise ValueError(sample)


def mouse_matrix() -> tuple[pd.DataFrame, pd.DataFrame]:
    folder = ROOT / "data" / "geo" / "GSE221358" / "suppl"
    frames = []
    for path in sorted(folder.glob("*.txt.gz")):
        model = path.name.split("_htseq_")[1].replace(".txt.gz", "")
        frame = pd.read_csv(path, sep="\t")
        frame["symbol"] = frame["GeneSymbol"].str.upper()
        frame = frame.groupby("symbol").sum(numeric_only=True)
        frame.columns = [f"{model}|{col}" for col in frame.columns]
        frames.append(frame)
    counts = pd.concat(frames, axis=1).fillna(0)
    library = counts.sum(axis=0)
    log_cpm = np.log2(counts.div(library, axis=1) * 1e6 + 1)
    meta_rows = []
    for col in log_cpm.columns:
        model, sample = col.split("|", 1)
        meta_rows.append({"sample_id": col, "model": model, "derivative": classify_mouse(sample)})
    return log_cpm, pd.DataFrame(meta_rows).set_index("sample_id")


def within_model_score(log_cpm: pd.DataFrame, meta: pd.DataFrame, genes: list[str]) -> pd.Series:
    present = [gene for gene in genes if gene in log_cpm.index]
    scores = []
    for model, rows in meta.groupby("model"):
        block = log_cpm.loc[present, rows.index]
        sd = block.std(axis=1, ddof=0)
        block = block.loc[sd > 0]
        z = block.sub(block.mean(axis=1), axis=0).div(block.std(axis=1, ddof=0), axis=0)
        scores.append(z.mean(axis=0))
    return pd.concat(scores)


def fit_mouse(scores: pd.Series, meta: pd.DataFrame) -> dict:
    frame = meta.copy()
    frame["score"] = scores.loc[frame.index]
    frame["derivative"] = pd.Categorical(frame["derivative"], ["Parental", "Parenchyma", "Leptomeningeal"])
    model = smf.ols("score ~ C(derivative, Treatment(reference='Parental')) + C(model)", data=frame).fit()
    means = frame.groupby("derivative")["score"].agg(["mean", "count"]).reset_index()
    return {
        "n": int(len(frame)),
        "n_genes_scored": None,
        "means": means.to_dict(orient="records"),
        "params": {key: float(value) for key, value in model.params.items()},
        "pvalues": {key: float(value) for key, value in model.pvalues.items()},
        "conf_int": {
            key: [float(low), float(high)]
            for key, (low, high) in model.conf_int().iterrows()
        },
    }


def bone_control(genes: dict[str, list[str]]) -> dict:
    path = ROOT / "data" / "geo" / "GSE112445" / "suppl" / "GSE112445_RNAseqAllTable-count.txt.gz"
    frame = pd.read_csv(path, sep="\t")
    frame["symbol"] = frame["official_symbol"].str.upper()
    keep = ["231-s1", "231-s2", "231-s3", "1833-s1", "1833-s2", "1833-s3"]
    counts = frame.groupby("symbol")[keep].sum()
    log_cpm = np.log2(counts.div(counts.sum(axis=0), axis=1) * 1e6 + 1)
    out = {}
    for name, gene_list in genes.items():
        present = [gene for gene in gene_list if gene in log_cpm.index]
        block = log_cpm.loc[present]
        sd = block.std(axis=1, ddof=0)
        block = block.loc[sd > 0]
        z = block.sub(block.mean(axis=1), axis=0).div(block.std(axis=1, ddof=0), axis=0)
        score = z.mean(axis=0)
        parental = score[["231-s1", "231-s2", "231-s3"]]
        bone = score[["1833-s1", "1833-s2", "1833-s3"]]
        _t, p = stats.ttest_ind(bone, parental, equal_var=False)
        out[name] = {
            "n_genes": len(gene_list),
            "n_present": int(block.shape[0]),
            "mean_parental": float(parental.mean()),
            "mean_bone_tropic": float(bone.mean()),
            "welch_p": float(p),
        }
    return out


def score_smid() -> dict:
    genes = score_mod.read_genes("smid_brain_up_genes.txt")
    results = {"n_genes": len(genes)}
    expr, meta = score_mod.prepare_expression("GSE76714")
    meta = meta.set_index("gsm")
    scores, coverage = score_mod.signature_scores(expr, genes)
    groups = meta.loc[scores.index, "chars"].map(
        lambda chars: score_mod.char_map(chars).get("patient group", "")
    )
    y = groups.map(
        lambda value: 1 if "with brain metastasis" in value.lower() else 0 if "without brain metastasis" in value.lower() else np.nan
    )
    results["GSE76714_external"] = score_mod.test_binary(y, scores, "GSE76714:smid_brain_up")
    results["GSE76714_external"]["coverage"] = coverage

    expr, meta = score_mod.prepare_expression("GSE12276")
    meta = meta.set_index("gsm")
    scores, coverage = score_mod.signature_scores(expr, genes)
    relapse = meta.loc[scores.index, "chars"].map(
        lambda chars: score_mod.char_map(chars).get("site of relapse (brain or other)", "")
    )
    y = relapse.map(lambda value: 1 if "brain" in value.lower() else 0 if value.lower() == "other" else np.nan)
    results["GSE12276_circular"] = score_mod.test_binary(y, scores, "GSE12276:smid_brain_up")
    results["GSE12276_circular"]["coverage"] = coverage
    results["GSE12276_circular"]["note"] = "This gene set was derived from GSE12276, so agreement here only checks the scoring pipeline."
    return results


def main() -> None:
    step_a = score_mod.read_genes("stepA_genes.txt")
    neural = score_mod.read_genes("neural_genes.txt")
    ortho = pd.read_csv(LOCK / "stepA_mouse_orthologs.csv")
    mouse_genes = ortho["mouse_symbol"].str.upper().tolist()
    log_cpm, meta = mouse_matrix()
    scores = within_model_score(log_cpm, meta, mouse_genes)
    mouse = fit_mouse(scores, meta)
    mouse["n_orthologs"] = len(set(mouse_genes))
    mouse["n_orthologs_in_matrix"] = int(sum(gene in log_cpm.index for gene in set(mouse_genes)))
    mouse["derivative_counts"] = meta["derivative"].value_counts().to_dict()
    bone = bone_control({"stepA": step_a, "neural": neural})
    smid = score_smid()
    payload = {"mouse_stepA": mouse, "bone_tropic_1833_vs_parental": bone, "smid_brain_up": smid}
    (OUT / "orthogonal.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
