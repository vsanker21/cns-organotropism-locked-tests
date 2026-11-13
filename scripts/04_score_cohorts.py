"""Score locked signatures on held-out cohorts.

Step A genes and the benchmark sets are read from results/locked.
No gene is added or removed based on these cohorts.
"""
from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from sklearn.metrics import roc_auc_score

ROOT = Path(r"g:\Cancer Metastasis")
DATA = ROOT / "data" / "geo"
LOCK = ROOT / "results" / "locked"
OUT = ROOT / "results" / "scores"
OUT.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(20261112)


def read_genes(name: str) -> list[str]:
    genes = []
    for line in (LOCK / name).read_text(encoding="utf-8").splitlines():
        gene = line.strip().upper()
        if gene:
            genes.append(gene)
    return list(dict.fromkeys(genes))


def parse_soft(acc: str) -> tuple[pd.DataFrame, dict[str, str]]:
    path = DATA / acc / "soft" / f"{acc}_family.soft.gz"
    samples: list[dict] = []
    current: dict | None = None
    platform_map: dict[str, str] = {}
    in_platform = False
    platform_cols: list[str] | None = None
    id_idx = symbol_idx = None

    def flush() -> None:
        nonlocal current
        if current is not None:
            samples.append(current)
            current = None

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("^SAMPLE"):
                flush()
                current = {"gsm": line.split("=", 1)[1].strip(), "title": "", "chars": []}
            elif current is not None and line.startswith("!Sample_title"):
                current["title"] = line.split("=", 1)[1].strip()
            elif current is not None and line.startswith("!Sample_characteristics"):
                current["chars"].append(line.split("=", 1)[1].strip())
            elif current is not None and line.startswith("!Sample_source_name"):
                current["chars"].append(line.split("=", 1)[1].strip())
            elif line.startswith("!platform_table_begin"):
                flush()
                in_platform = True
                platform_cols = None
            elif line.startswith("!platform_table_end"):
                in_platform = False
            elif in_platform:
                parts = line.rstrip("\n").split("\t")
                if platform_cols is None:
                    platform_cols = parts
                    lower = [col.strip().strip('"').lower() for col in platform_cols]
                    id_idx = 0
                    symbol_idx = next(
                        (i for i, col in enumerate(lower) if col in {"gene symbol", "symbol", "ilmn_gene", "gene_symbol"}),
                        None,
                    )
                    if symbol_idx is None:
                        symbol_idx = next((i for i, col in enumerate(lower) if "symbol" in col), None)
                elif symbol_idx is not None and len(parts) > symbol_idx:
                    probe = parts[id_idx].strip().strip('"')
                    symbol = parts[symbol_idx].strip().strip('"').split("///")[0].strip().upper()
                    if probe and symbol and symbol not in {"NA", "---", "NULL"}:
                        platform_map.setdefault(probe, symbol)
    flush()
    meta = pd.DataFrame(samples)
    return meta, platform_map


def read_series_matrix(acc: str) -> pd.DataFrame:
    path = next((DATA / acc / "matrix").glob("*series_matrix.txt.gz"))
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
        frame = pd.read_csv(handle, sep="\t", quotechar='"', index_col=0)
    frame = frame.loc[~frame.index.astype(str).str.startswith("!")]
    frame = frame.apply(pd.to_numeric, errors="coerce")
    return frame


def prepare_expression(acc: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta, platform_map = parse_soft(acc)
    expr = read_series_matrix(acc)
    values = expr.to_numpy(dtype=float)
    if np.nanmin(values) >= 0 and np.nanpercentile(values, 99) > 50:
        expr = np.log2(expr.clip(lower=0) + 1)
    symbols = pd.Series({probe: platform_map.get(str(probe), "") for probe in expr.index})
    expr = expr.loc[symbols.ne("")]
    symbols = symbols.loc[symbols.ne("")]
    # One probe per gene: the probe with the highest mean log expression.
    means = expr.mean(axis=1)
    chosen_rows = []
    chosen_genes = []
    for gene, probes in symbols.groupby(symbols).groups.items():
        probe_means = means.loc[list(probes)]
        chosen_rows.append(probe_means.idxmax())
        chosen_genes.append(gene)
    collapsed = expr.loc[chosen_rows]
    collapsed.index = chosen_genes
    collapsed = collapsed.groupby(level=0).mean()
    return collapsed, meta


def signature_scores(expr: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, dict]:
    present = [gene for gene in genes if gene in expr.index]
    coverage = len(present) / len(genes) if genes else 0
    if len(present) < 15:
        scores = pd.Series(np.nan, index=expr.columns, name="score")
        return scores, {"n_present": len(present), "n_genes": len(genes), "coverage": coverage, "valid": False}
    sub = expr.loc[present]
    sd = sub.std(axis=1, ddof=0)
    sub = sub.loc[sd > 0]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1, ddof=0), axis=0)
    scores = z.mean(axis=0)
    return scores, {
        "n_present": int(sub.shape[0]),
        "n_genes": len(genes),
        "coverage": coverage,
        "valid": True,
    }


def bootstrap_auc(y: np.ndarray, score: np.ndarray) -> tuple[float, float, float]:
    point = float(roc_auc_score(y, score))
    n = len(y)
    draws = []
    for _ in range(2000):
        idx = RNG.integers(0, n, n)
        if len(np.unique(y[idx])) < 2:
            continue
        draws.append(roc_auc_score(y[idx], score[idx]))
    low, high = np.percentile(draws, [2.5, 97.5])
    return point, float(low), float(high)


def test_binary(y: pd.Series, score: pd.Series, label: str) -> dict:
    frame = pd.DataFrame({"y": y, "score": score}).dropna()
    frame = frame.loc[frame["y"].isin([0, 1])]
    out = {
        "label": label,
        "n": int(len(frame)),
        "n_positive": int(frame["y"].sum()) if len(frame) else 0,
    }
    if frame["y"].nunique() < 2 or len(frame) < 10:
        out["error"] = "not_enough_outcome_variation"
        return out
    y_arr = frame["y"].to_numpy(dtype=float)
    x = frame["score"].to_numpy(dtype=float)
    x_std = (x - x.mean()) / x.std(ddof=0)
    pos = x[y_arr == 1]
    neg = x[y_arr == 0]
    u_stat, u_p = stats.mannwhitneyu(pos, neg, alternative="two-sided")
    auc, auc_low, auc_high = bootstrap_auc(y_arr, x)
    model = sm.Logit(y_arr, sm.add_constant(x_std)).fit(disp=False)
    conf = np.asarray(model.conf_int())
    log_or = float(np.asarray(model.params)[1])
    log_p = float(np.asarray(model.pvalues)[1])
    out.update(
        {
            "mean_score_positive": float(pos.mean()),
            "mean_score_negative": float(neg.mean()),
            "mannwhitney_p": float(u_p),
            "auc": auc,
            "auc_ci95": [auc_low, auc_high],
            "log_or_per_sd": log_or,
            "or_per_sd": float(np.exp(log_or)),
            "or_ci95": [float(np.exp(conf[1, 0])), float(np.exp(conf[1, 1]))],
            "logit_p": log_p,
        }
    )
    return out


def bh(pvalues: list[float]) -> list[float]:
    arr = np.array(pvalues, dtype=float)
    out = np.full(len(arr), np.nan)
    valid = np.isfinite(arr)
    if valid.sum() == 0:
        return out.tolist()
    order = np.argsort(arr[valid])
    ranked = arr[valid][order]
    m = len(ranked)
    adjusted = ranked * m / (np.arange(1, m + 1))
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    placed = np.empty(m)
    placed[order] = adjusted
    out[np.where(valid)[0]] = placed
    return out.tolist()


def char_map(chars: list[str]) -> dict[str, str]:
    mapped = {}
    for item in chars:
        if ":" in item:
            key, value = item.split(":", 1)
            mapped[key.strip().lower()] = value.strip()
    return mapped


def organ_from_title(title: str) -> str | None:
    match = re.search(r"(Brain|Lung|Bone|Liver) metastasis", title, flags=re.I)
    return match.group(1).capitalize() if match else None


def score_cohort(acc: str, expr: pd.DataFrame, sets: dict[str, list[str]]) -> tuple[pd.DataFrame, dict]:
    columns = {}
    coverage = {}
    for name, genes in sets.items():
        scores, info = signature_scores(expr, genes)
        columns[name] = scores
        coverage[name] = info
    frame = pd.DataFrame(columns)
    frame.index.name = "gsm"
    return frame, coverage


def main() -> None:
    sets = {
        "stepA": read_genes("stepA_genes.txt"),
        "neural": read_genes("neural_genes.txt"),
        "bos": read_genes("bos_genes.txt"),
        "proliferation": read_genes("proliferation_genes.txt"),
    }
    report: dict = {"signatures": {name: {"n_genes": len(genes)} for name, genes in sets.items()}}
    report["stepA_rule"] = (LOCK / "stepA_rule.txt").read_text(encoding="utf-8").strip()
    report["neural_set"] = (LOCK / "neural_set_name.txt").read_text(encoding="utf-8").strip()
    report["bos_set"] = (LOCK / "bos_set_name.txt").read_text(encoding="utf-8").strip()

    # Primary cohort: TNBC primaries, brain metastasis versus none.
    expr, meta = prepare_expression("GSE76714")
    meta = meta.set_index("gsm")
    scores, coverage = score_cohort("GSE76714", expr, sets)
    meta = meta.loc[scores.index]
    groups = meta["chars"].map(lambda chars: char_map(chars).get("patient group", ""))
    y = groups.map(lambda value: 1 if "with brain metastasis" in value.lower() else 0 if "without brain metastasis" in value.lower() else np.nan)
    primary_tests = []
    for name in sets:
        result = test_binary(y, scores[name], f"GSE76714:{name}")
        result["coverage"] = coverage[name]
        primary_tests.append(result)
    adjusted = bh([item.get("logit_p", np.nan) for item in primary_tests])
    for item, q in zip(primary_tests, adjusted):
        item["logit_q_bh"] = q
    saved = scores.copy()
    saved["brain_met"] = y.to_numpy()
    saved.to_csv(OUT / "gse76714_scores.csv")
    report["GSE76714"] = {
        "n_with_brain_met": int(np.nansum(y.to_numpy() == 1)),
        "n_without_brain_met": int(np.nansum(y.to_numpy() == 0)),
        "tests": primary_tests,
    }

    # Replication cohort. The Bos genes were historically tied to this series.
    expr, meta = prepare_expression("GSE12276")
    meta = meta.set_index("gsm")
    scores, coverage = score_cohort("GSE12276", expr, sets)
    meta = meta.loc[scores.index]
    relapse = meta["chars"].map(lambda chars: char_map(chars).get("site of relapse (brain or other)", ""))
    y = relapse.map(lambda value: 1 if "brain" in value.lower() else 0 if value.lower() == "other" else np.nan)
    tests = []
    for name in sets:
        result = test_binary(y, scores[name], f"GSE12276:{name}")
        result["coverage"] = coverage[name]
        result["bos_is_circular_on_this_cohort"] = name == "bos"
        tests.append(result)
    saved = scores.copy()
    saved["brain_relapse"] = y.to_numpy()
    saved["relapse_label"] = relapse.to_numpy()
    saved.to_csv(OUT / "gse12276_scores.csv")
    report["GSE12276"] = {
        "label_counts": relapse.value_counts(dropna=False).to_dict(),
        "tests": tests,
    }

    # Established metastases by organ. Neural scores here include glial RNA.
    organ_frames = []
    for acc in ("GSE14017", "GSE14018"):
        expr, meta = prepare_expression(acc)
        meta = meta.set_index("gsm")
        scores, coverage = score_cohort(acc, expr, sets)
        meta = meta.loc[scores.index]
        organ = meta["title"].map(organ_from_title)
        frame = scores.copy()
        frame["organ"] = organ.to_numpy()
        frame["acc"] = acc
        frame["tumor"] = meta["title"].str.extract(r"^([^:]+)", expand=False).str.strip().to_numpy()
        organ_frames.append((frame, coverage))
    organs = pd.concat([frame for frame, _ in organ_frames], axis=0)
    organs = organs.dropna(subset=["organ"])
    # Prefer the U133 Plus 2 array when the same tumor was profiled twice.
    organs["platform_rank"] = organs["acc"].map({"GSE14017": 0, "GSE14018": 1})
    organs = organs.sort_values("platform_rank").drop_duplicates("tumor", keep="first")
    organ_tests = {}
    for name in ("stepA", "neural"):
        grouped = [organs.loc[organs["organ"] == organ, name].dropna() for organ in ["Brain", "Lung", "Bone", "Liver"]]
        if min(len(item) for item in grouped) < 3:
            continue
        h_stat, h_p = stats.kruskal(*grouped)
        pairs = []
        raw_p = []
        for other in ["Lung", "Bone", "Liver"]:
            a = organs.loc[organs["organ"] == "Brain", name].dropna()
            b = organs.loc[organs["organ"] == other, name].dropna()
            _u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            pairs.append({"contrast": f"Brain_vs_{other}", "p": float(p), "median_brain": float(a.median()), "median_other": float(b.median())})
            raw_p.append(p)
        for pair, q in zip(pairs, bh(raw_p)):
            pair["q_bh"] = q
        organ_tests[name] = {
            "kruskal_p": float(h_p),
            "kruskal_h": float(h_stat),
            "pairwise": pairs,
            "n": organs.groupby("organ")[name].size().to_dict(),
        }
    organs.to_csv(OUT / "metastasis_site_scores.csv", index=False)
    report["metastasis_sites"] = organ_tests

    # Independent cell-line check: CN34 was not used to learn Step A.
    # MDA-MB-231 on this array is the same model system as the discovery RNA-seq.
    expr, meta = prepare_expression("GSE12237")
    meta = meta.set_index("gsm")
    scores, coverage = score_cohort("GSE12237", expr, {"stepA": sets["stepA"], "neural": sets["neural"]})
    meta = meta.loc[scores.index]
    titles = meta["title"]
    records = []
    for background in ("MDA-MB-231", "CN34"):
        sub_title = titles[titles.str.startswith(background)]
        sub = scores.loc[sub_title.index].copy()
        sub["class"] = np.where(sub_title.str.contains("Parental"), "Parental", "BrM")
        for name in ("stepA", "neural"):
            a = sub.loc[sub["class"] == "BrM", name]
            b = sub.loc[sub["class"] == "Parental", name]
            p = float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue) if len(a) and len(b) else np.nan
            records.append({
                "background": background,
                "signature": name,
                "n_brm": int(len(a)),
                "n_parental": int(len(b)),
                "mean_brm": float(a.mean()),
                "mean_parental": float(b.mean()),
                "mannwhitney_p": p,
                "independent_of_discovery": background == "CN34",
            })
    report["GSE12237"] = {"coverage": coverage, "contrasts": records}

    (OUT / "adjudication.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({
        "GSE76714": report["GSE76714"],
        "GSE12276_labels": report["GSE12276"]["label_counts"],
        "GSE12276_tests": report["GSE12276"]["tests"],
        "sites": {k: {"kruskal_p": v["kruskal_p"], "n": v["n"], "pairwise": v["pairwise"]} for k, v in organ_tests.items()},
        "GSE12237": records,
    }, indent=2))


if __name__ == "__main__":
    main()
