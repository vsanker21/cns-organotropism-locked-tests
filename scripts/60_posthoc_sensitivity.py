"""Post hoc checks on stored scores. These checks were not part of the lock."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
SEED = 20261112


def sigmoid(eta: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(eta, -30.0, 30.0)))


def penalized_loglik(beta: np.ndarray, x: np.ndarray, y: np.ndarray) -> float:
    p = np.clip(sigmoid(x @ beta), 1e-12, 1.0 - 1e-12)
    w = p * (1.0 - p)
    ll = float(np.sum(y * np.log(p) + (1.0 - y) * np.log(1.0 - p)))
    sign, logdet = np.linalg.slogdet(x.T @ (x * w[:, None]))
    if sign <= 0:
        return -np.inf
    return ll + 0.5 * logdet


def firth_fit(x: np.ndarray, y: np.ndarray, max_iter: int = 80) -> np.ndarray:
    beta = np.zeros(x.shape[1])
    current = penalized_loglik(beta, x, y)
    for _ in range(max_iter):
        p = sigmoid(x @ beta)
        w = np.clip(p * (1.0 - p), 1e-12, None)
        fisher = x.T @ (x * w[:, None])
        cov = np.linalg.pinv(fisher)
        leverage = np.sum((x @ cov) * (x * w[:, None]), axis=1)
        score = x.T @ (y - p + leverage * (0.5 - p))
        step = cov @ score
        accepted = beta
        accepted_ll = current
        for shrink in (1.0, 0.5, 0.25, 0.1, 0.01):
            trial = beta + shrink * step
            trial_ll = penalized_loglik(trial, x, y)
            if trial_ll >= current - 1e-10:
                accepted = trial
                accepted_ll = trial_ll
                break
        if np.max(np.abs(accepted - beta)) < 1e-12:
            return accepted
        beta = accepted
        current = accepted_ll
    return beta


def profile_limit(x: np.ndarray, y: np.ndarray, mle: np.ndarray, target: float, upper: bool) -> float:
    """Profile-likelihood bound for the slope. target is the chi-square cutoff."""
    maximum = penalized_loglik(mle, x, y)

    def gap(slope: float) -> float:
        def objective(intercept: float) -> float:
            return -penalized_loglik(np.array([intercept, slope]), x, y)

        best = None
        for start in (mle[0], 0.0, -2.0, 2.0, -4.0, 4.0):
            grid = np.linspace(start - 3.0, start + 3.0, 25)
            pick = grid[np.argmin([objective(value) for value in grid])]
            # Newton on the intercept at a fixed slope.
            intercept = pick
            for _ in range(40):
                beta = np.array([intercept, slope])
                p = sigmoid(x @ beta)
                w = np.clip(p * (1.0 - p), 1e-12, None)
                fisher = x.T @ (x * w[:, None])
                cov = np.linalg.pinv(fisher)
                leverage = np.sum((x @ cov) * (x * w[:, None]), axis=1)
                score = x.T @ (y - p + leverage * (0.5 - p))
                intercept = intercept + float(cov[0] @ score)
            value = objective(intercept)
            if best is None or value < best[0]:
                best = (value, intercept)
        return 2.0 * (maximum + best[0]) - target

    center = float(mle[1])
    if upper:
        lo, hi = center, center + 0.5
        while gap(hi) < 0 and hi < center + 20:
            hi += 0.5
    else:
        hi, lo = center, center - 0.5
        while gap(lo) < 0 and lo > center - 20:
            lo -= 0.5
    return float(brentq(gap, lo, hi))


def firth_or(y: np.ndarray, score: np.ndarray) -> dict:
    x_std = (score - score.mean()) / score.std(ddof=0)
    design = sm.add_constant(x_std)
    wald = sm.Logit(y, design).fit(disp=False)
    wald_ci = np.asarray(wald.conf_int())
    mle = firth_fit(design, y)
    cutoff = float(stats.chi2.ppf(0.95, 1))
    low = profile_limit(design, y, mle, cutoff, upper=False)
    high = profile_limit(design, y, mle, cutoff, upper=True)
    return {
        "wald_or": float(np.exp(wald.params[1])),
        "wald_low": float(np.exp(wald_ci[1, 0])),
        "wald_high": float(np.exp(wald_ci[1, 1])),
        "firth_or": float(np.exp(mle[1])),
        "firth_low": float(np.exp(low)),
        "firth_high": float(np.exp(high)),
        "n": int(len(y)),
        "events": int(y.sum()),
    }


def load_relapse() -> list[tuple[str, str, np.ndarray, np.ndarray]]:
    rows = []
    specs = (
        ("GSE76714", "stepA", "brain_met", "results/scores/gse76714_scores.csv"),
        ("GSE76714", "neural", "brain_met", "results/scores/gse76714_scores.csv"),
        ("GSE12276", "stepA", "brain_relapse", "results/scores/gse12276_scores.csv"),
        ("GSE12276", "neural", "brain_relapse", "results/scores/gse12276_scores.csv"),
    )
    for cohort, column, outcome, path in specs:
        frame = pd.read_csv(ROOT / path)
        frame = frame.loc[frame[outcome].isin([0, 1, 0.0, 1.0])].copy()
        y = frame[outcome].to_numpy(dtype=float)
        score = frame[column].to_numpy(dtype=float)
        rows.append((cohort, column, y, score))
    return rows


def g2m_adjustment() -> dict:
    frame = pd.read_csv(ROOT / "results/compartment/human_sample_scores.csv")
    fluid = (frame["compartment"] == "CSF").to_numpy(dtype=float)
    par = frame["parenchyma"].to_numpy(dtype=float)
    g2m = frame["proliferation_control"].to_numpy(dtype=float)
    design = np.column_stack([np.ones(len(frame)), g2m, fluid])
    model = sm.OLS(par, design).fit()
    rng = np.random.default_rng(SEED)
    observed = float(model.params[2])
    count = 0
    for _ in range(10000):
        shuffled = rng.permutation(fluid)
        fit = sm.OLS(par, np.column_stack([np.ones(len(frame)), g2m, shuffled])).fit()
        if abs(float(fit.params[2])) >= abs(observed):
            count += 1
    return {
        "n": int(len(frame)),
        "n_fluid": int(fluid.sum()),
        "coefficient": observed,
        "ols_p": float(model.pvalues[2]),
        "permutation_p": (1 + count) / 10001,
        "mean_fluid": float(par[fluid == 1].mean()),
        "mean_parenchyma": float(par[fluid == 0].mean()),
        "g2m_coefficient": float(model.params[1]),
    }


def depth_adjustment() -> dict:
    """Parenchymal score on G2M, median genes detected, and the compartment label.

    Median genes are the per-sample median, across malignant cells, of the
    number of symbols with a raw count above zero. The compartment label is
    permuted. The gene medians and the G2M scores stay with the sample.
    """
    scores = pd.read_csv(ROOT / "results/compartment/human_sample_scores.csv")
    depth = pd.concat(
        [
            pd.read_csv(ROOT / "results/posthoc/csf_median_genes.csv"),
            pd.read_csv(ROOT / "results/posthoc/parenchyma_median_genes.csv"),
        ],
        ignore_index=True,
    )
    frame = scores.merge(depth, on=["sample", "compartment", "n_malignant"], how="inner")
    if len(frame) != 16:
        raise RuntimeError(f"depth model has {len(frame)} samples, expected 16")
    fluid = (frame["compartment"] == "CSF").to_numpy(dtype=float)
    par = frame["parenchyma"].to_numpy(dtype=float)
    g2m = frame["proliferation_control"].to_numpy(dtype=float)
    genes = frame["median_genes"].to_numpy(dtype=float)
    design = np.column_stack([np.ones(len(frame)), g2m, genes, fluid])
    model = sm.OLS(par, design).fit()
    rng = np.random.default_rng(SEED)
    observed = float(model.params[3])
    count = 0
    for _ in range(10000):
        shuffled = rng.permutation(fluid)
        fit = sm.OLS(par, np.column_stack([np.ones(len(frame)), g2m, genes, shuffled])).fit()
        if abs(float(fit.params[3])) >= abs(observed):
            count += 1
    fluid_genes = genes[fluid == 1]
    solid_genes = genes[fluid == 0]
    return {
        "n": int(len(frame)),
        "coefficient": observed,
        "g2m_coefficient": float(model.params[1]),
        "genes_coefficient": float(model.params[2]),
        "ols_p": float(model.pvalues[3]),
        "permutation_p": (1 + count) / 10001,
        "permutation_count": int(count),
        "fluid_median_genes_min": float(fluid_genes.min()),
        "fluid_median_genes_max": float(fluid_genes.max()),
        "parenchyma_median_genes_min": float(solid_genes.min()),
        "parenchyma_median_genes_max": float(solid_genes.max()),
        "ranges_overlap": bool(fluid_genes.min() <= solid_genes.max() and solid_genes.min() <= fluid_genes.max()),
    }


def main() -> None:
    for cohort, column, y, score in load_relapse():
        out = firth_or(y, score)
        print(
            cohort,
            column,
            f"wald {out['wald_or']:.4f} ({out['wald_low']:.4f}-{out['wald_high']:.4f})",
            f"firth {out['firth_or']:.4f} ({out['firth_low']:.4f}-{out['firth_high']:.4f})",
            f"n {out['n']} events {out['events']}",
        )
    print(g2m_adjustment())
    depth_path = ROOT / "results/posthoc/parenchyma_median_genes.csv"
    if depth_path.exists() and (ROOT / "results/posthoc/csf_median_genes.csv").exists():
        print(depth_adjustment())


if __name__ == "__main__":
    main()
