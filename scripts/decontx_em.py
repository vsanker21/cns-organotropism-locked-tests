"""DecontX EM, matching decontX 1.8 (Yang et al., Genome Biology 2020).

theta is the native fraction of each cell. phi is the native gene distribution
of a cluster. eta is that cluster's contamination distribution, the native
expression of the other clusters. The Beta prior on theta is updated by maximum
likelihood, the same role as MCMCprecision::fit_dirichlet in the package.
"""

import numpy as np
from scipy import sparse
from scipy.special import gammaln
from scipy.optimize import minimize


def _beta_mle(native_proportion):
    p = np.clip(native_proportion, 1e-8, 1.0 - 1e-8)
    mean = float(p.mean())
    variance = float(p.var(ddof=1)) if len(p) > 1 else 0.0
    start = np.array([10.0, 10.0])
    if variance > 1e-8:
        common = mean * (1.0 - mean) / variance - 1.0
        if common > 0:
            start = np.array([max(mean * common, 1e-3), max((1.0 - mean) * common, 1e-3)])

    def objective(ab):
        a, b = ab
        return -(
            (a - 1.0) * np.log(p).sum()
            + (b - 1.0) * np.log(1.0 - p).sum()
            - len(p) * (gammaln(a) + gammaln(b) - gammaln(a + b))
        )

    fit = minimize(objective, start, method="L-BFGS-B", bounds=[(1e-3, 1e6), (1e-3, 1e6)])
    if not fit.success:
        return start
    return fit.x


def decontx_em(counts, labels, max_iter=100, convergence=0.001, delta=None, seed=12345):
    if delta is None:
        delta = np.array([10.0, 10.0])
    delta = np.asarray(delta, dtype=float)
    if sparse.issparse(counts):
        counts = counts.tocsr()
    labels = np.asarray(labels)
    n_genes, n_cells = counts.shape
    clusters = []
    for label in labels:
        if label not in clusters:
            clusters.append(label)
    index = {cluster: np.flatnonzero(labels == cluster) for cluster in clusters}
    if len(clusters) < 2:
        return {"contamination": np.zeros(n_cells), "iterations": 0, "native": counts, "delta": delta}

    rng = np.random.default_rng(seed)
    theta = rng.beta(delta[0], delta[1], size=n_cells)
    eps = 1e-20

    def cluster_block(ix):
        block = counts[:, ix]
        if sparse.issparse(block):
            block = block.toarray()
        return np.asarray(block, dtype=float)

    phi = np.full((n_genes, len(clusters)), eps)
    for k, cluster in enumerate(clusters):
        block = cluster_block(index[cluster])
        phi[:, k] += block @ theta[index[cluster]]
    eta = phi.sum(axis=1, keepdims=True) - phi
    phi = phi / phi.sum(axis=0, keepdims=True)
    eta = eta / eta.sum(axis=0, keepdims=True)
    library = np.asarray(counts.sum(axis=0)).ravel().astype(float)

    iterations = 0
    for iterations in range(1, max_iter + 1):
        previous = theta.copy()
        native_mass = np.zeros((n_genes, len(clusters)))
        native_total = np.zeros(n_cells)
        for k, cluster in enumerate(clusters):
            ix = index[cluster]
            block = cluster_block(ix)
            native_weight = (phi[:, k] + eps)[:, None] * (theta[ix] + eps)[None, :]
            contam_weight = (eta[:, k] + eps)[:, None] * (1.0 - theta[ix] + eps)[None, :]
            native_probability = native_weight / (native_weight + contam_weight)
            native = block * native_probability
            native_mass[:, k] = native.sum(axis=1)
            native_total[ix] = native.sum(axis=0)
        eta = native_mass.sum(axis=1, keepdims=True) - native_mass
        phi_sums = native_mass.sum(axis=0)
        eta_sums = eta.sum(axis=0)
        phi = np.divide(native_mass, phi_sums, out=np.full_like(native_mass, 1.0 / n_genes), where=phi_sums > 0)
        eta = np.divide(eta, eta_sums, out=np.full_like(eta, 1.0 / n_genes), where=eta_sums > 0)
        positive = library > 0
        native_proportion = np.zeros(n_cells)
        native_proportion[positive] = native_total[positive] / library[positive]
        delta = _beta_mle(native_proportion[positive])
        theta = np.zeros(n_cells)
        theta[positive] = (native_total[positive] + delta[0]) / (library[positive] + delta.sum())
        if np.max(np.abs(theta - previous)) < convergence:
            break

    native_counts = np.zeros((n_genes, n_cells), dtype=float)
    for k, cluster in enumerate(clusters):
        ix = index[cluster]
        block = cluster_block(ix)
        native_weight = (phi[:, k] + eps)[:, None] * (theta[ix] + eps)[None, :]
        contam_weight = (eta[:, k] + eps)[:, None] * (1.0 - theta[ix] + eps)[None, :]
        native_probability = native_weight / (native_weight + contam_weight)
        native_counts[:, ix] = block * native_probability
    contamination = np.zeros(n_cells)
    contamination[library > 0] = 1.0 - native_total[library > 0] / library[library > 0]
    # The last native_total is from the final E-step, before the last theta update.
    # Recompute contamination from the final theta, which is what decontX stores
    # only after another E-step. Use 1 - theta, the fitted native fraction.
    contamination = np.clip(1.0 - theta, 0.0, 1.0)
    return {
        "contamination": contamination,
        "iterations": iterations,
        "native": native_counts,
        "delta": delta,
        "theta": theta,
    }


def _simulate():
    rng = np.random.default_rng(20261112)
    n_genes = 100
    native = np.zeros(n_genes)
    native[:10] = 0.1
    foreign = np.zeros(n_genes)
    foreign[10:20] = 0.1
    n_cells = 60
    library = 2000
    beta_true = 0.25
    mixture = (1 - beta_true) * native + beta_true * foreign
    contaminated = np.vstack([rng.multinomial(library, mixture) for _ in range(n_cells)]).T
    clean = np.vstack([rng.multinomial(library, native) for _ in range(n_cells)]).T
    other = np.vstack([rng.multinomial(library, foreign) for _ in range(n_cells)]).T
    counts = np.concatenate([contaminated, clean, other], axis=1)
    labels = np.array(["contam"] * n_cells + ["clean"] * n_cells + ["source"] * n_cells)
    fit = decontx_em(counts, labels)
    recovered = float(np.median(fit["contamination"][:n_cells]))
    clean_level = float(np.median(fit["contamination"][n_cells : 2 * n_cells]))
    foreign_left = fit["native"][10:20, :n_cells].sum() / counts[10:20, :n_cells].sum()
    kept_native = fit["native"][:10, :n_cells].sum() / counts[:10, :n_cells].sum()
    print(
        f"SIM_CONTAM {recovered:.3f} CLEAN_CONTAM {clean_level:.3f} "
        f"FOREIGN_LEFT {foreign_left:.3f} NATIVE_KEPT {kept_native:.3f} ITER {fit['iterations']}"
    )
    if not (0.10 <= recovered <= 0.40 and clean_level < 0.15 and foreign_left < 0.55 and kept_native > 0.7):
        raise SystemExit("DecontX EM failed the simulation check")


if __name__ == "__main__":
    _simulate()
