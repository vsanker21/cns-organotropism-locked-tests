"""Patient-level CNS associations for a pre-specified mutation panel.

One row per patient. A patient is CNS-positive if any sample is marked
DMETS_DX_CNS_BRAIN = Yes. A gene is mutant if any of that patient's samples
has a nonsilent variant. Tests are Fisher exact tests, BH-adjusted within
each cancer type. This repeats a question already mapped by Nguyen et al.
and is reported as context, not as a new genomic screen.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

ROOT = Path(r"g:\Cancer Metastasis")
CLIN = ROOT / "data" / "cbioportal" / "msk_met_2021" / "data_clinical_sample.txt"
MUT = ROOT / "data" / "cbioportal" / "msk_met_2021" / "data_mutations.txt"
STEP_A = ROOT / "results" / "locked" / "stepA_genes.txt"
OUT = ROOT / "results" / "scores" / "msk_cns_mutations.json"

NONSILENT = {
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "Splice_Site",
    "In_Frame_Del",
    "In_Frame_Ins",
    "Translation_Start_Site",
    "Nonstop_Mutation",
}
PANEL = [
    "TP53", "PIK3CA", "ERBB2", "EGFR", "KRAS", "STK11", "KEAP1",
    "CDKN2A", "PTEN", "MYC", "CDH1", "BRCA1", "BRCA2", "RB1",
]


def bh(pvalues: list[float]) -> list[float]:
    arr = np.asarray(pvalues, dtype=float)
    order = np.argsort(arr)
    ranked = arr[order]
    m = len(ranked)
    adjusted = np.minimum.accumulate((ranked * m / np.arange(1, m + 1))[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.clip(adjusted, 0, 1)
    return out.tolist()


def fisher_table(cns: pd.Series, mutant: pd.Series) -> dict:
    a = int(((cns == 1) & (mutant == 1)).sum())
    b = int(((cns == 0) & (mutant == 1)).sum())
    c = int(((cns == 1) & (mutant == 0)).sum())
    d = int(((cns == 0) & (mutant == 0)).sum())
    odds, p = fisher_exact([[a, b], [c, d]])
    return {
        "cns_mutant": a,
        "no_cns_mutant": b,
        "cns_wildtype": c,
        "no_cns_wildtype": d,
        "odds_ratio": None if not np.isfinite(odds) else float(odds),
        "fisher_p": float(p),
    }


def main() -> None:
    clinical = pd.read_csv(CLIN, sep="\t", comment="#", low_memory=False)
    clinical = clinical.rename(columns={clinical.columns[0]: "SAMPLE_ID"})
    patient = clinical.groupby("PATIENT_ID").agg(
        cancer=("CANCER_TYPE", "first"),
        oncotree=("ONCOTREE_CODE", "first"),
        cns=("DMETS_DX_CNS_BRAIN", lambda values: int((values == "Yes").any())),
    )
    usecols = ["Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode"]
    mutations = pd.read_csv(MUT, sep="\t", usecols=usecols, low_memory=False)
    mutations = mutations.loc[mutations["Variant_Classification"].isin(NONSILENT)]
    sample_to_patient = clinical.set_index("SAMPLE_ID")["PATIENT_ID"]
    mutations["PATIENT_ID"] = mutations["Tumor_Sample_Barcode"].map(sample_to_patient)
    mutations = mutations.dropna(subset=["PATIENT_ID"])
    mutant_pairs = mutations.groupby(["PATIENT_ID", "Hugo_Symbol"]).size().rename("n").reset_index()

    step_a = {line.strip().upper() for line in STEP_A.read_text(encoding="utf-8").splitlines() if line.strip()}
    cohorts = {
        "breast": patient.loc[patient["cancer"] == "Breast Cancer"],
        "lung_adenocarcinoma": patient.loc[patient["oncotree"] == "LUAD"],
    }
    results = {}
    for name, cohort in cohorts.items():
        genes = list(PANEL)
        counts = mutant_pairs.loc[mutant_pairs["PATIENT_ID"].isin(cohort.index), "Hugo_Symbol"].value_counts()
        extra = [gene for gene in sorted(step_a) if counts.get(gene, 0) >= 20 and gene not in genes]
        rows = []
        for gene in genes + extra:
            mutant_patients = set(mutant_pairs.loc[mutant_pairs["Hugo_Symbol"] == gene, "PATIENT_ID"])
            mutant = cohort.index.map(lambda patient_id: int(patient_id in mutant_patients))
            row = fisher_table(cohort["cns"], pd.Series(mutant, index=cohort.index))
            row["gene"] = gene
            row["in_stepA"] = gene in step_a
            row["pre_specified_panel"] = gene in PANEL
            rows.append(row)
        qvals = bh([row["fisher_p"] for row in rows])
        for row, q in zip(rows, qvals):
            row["q_bh"] = q
        results[name] = {
            "n_patients": int(len(cohort)),
            "n_cns": int(cohort["cns"].sum()),
            "tests": rows,
        }
    OUT.write_text(json.dumps(results, indent=2), encoding="utf-8")
    summary = {
        name: {
            "n_patients": payload["n_patients"],
            "n_cns": payload["n_cns"],
            "significant": [
                {key: row[key] for key in ("gene", "odds_ratio", "fisher_p", "q_bh", "cns_mutant", "no_cns_mutant")}
                for row in payload["tests"]
                if row["q_bh"] < 0.05
            ],
        }
        for name, payload in results.items()
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
