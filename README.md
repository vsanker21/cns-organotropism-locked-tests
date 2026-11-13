# Locked tests of central nervous system organotropism

Vivek Sanker

Analysis code and the result tables for pre-specified scores of a culture program, a neural program, three mediators, a mouse leptomeningeal-versus-parenchymal program, and a macrophage contrast.

`results/locked/` holds the gene lists and the rules written before those scores. Scripts `01` through `22`, and `decontx_em.py`, score and test those rules. `60_posthoc_sensitivity.py` fits a Firth penalized logistic model to the stored culture and neural scores, a linear G2M adjustment to the stored parenchymal scores, and the same model with the median number of genes detected per malignant cell. `60_detection_depth.py` counts those medians. Both were run after the locked estimates and are not part of the lock. The median-gene tables are in `results/posthoc/`.

The expression matrices are the public GEO series named with the results. They are not included here.

https://github.com/vsanker21/cns-organotropism-locked-tests

https://doi.org/10.5281/zenodo.24197469

The git history contains the original deposit commit and a later commit that adds the detection-depth check. The original commit does not itself show that the locked files were written before the scores.
