# Locked tests of central nervous system organotropism

Vivek Sanker

Analysis code and the result tables for pre-specified scores of a culture program, a neural program, three mediators, a mouse leptomeningeal-versus-parenchymal program, and a macrophage contrast.

`results/locked/` holds the gene lists and the rules written before those scores. Scripts `01` through `22`, and `decontx_em.py`, score and test those rules. `60_posthoc_sensitivity.py` fits a Firth penalized logistic model to the stored culture and neural scores and a linear G2M adjustment to the stored parenchymal scores. That script was run after the locked estimates and is not part of the lock.

The expression matrices are the public GEO series named with the results. They are not included here.

https://github.com/vsanker21/cns-organotropism-locked-tests

This archive records the files at the time of deposit. One commit made for the deposit does not itself show that the locked files were written before the scores.
