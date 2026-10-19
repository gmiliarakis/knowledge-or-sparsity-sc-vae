# Review tracker

Living list of every weakness raised by the referee-style review (2026-10-14) and what was done about it. Goal: after all
items are closed, a fresh review should find nothing that can sink the project. Update the status, evidence and date
whenever something changes. Statuses: **open**, **running**, **done** (analysis or text change made and checked), **accepted**
(cannot be fixed with this data; stated as a limitation), **needs decision** (waiting on the author).

Last updated: 2026-10-19.

## Summary

| Status | Items |
|---|---|
| done | R1, R2, R3, R4, R5, R6, R7, R8, R9, R11, R12, R13, R16, R17, R18 |
| future work | R10 (Reactome) |
| accepted | R14, R15 |

## Major risks

| ID | Risk | Status | Evidence so far | Next step |
|---|---|---|---|---|
| R1 | The primary comparison (Hallmark vs random, fixed named latent) is close to built in: in Hallmark the latent holds the interferon genes, in the random mask it holds other genes, and any latent can separate stimulated from control cells | done (README rewritten 2026-10-14) | Best-latent comparison already existed (+0.008). New set-only random mask (`setshuffle`, free-latent pool identical to real) scores 0.815 on the named latent; real − setshuffle +0.182 (95% CI −0.020 to +0.383, uncorrected Wilcoxon p = 0.002) | Report setshuffle next to the old random mask; compute best-latent comparison for setshuffle; reword takeaways 1 and 2 to what both nulls show |
| R2 | Ceiling: every method scores 0.97 to 0.998, so the test cannot discriminate | done (README rewritten 2026-10-14) | Best principal component 0.997; mean log expression of the 95 interferon-α genes 0.996 (no model, no labels) | Put the reference scores in the README and in the learning-curve figure; reword claims to "all methods separate the conditions"; add a complementary measure that does not saturate (R5, R8) |
| R3 | Free-latent pool confound in the random mask | done | Old random mask: 43 to 59 of 95 interferon-α genes in the free pool (0 in the real mask); best latent was a free latent in 3 of 10 seeds. `setshuffle` keeps the free pool identical (0 of 95 in all seeds) | Describe both nulls in the README (Methods and Limitations) |
| R4 | Stimulation is confounded with the 10x run (each condition is one library; no ambient-RNA correction) | done (limitation stays) | Best component after removing genes: 0.997 → 0.945 (216 annotated interferon-α/γ genes), 0.744 (top 500), 0.649 (top 1,000), 0.573 (top 2,000). The separation rests on a few hundred genes, not on a diffuse shift; run effects among those genes cannot be excluded | State the confound and this result in the README | Score the best principal component after removing (a) the annotated interferon-α/γ genes and (b) the most stimulation-responsive training genes; state the confound as a limitation (the data cannot fully exclude it) and add an independent dataset to Future work |
| R5 | The score measures concentration of the response in one latent, not capture; the no-mask VAE can spread it over many latents | done | Probe on all 58 latents: 0.998 to 0.999 for every variant (`diagnostics_probe.csv`), so every variant holds the response; the masks differ in whether one latent carries it | Reword "captures" to "in one latent" in the README |
| R6 | No simple baseline | done | `reports/tables/diagnostics_reference.csv`: best PC 0.997, mean interferon-α gene expression 0.996 | Cite in README (see R2) |

## Biology and wording

| ID | Risk | Status | Next step |
|---|---|---|---|
| R7 | The data-efficiency claim cannot be tested at a ceiling (no-mask VAE already 0.992 at 500 cells; random mean 0.48 at 500 cells); expiMap's sample-efficiency result was on a different task, so "tests this directly" overreaches | done (README rewritten 2026-10-14) | Add reference curves (PC, gene-set mean); reword to "cannot be detected on this task"; optionally add 100 and 250-cell sizes |
| R8 | Only 1 of 50 named latents is validated; "captures the biology" overgeneralises | done (README rewritten 2026-10-14) | Reword takeaway 2; optionally add a descriptive gene-loading enrichment for the best latents |
| R9 | Hallmark is not "pure knowledge": its sets are refined for coherent expression and were built for enrichment analysis; co-expression is therefore an alternative, not a null | done (README rewritten 2026-10-14) | Reframe as pre-specified named sets vs sets found in this data; drop "null" for co-expression and "expert-curated"; cite Liberzon 2015 after checking the abstract wording |
| R10 | A knowledge source independent of expression coherence (Reactome) | future work (decided 2026-10-19; stays in README Future work) | `masks_reactome.npz` exists. Needs a rule to pick about 50 non-redundant pathways of at least 12 genes, then 10 seeds per variant |
| R16 | The 'interferon-α latent' for an IFN-β stimulus; the interferon-α and -γ sets share 71 of 95 genes; lupus donors have raised baseline interferon activity | done (README rewritten 2026-10-14) | State specificity limits next to the claim, not only in Limitations |

## Statistics

| ID | Risk | Status | Next step |
|---|---|---|---|
| R11 | Equivalence of Hallmark and co-expression was removed from the script and its table deleted; "very similar" has no margin | done | Reinstated in `src/11_tests.py` with the margin ±0.005 AUROC, committed in 8599a48 (2026-10-06) before the review, so it was pre-specified. Equivalent (Holm-adjusted p = 0.0098; Hodges–Lehmann −0.0006, interval −0.0028 to +0.0002). Justification and the ceiling caveat are in the script header | Quote the margin, its justification and the ceiling caveat in the README Statistics section |
| R12 | Multiple testing: every p-value is uncorrected | done (README rewritten 2026-10-14) | Implemented in `src/11_tests.py` → `reports/tables/tests.csv` (10 tests in family Q1 with Holm; 18 in Q1b with Benjamini–Hochberg; distribution-free intervals at the matching level). Uncorrected p-values and intervals were removed from `08_evaluate.py`, `10_robustness.py` and their tables. Still to do: make the README quote `tests.csv` only. Rationale: Holm for the Q1 family (small confirmatory family, family-wise error); Benjamini–Hochberg for the learning curves (about 21 tests, same seeds across sizes so positively dependent). With 10 seeds the smallest Wilcoxon p is 0.00195, so Holm over more than 25 tests can never reject. Descriptive results (per cell type, per donor) get no p-values |
| R13 | The free-latent restriction was adopted after seeing results; the README says the design was fixed before the 10-seed runs | done (README rewritten 2026-10-14) | Say plainly in the README which choices were made after pilot runs |

| R17 | Figures (`overview.png`, `scores.png`, `learning_curves.png`) show four variants and no reference scores | done (2026-10-19) | All figures regenerated with five variants. Reference scores (best principal component 0.997, interferon-α gene mean 0.996) are quoted in the text and `diagnostics_reference.csv`, not drawn: as lines they sit under the Hallmark and co-expression dots and add nothing. The within-set random mask is full size only, so it is absent from the learning curves by design |

## Protocol: weaker-signal benchmark (fixed 2026-10-14, before any run)

Why: every method, with or without a VAE, scores 0.97 to 0.998 on the full data, so the main benchmark cannot separate
variants. Baselines as the strongest stimulation genes are removed (test cells): best principal component 0.863 and mean of
the remaining interferon-α genes 0.689 (0.724 when the library size is recomputed without the removed genes, as the models do) with the top 250 genes removed; 0.744 and 0.536 with the top 500 removed.
The question this benchmark asks: does a mask built on annotated sets help or hurt when the annotated set is incomplete and
the signal is weak?

1. **Genes removed:** the 250 genes with the largest |t| of stimulated vs control log expression in training cells (megakaryocytes excluded), chosen from data, never from the gene sets. List saved in `reports/tables/removed_genes_drop250.csv`.
2. **Masks rebuilt on the remaining genes:** the same 50 Hallmark sets (kept even if they shrink below 12 genes, to keep 58 latents), both shuffles, co-expression modules.
3. **Unchanged:** architecture, training settings, scoring rules (named latent for Hallmark and the two random masks, best latent chosen on validation labels for no mask and co-expression), donor covariate.
4. **Baselines at the same difficulty:** best principal component (chosen on validation labels) and mean of the remaining interferon-α genes.
5. **Pilot:** 5 variants × seeds 0 and 1 = 10 runs, drop250 only. Run with `SC_TAG=drop250`; outputs go to `results_drop250/` and `reports/tables_drop250/`, so the main experiment is untouched.
6. **Go/no-go for 10 seeds:** go if the spread among VAE variants exceeds their seed-to-seed range, or one variant is clearly apart from both baselines. If all variants sit at a baseline, report that and stop.
7. **If go:** 10 seeds at drop250 (and drop500 as a second level only if the pilot shows room), tests added to `src/11_tests.py` as their own family with the same correction rules. Reactome and held-out donors come after, on whichever benchmark discriminates.

| R18 | Both headline scores depend on choosing a latent, and the score saturates; the project needs measures that choose no latent and that show whether knowledge makes the network learn something | done (panel run, README rewritten around it, 2026-10-14) | Latent-free panel in `src/10_checks.py`: held-out loss (H), decoded stimulation shift correlated across genes (I), cell-type probe on the latents, full / 10 / 50 / 200 per type and two hard pairs (J), stimulation probe (G). Confirmatory family Q1c in `src/11_tests.py`: 4 measures × 3 contrasts (Hallmark minus each random mask and minus co-expression) = 12 tests, Holm. No-mask VAE, few-shot and hard-pair scores are descriptive. Fixed before any of it was run | Reactome (R10) and held-out donors remain open decisions |

## Accepted limitations

| ID | Limitation | Where stated |
|---|---|---|
| R14 | Test cells come from the same 8 donors as training cells; seeds cover one split only | README Limitations |
| R15 | One dataset, one stimulus, one target | README Limitations and Summary |

## Open decisions for the author

None. Reactome (R10) is future work.

## Change log

- 2026-10-14: review written; tracker created. Added `src/04b_build_set_shuffle.py`, `src/07b_setshuffle.sh` and `src/11_diagnostics.py`. Trained 10 seeds of `setshuffle`; scored with `src/08_evaluate.py` (variant added to 08 and 10). Committed the pre-review work as `e8a78ce`.
- 2026-10-14 (later): added `src/11_tests.py` and `reports/tables/tests.csv`; removed uncorrected p-values from 08 and 10; set-only random mask results added to the tables.
- 2026-10-14 (later): decided against 20 seeds; two-family correction stands. Added R4 (gene removal) and R5 (probe) results to `src/11_diagnostics.py` and `reports/tables/diagnostics_*.csv`. README rewrite pending; figures still show four variants (new item R17).
- 2026-10-14 (later): README rewritten for five variants, corrected statistics (quotes `tests.csv` only), reference scores, gene-removal and probe results, new Limitations; `docs/methods.md` updated. Equivalence margin closed at ±0.005 (pre-specified in 8599a48). Remaining: regenerate figures (R17), Reactome decision (R10), commit.
- 2026-10-14 (later): protocol for the weaker-signal benchmark fixed (above). Baseline ladder computed: see protocol.
- 2026-10-14 (pilot, drop250, seeds 0 and 1, 10 runs; tables in `reports/tables_drop250/`, runs in `results_drop250/`): baselines on this data are best PC 0.864 and mean of the 32 remaining interferon-α genes 0.724. In all 6 Hallmark, all-gene-shuffle and sets-only-shuffle runs the named interferon-α latent collapsed (posterior-mean variance ≤ 0.0001, so the named-latent score is 0.5) and the response sat in a free latent (best latent 0.855 to 0.943). Co-expression: best latent a free latent in both seeds (0.850, 0.924). No mask: best latent 0.820 and 0.898. Probe on all latents: 0.967 to 0.977 for the masked variants, 0.988 for no mask. Go/no-go rule is ambiguous here: go by the named-latent score (0.5 vs 0.82 to 0.92), no-go by the best-latent score (all overlap with each other and with the PC baseline). Awaiting the author's decision on the full 10 seeds.
- 2026-10-14 (consolidation): `src/` reduced to `scoring.py` (shared setup), `04_build_masks.py` (now also builds the sets-only shuffle), `07_experiment.sh` (all five variants; seeds and `SC_TAG` as options), `08_evaluate.py`, `10_checks.py` (merged robustness, diagnostics and baselines), `11_tests.py` (renamed from 12). Removed `04b`, `07b`, `07c`, `10_robustness`, `11_diagnostics`, `13_baselines`. Numbers must match the tables from before the merge; see the verification commands given to the author.
- 2026-10-14 (panel): wrote sections H, I, J in `src/10_checks.py` and family Q1c in `src/11_tests.py`; nothing run. Reactome rule and weaker-signal seeds wait for the panel results.
- 2026-10-14 (panel run and README): panel results in (Q1c tests: loss, decoded response, cell-type probe, stimulation probe). README rewritten around the panel; `decisions.md` D19 to D23 added; `09_figures.py` updated for five variants plus a panel figure. Remaining: author reruns `11_tests.py` (family label Q1c) and `09_figures.py`; Reactome rule; optional weaker-signal seeds and held-out donors.
- 2026-10-14 (weaker-signal run): purpose revised to test whether the Hallmark advantage in held-out loss and decoded response over the random masks persists without the 250 strongest stimulation genes (criterion in decisions.md D20). `11_tests.py` accepts `SC_TAG` and then tests only the latent-free family on that benchmark. Reactome and held-out donors parked.
- 2026-10-19: the two random masks renamed *fully random* (all genes shuffled, the control fixed before the runs) and *within-set random* (only set genes shuffled; the stricter control added after review); README states the reason for two masks. Figures regenerated; overview stage 3 now a 2 × 2 block of masks with no mask in a third column; stage 4 axis uses full names. `cell_types.png` and `training.png` removed (cited nowhere), `scores.png` embedded in supplement S1. Dashed reference lines tried and dropped.
