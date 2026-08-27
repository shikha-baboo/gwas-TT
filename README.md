---
title: GWAS Explorer
emoji: 🧬
colorFrom: green
colorTo: blue
sdk: gradio
sdk_version: 4.44.0
app_file: app.py
pinned: false
---

# GWAS Explorer — User Manual: Settings, What They Mean, What to Pick

This is a practical guide, not a code walkthrough. For every setting in the
app, it explains **what it controls**, **what happens if you set it too
high/too low**, and **what range to actually use** for a typical rice
breeding panel.

---

## 1. Before you touch any settings — data checks

- Genotype coding must be consistent across the whole file: either all
  `0/1/2` or all `-1/0/1`. Don't hand-mix them.
- Missing data: only `-9`, `-99`, `-999` are read as missing. If your file
  uses a different missing code, recode it before uploading, or it'll be
  silently treated as a real genotype call.
- MAF filtering happens *before* any model runs — a marker dropped here
  never appears in any plot or the Excel output.

---

## 2. Tab 1 — Core run settings

### MAF Threshold — **default 0.05, typical range 0.03–0.10**
Drops markers whose minor allele frequency is below this value.
- **Too low (e.g. 0.01)**: keeps very rare-allele markers, which are
  statistically unstable — a handful of individuals can swing the p-value
  a lot, inflating false positives.
- **Too high (e.g. 0.15+)**: throws away real, usable markers, especially
  in panels with several minor subpopulations.
- **Recommendation**: 0.05 is a safe default for most diversity/breeding
  panels (100+ lines). Drop to 0.03 if your panel is small (<80 lines) and
  you don't want to lose markers; raise to 0.08–0.10 for very large,
  diverse panels (300+ lines) where rare-allele noise is more of a
  problem.

### Fallback Top-N — **default 15, range 10–20**
Only used if the permutation threshold fails to compute. It's a safety
net, not a real significance rule — don't rely on it to decide anything.
Leave at default; if it's actually being used, something upstream (usually
too few markers or too few permutations) needs fixing instead.

### Permutations for LOD threshold — **default 100, range 100–200**
Number of phenotype shuffles used to build the empirical significance
cutoff (Churchill & Doerge method). This is the real basis for calling a
marker significant.
- **Too low (<50)**: threshold estimate is noisy — rerunning the pipeline
  can give you a visibly different threshold each time.
- **Too high (>200, capped internally anyway)**: diminishing returns, just
  slower.
- **Recommendation**: 100 is fine for exploratory work. Use 200 (the max)
  for a threshold you're going to report/publish — it's more stable and
  the app caps compute anyway, so there's no real reason to use fewer.

### Manual LOD threshold — **default 0, i.e. off (automatic)**
Overrides the permutation threshold with a fixed LOD value you set
yourself.
- **When to use**: only when you have a specific external reason to fix
  the threshold — e.g. matching a threshold used in a companion paper/QTL
  study, or comparing runs on unequal marker sets where the permutation
  threshold isn't comparable across them.
- **Typical values if you do use it**: LOD 3–4 is the conventional QTL
  mapping range; for GWAS with thousands of markers you'd generally want
  it higher (LOD ~4–6) to stay close to a Bonferroni-equivalent stringency.
- **Default advice**: leave at 0 and let permutation set it — that's the
  statistically defensible choice for a given dataset.

### Minimum PVE % — **default 0, i.e. off**
Extra filter: a marker must explain at least this % of phenotypic variance
*in addition to* passing the LOD threshold.
- **When to use**: when you specifically want to report only
  large-effect MTAs (e.g. for marker development / MAS candidates), not
  every marker that's merely statistically significant.
- **Typical range if enabling**: 5–10% for a moderately polygenic trait;
  10–15%+ if you specifically want major-effect QTNs only.
- **Caution**: setting this too high (>20%) on a polygenic trait can wipe
  out every hit — most real QTNs for complex traits explain a few percent
  each, not 20%.

---

## 3. Model selection — which of the 8 to run

You don't need all 8 every time.

- **Always include one baseline + one kinship-corrected model** so you can
  see the difference structure correction makes: **OLS + EMMAX** (or MLM)
  is the minimum sensible pair.
- **MLM vs EMMAX**: MLM adds fixed kinship-PC covariates on top of the
  kinship random effect; EMMAX uses kinship alone. If your panel has clear
  subpopulation structure (e.g. indica/japonica mix), MLM's extra Q
  covariates usually help. If it's a fairly uniform breeding population,
  EMMAX alone is often enough and slightly less conservative.
- **GEMMA**: the most rigorous single-locus model (per-marker REML), but
  the slowest. Use it as a confirmatory check on your top hits from
  MLM/EMMAX rather than as a first-pass scan on a huge marker set.
- **FarmCPU / BLINK**: multi-locus methods — better power when several
  QTNs are linked or when polygenic background masks individual effects
  under single-locus models. Slower, and iterative, so results can shift
  slightly run-to-run depending on cofactor settings (see below). Good
  second-pass once you know roughly what OLS/MLM/EMMAX are finding.
- **mrMLM / FASTmrMLM**: also multi-locus, but with an explicit two-stage
  screen → joint-fit → backward-elimination design. Useful if you want an
  independent multi-locus method to cross-check FarmCPU/BLINK hits against
  (different QTN-selection logic → different false-positive profile).
  FASTmrMLM trades some rigor (drops kinship in the joint-fit stage) for
  speed — use it for quick multi-locus scans, mrMLM for the more careful
  version.

**Practical recommendation for a standard run**: OLS + MLM + EMMAX +
FarmCPU + GEMMA (the app's default selection) gives you baseline,
single-locus-corrected, and multi-locus coverage without running
everything. Add BLINK/mrMLM/FASTmrMLM only when you want to
cross-validate a specific set of candidate MTAs.

---

## 4. Cofactors / Advanced tab

### FarmCPU & BLINK — pseudo-QTN iteration

**Max FEM/REM iterations — default 8, range 5–12**
How many rounds of "scan genome → pick pseudo-QTNs → re-check them" the
method runs before stopping (it also stops early if the QTN set stops
changing).
- **Too few (<4)**: may stop before the pseudo-QTN set has converged,
  under-correcting for background genetic effects.
- **Too many (>12)**: rarely helps — if it hasn't converged by ~10
  iterations it usually won't, and you're just adding runtime.
- **Recommendation**: 8 is a good default; drop to 5–6 for a quick look on
  a large marker set, raise to 10–12 only if you notice the results
  changing meaningfully when you bump it up from 8 (a sign it hadn't
  converged).

**FarmCPU bin size (bp) — default 1,000,000 (1 Mb), range 500 kb – 2 Mb**
Genome is chopped into bins this wide; only the single best candidate per
bin becomes a pseudo-QTN each round.
- **Too small (e.g. 100 kb)**: can select several near-duplicate
  pseudo-QTNs from the same LD block as if they were independent —
  inflates the covariate set with redundant markers.
- **Too large (e.g. 5 Mb)**: can merge two genuinely distinct, linked QTLs
  into one bin and only keep one of them.
- **Recommendation**: 1 Mb is a reasonable default for rice (LD generally
  decays within a few hundred kb in most cultivated panels). If your LD
  decay plot (Tab 4) shows LD extending past ~1 Mb, consider widening the
  bin to 1.5–2 Mb; if LD decays very fast (<200 kb), you can narrow it to
  500 kb.

**BLINK LD r² pruning cutoff — default 0.7, range 0.5–0.9**
Two candidates in the same local LD block are treated as redundant (and
the weaker dropped) once their r² exceeds this.
- **Too low (e.g. 0.3)**: over-prunes — treats markers that aren't really
  redundant as duplicates, potentially discarding a real second QTN that
  happens to sit near a stronger one.
- **Too high (e.g. 0.95)**: under-prunes — keeps near-identical markers as
  separate pseudo-QTNs, which can waste degrees of freedom and slow
  convergence.
- **Recommendation**: 0.7 is the standard LD-pruning cutoff used in most
  GWAS pipelines (similar to PLINK's common defaults). Only move it if you
  have a specific reason — e.g. very high average LD in your panel, where
  0.8–0.9 avoids over-pruning.

### mrMLM & FASTmrMLM — two-stage screening

**Stage-1 screening p-value threshold — default 0.01, range 0.001–0.05**
How lenient the initial scan is about calling something a "candidate"
worth testing jointly in Stage 2.
- **Too strict (0.001)**: may exclude real but modest-effect QTNs before
  they even get a chance in the joint-fit stage.
- **Too lenient (0.05)**: floods Stage 2 with weak candidates, most of
  which get backward-eliminated anyway — just adds runtime and noise.
- **Recommendation**: 0.01 is a sensible middle ground; loosen to 0.02–0.05
  if you suspect real signals are being screened out (e.g. a highly
  polygenic trait with many small effects), tighten to 0.001–0.005 for a
  cleaner, more conservative candidate list.

**Max candidate QTNs carried into Stage 2 — default 20, range 10–30**
Caps how many Stage-1 candidates get fit jointly.
- **Too low (<10)**: may cut off real candidates before backward
  elimination gets a chance to sort them out, especially on a trait with
  many markers passing Stage 1.
- **Too high (>30)**: joint model becomes large and can get numerically
  unstable, especially with a small sample size relative to candidate
  count.
- **Recommendation**: 20 works for most panels. Scale down toward 10–15 if
  your sample size is small (<100 lines) — you don't have the degrees of
  freedom to support a 20-marker joint model reliably.

**Stage-2 backward-elimination drop threshold — default 0.05, range
0.01–0.1**
In the joint fit, the weakest surviving candidate is dropped (and the
model refit) if its p-value exceeds this, repeated until everyone left
clears the bar.
- **Too strict (0.01)**: can eliminate real co-segregating QTNs that only
  look weak because they're sharing variance with a correlated candidate.
- **Too lenient (0.1)**: leaves weaker/noisier candidates in the final
  set.
- **Recommendation**: 0.05 (standard significance level) is the sensible
  default; there's rarely a good reason to move it unless you're
  deliberately being more/less conservative about the final candidate
  list.

### Kinship-PC covariates (per model) — default 3 for MLM/mrMLM/FASTmrMLM, 0 for GEMMA; range 0–6

Number of top kinship eigenvectors used as fixed "Q" structure covariates.
- **0**: disables it — the model relies on kinship (K) alone, or on the
  method's own structure handling (e.g. GEMMA's per-marker REML).
- **Too many PCs (8–10)**: starts absorbing real trait-associated genetic
  variance into the covariates, reducing power to detect true
  associations (over-correction).
- **Recommendation**: 3 PCs is the standard starting point and works for
  most panels with 2–4 broad subpopulations. If a PCA plot (Tab 4) shows
  more than ~3–4 distinct clusters, consider 4–5 PCs. For GEMMA, leaving
  it at 0 is usually right — GEMMA's per-marker REML already handles
  structure, and adding PCs on top is somewhat redundant. Never go above
  6 unless you have a very large, strongly structured panel and have
  checked that it doesn't kill all your signal.

---

## 5. Quick-reference table

| Setting | Default | Typical range | Raise it when… | Lower it when… |
|---|---|---|---|---|
| MAF threshold | 0.05 | 0.03–0.10 | large/diverse panel, noisy rare alleles | small panel, don't want to lose markers |
| Permutations | 100 | 100–200 | reporting/publishing a threshold | quick exploratory look |
| Manual LOD | 0 (off) | 3–6 if used | matching an external threshold | (leave off — use automatic) |
| Min PVE % | 0 (off) | 5–15 if used | want major-effect QTNs only | polygenic trait, don't want to lose real hits |
| FarmCPU/BLINK iterations | 8 | 5–12 | results still shifting at 8 | quick pass on large data |
| FarmCPU bin size | 1 Mb | 0.5–2 Mb | slow LD decay | fast LD decay |
| BLINK LD r² cutoff | 0.7 | 0.5–0.9 | very high panel-wide LD | rarely — 0.7 is standard |
| mrMLM/FASTmrMLM screen p | 0.01 | 0.001–0.05 | suspect missed weak QTNs | want a cleaner candidate list |
| Max Stage-2 candidates | 20 | 10–30 | large sample, many Stage-1 hits | small sample (<100 lines) |
| Stage-2 drop threshold | 0.05 | 0.01–0.1 | rarely — 0.05 is standard | want stricter final QTN set |
| Kinship PCs | 3 | 0–6 | >3–4 clusters in PCA plot | simple/uniform panel, or using GEMMA |

If in doubt: run once with all defaults, look at the QQ plots (λ close to
1 = well-controlled) and the PCA plot (how many clusters?), then adjust
kinship-PC count and MAF threshold based on what you see — those two have
the biggest practical effect on your results.
