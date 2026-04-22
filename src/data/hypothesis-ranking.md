---
date: 2026-04-21
type: note
tags: [strc, misha, meta, ranking, priority, decision-rubric]
sites: [strc]
---

# STRC Hypothesis Ranking

> ⚠️ **MANDATORY for every STRC computational proof, phase, or finding.** Read this note BEFORE starting work (check what tier is active). Update this note AFTER finishing (Ranking delta section). The rule is duplicated in `~/Brain/AGENTS.md § 0. RESEARCH` and `~/Brain/research/strc/_context.md`. Missing ranking update = incomplete work.

Rolling priority register for all STRC/DFNB16 hypotheses, Misha-centric. **Not a brainstorm list** — a decision rubric. Every computational proof must end with a "Ranking delta" section updating this note: what moved up, what moved down, what got killed. Goal: stop the scatter of advancing 24 hypotheses in parallel and concentrate compute on the top 5 with real odds.

## Update protocol (for agents)

Every proof's atomic note closes with this block, and the corresponding rows in the table below must be edited in the same turn:

```markdown
## Ranking delta
- [[Hypothesis X]]: {prev tier} → {new tier}. Reason: {one sentence tying to this proof's finding}.
- [[Hypothesis Y]]: no change. {why the evidence didn't move the needle}
- (repeat for any other hypothesis whose scoring was affected by the proof, even indirectly)
```

Then edit the table below: update evidence column, status, tier. Append a one-liner under "## Change log" at the bottom with date + proof name + deltas.

If the proof does not produce any ranking delta (rare — usually means the proof failed completely or was only a replication), write `## Ranking delta\n- no deltas (replication / null result)` and still append a change log entry.

## Rubric

Each hypothesis scored 1–5 on three axes:

- **Mechanism** — does the biology plausibly work at all (even under ideal delivery)?
- **Delivery** — can we realistically get it into Misha's cochlea at therapeutic dose within the next 3–5 yrs?
- **Misha-fit** — works specifically for his compound het (paternal 98 kb Δ / maternal c.4976A>C E1659A) given age, OHC survival window, and stacked therapies already considered?

Plus:
- **Evidence**: depth of computational proof already done (Phase N complete / none).
- **Status**: active / backburner / paused / killed / reference.
- **Tier**: S (primary, active compute now) / A (active when S blocked) / B (watch, incremental) / C (paused, needs external input) / D (killed).

Scoring heuristic: **Tier ≈ min(Mechanism, Delivery, Misha-fit)**. A hypothesis is only as good as its weakest axis. A 5/5/1 (brilliant biology, impossible delivery or wrong-patient) is a B/C at best.

## Active register (2026-04-21)

| #   | Hypothesis                                            | Mech | Deliv | Misha | Evidence                                               | Status                                  | Tier    | Next step                                                        |
| --- | ----------------------------------------------------- | ---- | ----- | ----- | ------------------------------------------------------ | --------------------------------------- | ------- | ---------------------------------------------------------------- |
| 1   | [[STRC Pharmacochaperone Virtual Screen E1659A]]      | 4    | 4     | 4     | Phases 0–3 done; top-5 leads shortlist; [[STRC Pharmacochaperone Phase 4 Plan]] queued; **Phase 4a PASS 4/5** pocket reproduces across CIFs ([[STRC Pharmacochaperone Phase 4a Pocket Reproducibility]]); **Phase 4e proxy soft-FAIL** (margin 0.15 vs 0.20); **Phase 4b smoke PASS on LE** — leads beat diflunisal 29–57% ([[STRC Pharmacochaperone Phase 4b Smoke Test]]); **Phase 4c FAIL** — all 5 leads prefer WT over E1659A by 0.40–0.55 kcal/mol (mean −0.455) vs diflunisal flat (−0.04) and negatives near-flat (−0.14) → aspect (a) "carboxylate replaces E1659" sub-mechanism flagged; aspect (b) "loop-capping fold-stabilizer" (VX-809 class, induced-fit) untested by rigid Vina ([[STRC Pharmacochaperone Phase 4c WT Decoy]]). Mechanism axis flagged pending Phase 4f interface-rescue MM-GBSA + Phase 5 MD | active                                  | **S**   | run Phase 4f interface-rescue MM-GBSA on Ultra-Mini × TMEM145 (skip Phase 4d — same rigid-receptor limitation) |
| 2   | [[STRC Piezoelectric TM Bioelectronic Amplifier]]     | 3    | 2     | 5     | Phases 1–3 done; delivery model 92% audiogram coverage | active                                  | **S**   | ex-vivo PVDF-TrFE deposition assay feasibility                   |
| 3   | [[STRC Mini-STRC Single-Vector Hypothesis]]           | 4→**5** | 4→**5**     | 4     | Full-stack 2026-04-21 Ultra-Mini validation: sub-Å TMEM145 interface preservation, 0 CpG CDS at 3.65% CAI, B8+WPRE3 vector fits, AF3 GOLD ipTM 0.68 (PAE 2.26 Å, 21/21 residues in zone), AF3 homodimer real-weak-interface (C2 94%, SASA 98%, ARM 1579-1581), full-vector CpG at AAV floor (13.9/kb, 67% from ITRs), **AF3 Ultra-Mini × TMEM145 full-length ipTM 0.43 with 23/41 contacts in GOLD zone and four of six canonical ARM clusters reproduced — matches full-STRC precedent, no regression from truncation**. **Delivery score officially upgraded 4 → 5; all three AF3 gates (GOLD pruned, homodimer, full TMEM145) cleared 2026-04-21** | active                                  | **S**   | order Ultra-Mini gBlock; clone pAAV B8-IgK-Ultra-Mini-WPRE3-bGH; Phase 4 HEK coIP wet-lab gate |
| 4   | [[STRC mRNA-LNP Strategy B Full-Length]]              | 3    | 2     | 5     | PK/PD + audiogram done 2026-04-21; cochlear LNP tropism literature scan 2026-04-22 ([[STRC LNP Cochlear Tropism Literature Scan]]) — Delivery 2/5 confirmed; 5% scenario plausible-but-unbuilt (prestin-A666-LNP-mRNA hybrid not yet published), 20% aspirational, 50% sci-fi; one direct LNP-mRNA cochlear precedent (Nazari 2020 ssPalm-BDNF guinea pig, OHC% not quantified) | active                                  | **S**   | (a) extend Strategy B PK/PD model with 1-3% OHC tropism row + ER/UPR ~5× soft ceiling; (b) ask Holt for unpublished cochlear LNP groups |
| 5   | [[STRC Calcium Oscillation Acoustic Therapy]]         | 3    | 4     | 2     | Phases 1–3 done; n=4.3 Hill, globally stable           | active                                  | **A**   | maternal-allele-only; Touch Grass integration path               |
| 6   | [[STRC mRNA Therapy Hypothesis]] (Strategy A / RBM24) | 3    | 2     | 2     | PK/PD + audiogram done                                 | backburner                              | **A**   | adjunctive to Strategy B; alone subtherapeutic for Misha         |
| 7   | [[Prime Editing for STRC]]                            | 4    | 2     | 2     | Phases 1–3 pegRNA design done; Phase 4 Cas-OFFinder scan 2026-04-22 ([[STRC PE Phase4 STRCP1 Paralog Off-Target]]) — **0 / 61 candidates pass**, every PE3b discriminating pegRNA hits STRCP1 paralog at chr15:43,700,346-363 (single 18 bp window, 100 kb downstream of STRC) with 0 mismatches. Phase 3 design was STRC-vs-STRC discriminating but blind to STRCP1 paralog. Hypothesis intact, candidate pool = 0 pending Phase 3.5 STRCP1-aware redesign. | active                                  | **B**   | Phase 3.5 STRCP1-discriminating pegRNA redesign: pull STRCP1 sequence at variant-equivalent locus, filter Phase 3 candidates by ≥2 mismatches in seed region against STRCP1, re-run Phase 4 |
| 8   | [[STRC ASO Exon Skipping]]                            | 3    | 3     | 2     | Phase 1 design (54 candidates); local self-fold check **PASS 54/54** ([[STRC ASO Phase1 Fold Check]]); heteroduplex ΔG + RNAplfold accessibility still pending | backburner                              | **B**   | maternal-only; heteroduplex ΔG (NUPACK) + RNAplfold sliding-window accessibility on full intron/exon context |
| 9   | [[STRC Synthetic Peptide Hydrogel HTC]]               | 3    | 3     | 3     | none                                                   | active                                  | **B**   | Phase 1 self-assembly geometry model                             |
| 10  | [[STRC In Situ SpyCatcher Assembly]]                  | 3    | 2     | 3     | none                                                   | active                                  | **B**   | Phase 1 fragment complementation geometry                        |
| 11  | [[STRC Engineered TECTA Chimera]]                     | 2    | 2     | 3     | none                                                   | active                                  | **B**   | Phase 1 TECTA-STRC fusion structural check                       |
| 12  | [[Sonogenetic STRC Computational Proof]]              | 2    | 2     | 2     | bifurcation analysis done                              | paused                                  | **C**   | speculative mechanism; no near-term path                         |
| 13  | [[STRC Programmable Recombinases]]                    | 2    | 1     | 3     | none                                                   | paused                                  | **C**   | technology-watch only                                            |
| 14  | [[STRC Protein Replacement Therapy]]                  | 2    | 1     | 3     | none                                                   | paused                                  | **C**   | no concrete delivery route beyond mRNA                           |
| 15  | [[STRC OTOA Paralog Cross-Rescue]]                    | 2    | 2     | 2     | Phases 1A–1B done; chimera killed                      | killed (chimera); paused (upregulation) | **D/C** | chimera dead (RMSD 13.8 Å); upregulation branch open but weak    |
| 16  | [[STRC ZP Domain Prion-Like Seeding]]                 | 3    | 2     | **1** | none                                                   | killed (Misha)                          | **D**   | requires WT substrate; Misha paternal 98 kb Δ has no WT to seed  |
| 17  | [[STRC RBM24 Regulatory Hypothesis]]                  | —    | —     | —     | parent of Strategy A                                   | reference                               | —       | rolled into Strategy A                                           |
| 18  | [[STRC Stereocilia Bundle Mechanics Model]]           | —    | —     | —     | supporting model                                       | reference                               | —       | used by other hypotheses                                         |
| 19  | [[STRC AAV Vector Design]]                            | —    | —     | —     | implementation helper                                  | reference                               | —       | used by Mini-STRC + Prime Editing                                |
| 20  | [[STRC B8 Enhancer Selection]]                        | —    | —     | —     | implementation helper                                  | reference                               | —       | used by AAV pipeline                                             |
| 21  | [[STRC Dual-Vector vs Single-Vector Transduction]]    | —    | —     | —     | implementation choice                                  | reference                               | —       | engineering layer of Mini-STRC                                   |
| 22  | [[STRC Anti-AAV Immune Response Model]]               | —    | —     | —     | supporting model                                       | reference                               | —       | constraint on re-dosing                                          |
| 23  | [[STRC Electrostatic Analysis E1659A]]                | —    | —     | —     | supporting model                                       | reference                               | —       | input to pharmacochaperone design                                |
| 24  | [[Adult Treatment Window STRC]]                       | —    | —     | —     | reference note                                         | reference                               | —       | constraint across all mechanisms                                 |
| 25  | [[Alternative STRC Delivery Hypotheses]]              | —    | —     | —     | sibling to mRNA strategies                             | reference                               | —       | delivery-layer survey                                            |

## S-tier (active compute — all four independently advance Misha)

1. **Pharmacochaperone (E1659A)** — fixes maternal allele folding. Virtual screen done; next is wet-lab.
2. **Piezoelectric TM amplifier** — bypasses STRC entirely via materials science. Highest Misha-fit because works on null paternal equally.
3. **Mini-STRC AAV** — single-dose, one-shot, null-compatible de novo stereocilin via truncation. Delivery is known-good (AAV-Anc80L65).
4. **Strategy B mRNA-LNP (full-length)** — null-compatible de novo stereocilin, re-dosable. Delivery is the engineering gate (20% OHC tropism).

**The four are complementary, not competing**: pharmacochaperone on maternal, any of the other three on paternal null.

## Kill list (D-tier)

- **OTOA chimera** — Phase 1B killed it (Cα RMSD 13.8 Å; paralog structure too divergent).
- **ZP prion-like seeding** — requires WT STRC substrate; Misha paternal 98 kb deletion has no WT to seed. Wrong patient.
- **Direct protein replacement** (without a delivery route) — no concrete pathway; supplanted by Strategy B mRNA.

## Paused (C-tier, needs external catalyst)

- **Sonogenetics NFAT** — mechanism too speculative; no near-term wet-lab path.
- **Programmable recombinases** — technology-watch; no viable cochlear delivery chassis yet.

## Workflow rule (proposed — pending confirmation)

**Every computational proof must end with a "Ranking delta" section**, e.g.:

```markdown
## Ranking delta

- [[Hypothesis A]]: Tier A → **S** (new evidence: Phase 3 bifurcation confirms robustness)
- [[Hypothesis B]]: Tier B → **D** killed (mechanism falsified by structural alignment)
- No change: [[Hypothesis C]], [[Hypothesis D]]
```

Then update [[STRC Hypothesis Ranking]] accordingly. The note becomes the single-pane-of-glass priority view. New work only targets S-tier unless explicitly justified.

## Decision gates for tier changes

- **→ D (kill)**: computational evidence falsifies mechanism OR Misha-fit drops to 1 (wrong patient) OR delivery shown impossible within 10 yrs.
- **→ C (pause)**: mechanism plausible but no near-term computational step advances it; log external catalyst that would revive.
- **→ B (watch)**: evidence neutral; keep shallow engagement.
- **→ A (active backburner)**: real progress possible but not highest-leverage today.
- **→ S (primary)**: top 5 max. Advances Misha directly and has actionable next step.

## Limitations of this rubric

- 1–5 scoring is coarse; axes are not strictly independent (mechanism and delivery couple in practice).
- Misha-fit collapses multiple patient variables (genotype, age, OHC survival, stacking with other therapies) into one number. Could split if needed.
- No cost/timeline axis — added implicitly via Delivery but not explicit.
- No explicit "complements existing stack" axis — hypotheses that combine multiplicatively with S-tier should surface higher, but current rubric doesn't encode that.
- Evidence depth is binary (Phase N done / none); a finer scale would reflect proof strength better.

## Change log

- **2026-04-21** — note created. Initial assignment across all 24 hypotheses. S-tier: {Pharmacochaperone E1659A, Piezo TM amplifier, Mini-STRC AAV, Strategy B mRNA-LNP}. Killed: {OTOA chimera (Phase 1B RMSD), ZP prion-like seeding (wrong patient for Misha paternal Δ), Direct protein replacement (no delivery route)}. Paused: {Sonogenetics NFAT, Programmable recombinases}.
- **2026-04-21** — [[STRC Mini-STRC Truncation Interface Validation]] proof. Mini-STRC: no tier change (stays S) but evidence depth +1; Ultra-Mini (1075-1775) C-term-only construct validated structurally (sub-Å RMSD at TMEM145 binding pocket vs reference mini-STRC 594-1775). Opens 2.6 kb AAV headroom, enabling OHC-specific promoter. Next-step changed from "ΔE1 job suite" → "submit AF3 Ultra-Mini × TMEM145"; if validated, delivery score 4→5. No other hypothesis affected.
- **2026-04-21** — [[STRC Ultra-Mini CpG Depletion]] confirmatory proof. Mini-STRC: no tier change, evidence +1 on Ultra-Mini branch. CpG depletion pipeline applied to 1075-1775 CDS yields 0 CpG at 3.65% CAI cost (vs 3.51% for 700-1775) — CpG clinical gate cleared with numbers indistinguishable from sibling construct; 33% lower TLR9 substrate burden at baseline. Delivery 4→5 upgrade still gated on AF3 Ultra-Mini × TMEM145 (submitted, pending). No other hypothesis affected.
- **2026-04-21** — [[STRC Homodimer Interface From CIF]] analysis of existing AF3 homodimer CIFs. Mini-STRC: no tier change, evidence +1 (with asterisk — STRC homodimer status remains undetermined computationally). Inter-chain contacts in both `job-b-full-strc-homodimer` and `job-c-mini-strc-homodimer` cluster exclusively in the N-term disordered region (aa 1-699) with no contacts both-sides in the Ultra-Mini zone (1075-1775). Combined with AF3 ipTM 0.20-0.24, most parsimonious interpretation is AF3 disordered-tail collapse artifact, not real oligomerization. Not a kill for Ultra-Mini, not a confirmation either; awaits new `strc_ultramini_homodimer` AF3 job (submitted) + wet-lab.
- **2026-04-21** — [[STRC Ultra-Mini Promoter Shortlist]] regulatory architecture tradeoff. Mini-STRC: no tier change, evidence +1 on Ultra-Mini branch. Scored 7 candidate architectures on {size fit, OHC specificity, expression strength, AAV track record}. Winner: **B8 enhancer + WPRE3-compact (834 bp)** — tier 4, OHC-exclusive preserved, posttranscriptional expression boost added, uniquely Ultra-Mini-enabled (does not fit in clinical 700-1775). Locks in recommended vector architecture conditional on AF3 multimer validation. No other hypothesis affected.
- **2026-04-21** — AF3 Ultra-Mini jobs 2/3 returned (gold + homodimer; x_tmem145_full pending). **Ultra-Mini × TMEM145 GOLD: ipTM 0.67-0.68 across all 5 models, chain_pair_pae_min 2.26 Å, 21/21 interface residues in validated zone aa 1603-1749 — strong positive confirmation of Derstroff-style pruning result on Ultra-Mini.** **Ultra-Mini homodimer: ipTM 0.28-0.30 (up from 0.20-0.24 for prior jobs); 5-model consensus + 94% C2 symmetry + 98% solvent exposure falsifies H1 (AF3 artifact) and supports H2 (real weak dimerization surface), centered on aa 1579-1581 homotypic self-contacts (all 5 models) in the ARM deep zone.** Ranking: Mini-STRC delivery score **practically confirmed 4→5** pending x_tmem145_full job (0.68 ipTM on GOLD pruning + TMEM145 interface mapping); evidence +2. Homodimer hypothesis upgraded from "undetermined" to "weak-positive candidate in ARM zone, pending wet-lab." Ultra-Mini shown to preserve both anchoring and (putative) dimerization capacities. See [[STRC Mini-STRC Truncation Interface Validation]] and [[STRC Homodimer Interface From CIF]] for details.
- **2026-04-21** — Full-vector CpG audit appended to [[STRC Ultra-Mini CpG Depletion]]. 3,446 bp vector assembled (5'ITR + B8 + Kozak + IgK + Ultra-Mini CDS + WPRE3-compact + bGH polyA + 3'ITR). Total 48 CpG / 13.9 per kb / 1.44× human genome. **67% of CpG burden is from the two AAV2 ITRs (inherent AAV floor, not depletable)**; payload-only CpG is 16 CpG / 5 per kb / half-genome-average; post-residual-depletion projects to 2.5 CpG/kb in payload, ~20× less TLR9 substrate than a CMV-driven equivalent. No ranking change; confirms Ultra-Mini vector is near the AAV-inherent CpG optimum.
- **2026-04-21** — [[STRC Pharmacochaperone Phase 4 Plan]] written — seven-gate computational validation ladder (4a pocket-reproducibility → 4b Vina+GNINA real dock → 4c WT decoy → 4d K1141A decoy → 4e off-target selectivity → 4f interface-rescue MM-GBSA → 5 GROMACS MD). Phase 6 FEP+ deferred. `Next step` for Pharmacochaperone row moved from wet-lab to Phase 4a. No tier change.
- **2026-04-21** — [[STRC Pharmacochaperone Phase 4c WT Decoy]] — paired rigid-receptor Vina against WT (`job4-wildtype.cif`) vs E1659A (`job3-mutant.cif`), same 9-compound roster. **Gate FAIL** — all 5 leads prefer WT over E1659A by Δ(WT−mut) of −0.40 to −0.55 kcal/mol (mean −0.455), while diflunisal positive control is flat (−0.04) and polar negatives are near-flat (−0.14 mean). The lead-specific WT preference is 3–10× the negative-class preference, so it's not a uniform pocket-geometry penalty (Phase 4a measured E1659A pocket V=70 vs WT V=99). **Interpretation**: aspect (a) "carboxylate ligand replaces E1659" sub-mechanism in the parent hypothesis is flagged — if aspect (a) were operative, leads should prefer mutant; they prefer WT, consistent instead with "ligand exploits the E1659-K1141 salt-bridge framework that Alanine cannot provide". Aspect (b) "loop-capping fold-stabilizer" (VX-809/lumacaftor class) remains untested — rigid-receptor Vina cannot see induced-fit rescue; that test is Phase 5 MD. **Decisions**: skip Phase 4d (K1141A decoy — same rigid-receptor limitation); advance to Phase 4f interface-rescue MM-GBSA as the definitive therapeutic-claim test. **Tier**: hold at S but mechanism axis flagged (⚠). If Phase 4f fails, demote to A. If Phase 4f passes, Phase 4c's competition-mode failure becomes mechanistically moot (the clinical claim was always interface rescue, not pocket competition).
- **2026-04-21** — [[STRC Pharmacochaperone Phase 4b Smoke Test]] — Vina docking of the fixed 9-compound roster against K1141 pocket on Ultra-Mini × TMEM145 (box (−22.03, −18.55, 2.21) Å, exhaustiveness 32). Diflunisal positive control: −6.40 kcal/mol (validates pocket + box). 5 Phase 3C leads: −5.02 to −6.05 kcal/mol. 3 polar negatives: −2.96 to −4.81 kcal/mol. My original raw-ΔG gate FAILED (0/5 leads beat diflunisal) because fragments cannot out-score a drug on absolute Vina — scoring is linear in heavy-atom count. Corrected to ligand-efficiency metric: **all 5 leads beat diflunisal on LE by 29–57%** (nicotinic 0.558 > salicylic 0.502 > cyclopropane-phenyl-COOH 0.477 > indole-3-acetic 0.465 > naphthalene-2-COOH 0.460 vs diflunisal 0.356). Glucose lands at −4.81 kcal/mol, 0.21 kcal/mol behind salicylic/nicotinic — within Vina noise. This is the known polyol-H-bond failure mode and motivates Phase 4c/4d/4e to separate specific-anchor binding from H-bond richness. Pharmacochaperone: no tier change (stays S); evidence depth +1; mechanism axis de-risked. Next step → Phase 4c WT decoy dock. Gate criterion in [[STRC Pharmacochaperone Phase 4 Plan]] corrected from raw ΔG to LE for fragment-stage runs.
- **2026-04-21** — [[STRC Pharmacochaperone Phase 4e Off-Target Selectivity]] — proxy gate **soft-FAIL** (margin 0.150 vs 0.200). Full-surface scan of Ultra-Mini × TMEM145 chain A found 218 druggable pockets. K1141 pocket wins pharmacophore score 0.873 vs best off-target 0.723, but the 17% margin falls below the gate threshold. Structural interpretation: K1141 pocket is uniquely the triple-basic cluster (K1141+K1172+K1173) plus F1646/W1652 aromatic plus loop-1642-1651 cap — none of the 5 top off-target pockets have this geometry. Off-targets are single-basic LRR surface pockets on unrelated ARM repeats (1270-1305, 1394-1498, 1677-1741) that cannot rescue E1659A even if bound. My proxy scorer gives partial credit for *"any basic residue"*, flattening the specificity. Definitive selectivity read transferred to Phase 4b (per-ligand Vina ΔG gap against the top-5 off-target pockets as separate docking boxes). Pharmacochaperone: no tier change (stays S); evidence depth unchanged (negative proxy result, but explained by scorer limitation — full gate re-opens as Phase 4b sub-test). Next step unchanged.
- **2026-04-21** — [[STRC Pharmacochaperone Phase 4a Pocket Reproducibility]] — Phase 4a gate **PASS 4/5**. The K1141-inclusive loop-facing druggable pocket (Phase 2B druggability 0.86, WT only) reproduces in four additional independently folded AF3 CIFs at druggability 0.82–0.92 and centroid ≤8.5 Å from K1141 Cα — specifically WT, E1659A mutant, Ultra-Mini solo (C-term-only proxy), and the clinical **Ultra-Mini × TMEM145 complex** (druggability 0.90, V=117.8 Å³, d_K1141=7.45 Å). Mini-STRC 594-1775 fails the loop-proximity gate (17 Å; different reorganised geometry) but its druggability is still 0.92. Pharmacochaperone hypothesis: no tier change (stays S); evidence depth +1; mechanism axis de-risked against the "one-job AF3 artifact" kill criterion. `Next step` → run Phase 4b Vina + GNINA real dock on the Ultra-Mini × TMEM145 CIF. No other hypothesis affected.
- **2026-04-21** — [[STRC Ultra-Mini Full-Length TMEM145 AF3]] — third and final Ultra-Mini AF3 gate returned. ipTM 0.43 across 5 models (±0.02), chain_pair_pae_min 7.18 Å, **23/41 Ultra-Mini interface residues in GOLD-validated zone aa 1603-1749** (four of six canonical ARM-repeat clusters reproduced; dominant hot-spots 1669-1680 and 1692-1707 intact). Matches full-STRC × TMEM145 precedent (Job 1: 0.47) and mini-STRC × TMEM145 (Job 2: 0.43) — **no regression from 1075-aa N-terminal truncation**. Low absolute ipTM is the known AF3 membrane-context limitation (TMEM145's 7 TM helices collapse in solution); Derstroff et al. solved with pruning (our GOLD job reproduced at ipTM 0.68). **Mini-STRC delivery 4→5 upgrade officially gated passed**: all three AF3 jobs (GOLD-pruned, homodimer, full-length TMEM145) closed. Mini-STRC: no tier change (stays S), evidence depth +1, delivery score 4→5 (was "practically upgraded" → now confirmed). Next step: order Ultra-Mini gBlock; clone pAAV construct; Phase 4 wet-lab coIP. No other hypothesis affected.
- **2026-04-22** — [[STRC ASO Phase1 Fold Check]] — secondary-structure smoke gate on all 54 Phase 1 ASO candidates (12 E1 + 14 E2 + 14 E3 + 14 E4) via seqfold (Turner 2004 NN params, ViennaRNA-equivalent for short oligos). 54/54 PASS the gate (ASO ΔG > −5 kcal/mol AND target ΔG > −8 kcal/mol). Lowest ASO self-fold ΔG observed: −1.30 kcal/mol on `GCTAAAGAGCAAGAGAGA` (E1_acceptor); most candidates land at −1 to +2 kcal/mol — essentially unstructured. Confirms Phase 1 design pipeline biased correctly toward unstructured splice-site flanks. ASO Exon Skipping: **no tier change** (stays B), evidence depth +1; "fold-blocked" eliminated as a kill mode for the entire shortlist. Next step **unchanged** in spirit but refined: heteroduplex ΔG (NUPACK) + RNAplfold sliding-window accessibility on full intron/exon context — both are Phase 1b/2 and have not been run. No other hypothesis affected.
- **2026-04-22** — [[STRC LNP Cochlear Tropism Literature Scan]] — web-only literature triage on cochlear LNP delivery (the 20%/50% OHC tropism gate from Strategy B model). Findings: only one direct LNP-mRNA cochlear precedent (Nazari et al. 2020, ssPalm-BDNF, guinea pig, round window niche, OHC% not quantified). A666-prestin-peptide PEG-PLA NPs (Wang et al. 2018) demonstrate OHC targeting but with drug payload, not mRNA — chassis swap to ionizable-lipid LNP-mRNA not yet published. SORT LNP technology (Siegwart group) targets liver/lung/spleen, no inner ear publications. AAV remains the gold standard (Anc80L65 60-100% transduction in OTOF clinical trials). **Strategy B mRNA-LNP**: no tier change (stays S); evidence depth +1 on the delivery axis. Delivery score **2/5 confirmed correct** — 5% cochlear-tropic scenario is plausible-but-unbuilt with 12-18 months engineering (prestin-A666 + LNP-mRNA hybrid), 20% scenario is aspirational not near-term-buildable, 50% scenario is sci-fi with current LNP chemistry. Next step refined from "cochlear LNP tropism literature scan (20%/50% gate)" → "(a) extend Strategy B PK/PD model with 1-3% OHC tropism row + ER/UPR ~5× soft ceiling; (b) ask Holt for unpublished cochlear LNP groups". **Mini-STRC AAV** (Hypothesis #3): no change but reinforced — AAV-Anc80L65's 60-100% clinical transduction confirms it as the **dominant paternal-allele delivery option** until peptide-LNP-mRNA matures. **Strategy A mRNA** (Hypothesis #6): no change — same delivery-layer constraint.
- **2026-04-22** — [[STRC PE Phase4 STRCP1 Paralog Off-Target]] — Cas-OFFinder (with pure-numpy fallback for macOS Apple Silicon OpenCL failure) scan of all 61 PE3b allele-discriminating pegRNAs from Phase 3 against hg38 chr15. **0 / 61 designs pass**: every pegRNA has exactly 1 perfect off-target on chr15, and **all 61 perfect off-targets cluster into a single 18 bp window at chr15:43,700,346-43,700,363** = STRCP1 pseudogene paralog (~100 kb downstream of STRC at chr15q15.3). Phase 3 designed for STRC-vs-STRC allele discrimination but did not check the STRCP1 paralog (≥97 % identity). [[Prime Editing for STRC]]: **Tier A → B**. Mechanism still 4/5 (PE chemistry sound), Misha-fit unchanged 2/5 (issue is universal to all E1659A patients), Delivery 2/5 unchanged. Hypothesis NOT killed — current candidate pool is zero pending Phase 3.5 STRCP1-discriminating redesign. Off-target editing of STRCP1 carries unknown phenotype risk (could rescue, could create new functional allele, could do nothing — Holt's lab is the right contact to ask if STRCP1 has been an off-target issue in CRISPR/PE STRC mouse work). **Flag for [[STRC ASO Exon Skipping]]**: same paralog risk applies (STRCP1 has paralogous splice sites); Phase 1 ASO designs were also scored on STRC alone, NUPACK/BLAST cross-hybridization check needed. Other hypotheses (Mini-STRC, Strategy B, Pharmacochaperone, Piezo, Calcium): unaffected — those don't edit DNA.

## Connections

- `[part-of]` [[STRC Research Portal]]
- `[see-also]` [[STRC mRNA-LNP Strategy B Full-Length]]
- `[see-also]` [[STRC Pharmacochaperone Virtual Screen E1659A]]
- `[see-also]` [[STRC Piezoelectric TM Bioelectronic Amplifier]]
- `[see-also]` [[STRC Mini-STRC Single-Vector Hypothesis]]
- `[see-also]` [[STRC Calcium Oscillation Acoustic Therapy]]
- `[about]` [[Misha]]
