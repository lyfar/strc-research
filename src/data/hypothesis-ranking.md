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
| 1   | [[STRC Pharmacochaperone Virtual Screen E1659A]]      | 4    | 4     | 4     | Phases 0–3 done; top-5 leads shortlist; [[STRC Pharmacochaperone Phase 4 Plan]] queued (seven computational gates before any wet spend) | active                                  | **S**   | run Phase 4a pocket reproducibility gate (5-CIF fpocket-style scan) |
| 2   | [[STRC Piezoelectric TM Bioelectronic Amplifier]]     | 3    | 2     | 5     | Phases 1–3 done; delivery model 92% audiogram coverage | active                                  | **S**   | ex-vivo PVDF-TrFE deposition assay feasibility                   |
| 3   | [[STRC Mini-STRC Single-Vector Hypothesis]]           | 4→**5** | 4→**5**     | 4     | Full-stack 2026-04-21 Ultra-Mini validation: sub-Å TMEM145 interface preservation, 0 CpG CDS at 3.65% CAI, B8+WPRE3 vector fits, AF3 GOLD ipTM 0.68 (PAE 2.26 Å, 21/21 residues in zone), AF3 homodimer real-weak-interface (C2 94%, SASA 98%, ARM 1579-1581), full-vector CpG at AAV floor (13.9/kb, 67% from ITRs), **AF3 Ultra-Mini × TMEM145 full-length ipTM 0.43 with 23/41 contacts in GOLD zone and four of six canonical ARM clusters reproduced — matches full-STRC precedent, no regression from truncation**. **Delivery score officially upgraded 4 → 5; all three AF3 gates (GOLD pruned, homodimer, full TMEM145) cleared 2026-04-21** | active                                  | **S**   | order Ultra-Mini gBlock; clone pAAV B8-IgK-Ultra-Mini-WPRE3-bGH; Phase 4 HEK coIP wet-lab gate |
| 4   | [[STRC mRNA-LNP Strategy B Full-Length]]              | 3    | 2     | 5     | PK/PD + audiogram done 2026-04-21                      | active                                  | **S**   | cochlear LNP tropism literature scan (20%/50% gate)              |
| 5   | [[STRC Calcium Oscillation Acoustic Therapy]]         | 3    | 4     | 2     | Phases 1–3 done; n=4.3 Hill, globally stable           | active                                  | **A**   | maternal-allele-only; Touch Grass integration path               |
| 6   | [[STRC mRNA Therapy Hypothesis]] (Strategy A / RBM24) | 3    | 2     | 2     | PK/PD + audiogram done                                 | backburner                              | **A**   | adjunctive to Strategy B; alone subtherapeutic for Misha         |
| 7   | [[Prime Editing for STRC]]                            | 4    | 2     | 2     | Phase 1–3 pegRNA design done                           | active                                  | **A**   | maternal-allele only; Cas-OFFinder off-target scan               |
| 8   | [[STRC ASO Exon Skipping]]                            | 3    | 3     | 2     | Phase 1 design; lead ASO                               | backburner                              | **B**   | maternal-only; ViennaRNA fold check                              |
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
- **2026-04-21** — [[STRC Ultra-Mini Full-Length TMEM145 AF3]] — third and final Ultra-Mini AF3 gate returned. ipTM 0.43 across 5 models (±0.02), chain_pair_pae_min 7.18 Å, **23/41 Ultra-Mini interface residues in GOLD-validated zone aa 1603-1749** (four of six canonical ARM-repeat clusters reproduced; dominant hot-spots 1669-1680 and 1692-1707 intact). Matches full-STRC × TMEM145 precedent (Job 1: 0.47) and mini-STRC × TMEM145 (Job 2: 0.43) — **no regression from 1075-aa N-terminal truncation**. Low absolute ipTM is the known AF3 membrane-context limitation (TMEM145's 7 TM helices collapse in solution); Derstroff et al. solved with pruning (our GOLD job reproduced at ipTM 0.68). **Mini-STRC delivery 4→5 upgrade officially gated passed**: all three AF3 jobs (GOLD-pruned, homodimer, full-length TMEM145) closed. Mini-STRC: no tier change (stays S), evidence depth +1, delivery score 4→5 (was "practically upgraded" → now confirmed). Next step: order Ultra-Mini gBlock; clone pAAV construct; Phase 4 wet-lab coIP. No other hypothesis affected.

## Connections

- `[part-of]` [[STRC Research Portal]]
- `[see-also]` [[STRC mRNA-LNP Strategy B Full-Length]]
- `[see-also]` [[STRC Pharmacochaperone Virtual Screen E1659A]]
- `[see-also]` [[STRC Piezoelectric TM Bioelectronic Amplifier]]
- `[see-also]` [[STRC Mini-STRC Single-Vector Hypothesis]]
- `[see-also]` [[STRC Calcium Oscillation Acoustic Therapy]]
- `[about]` [[Misha]]
