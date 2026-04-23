---
date: 2026-04-21
type: note
tags: [strc, misha, meta, ranking, priority, decision-rubric]
sites: [strc]
---

# STRC Hypothesis Ranking

> MANDATORY before any STRC work. Read tier table. Update after every proof. Missing update = incomplete work.

**Update protocol:** see [[AGENTS]] §0. TL;DR — `## Ranking delta` in proof note → prepend to `h{N}/log.md` → edit `h{N}/hub.md` frontmatter if tier/score/next-step changed. Dataview re-renders this page automatically. No sync scripts. Architecture: [[Vault Stack and Schemas]].

## Rubric

Axes 1-5: **Mech** (biology plausible), **Deliv** (cochlear delivery 3-5 yr), **Misha-fit** (compound het, age, OHC window).

Tier heuristic: `min(Mech, Deliv, Misha-fit)`. S = top 5 active. A = backburner. B = watch. C = paused, needs external catalyst. D = killed.

## Active register (auto-rendered from `h{N}/hub.md` frontmatter)

Canonical view. To change a row — edit the hub's frontmatter (`tier`, `mech`, `deliv`, `misha_fit`, `next_step`), not this file. All 18 numbered hypotheses (h01–h16, h26, h27) are covered. Ref rows (supporting notes, clinical plans) are listed separately below.

```dataview
TABLE WITHOUT ID
  hypothesis_num AS "#",
  link(hypothesis_title) AS "Hypothesis",
  mech AS "Mech",
  deliv AS "Deliv",
  misha_fit AS "Misha",
  tier AS "Tier",
  next_step AS "Next step",
  file.link AS "Hub"
FROM "research/strc/hypotheses"
WHERE type = "hypothesis-hub"
SORT hypothesis_num ASC
```

### Supporting notes (ref — no hub, not in active register)

These hypotheses/models feed into the active register above but don't get their own hypothesis-hub folder. They're either engineering sub-layers, survey notes, or clinical strategy docs.

| Hypothesis / model                                 | Role                                         |
| -------------------------------------------------- | -------------------------------------------- |
| [[STRC RBM24 Regulatory Hypothesis]]               | Rolled into Strategy A                       |
| [[STRC Stereocilia Bundle Mechanics Model]]        | Supporting biophysics model                  |
| [[STRC AAV Vector Design]]                         | Mini-STRC implementation detail              |
| [[STRC B8 Enhancer Selection]]                     | AAV pipeline                                 |
| [[STRC Dual-Vector vs Single-Vector Transduction]] | Engineering layer                            |
| [[STRC Anti-AAV Immune Response Model]]            | Re-dosing constraint                         |
| [[STRC Electrostatic Analysis E1659A]]             | Pharmacochaperone input                      |
| [[Adult Treatment Window STRC]]                    | Constraint on all mechanisms                 |
| [[Alternative STRC Delivery Hypotheses]]           | Delivery-layer survey                        |
| [[Misha Compound-Het Therapy Stack Model]]         | Clinical plan: parallel h01+h03 stack        |

_All 18 numbered hypotheses (h01–h16, h26, h27) now live in `research/strc/hypotheses/h{N}/hub.md` and render via the Dataview block above. Manual table removed 2026-04-23 — edit hub frontmatter to change tier/score/next-step._

## S-tier (live)

```dataview
TABLE WITHOUT ID
  link(hypothesis_title) AS "Hypothesis",
  next_step AS "Next step",
  file.link AS "Hub"
FROM "research/strc/hypotheses"
WHERE type = "hypothesis-hub" AND tier = "S"
SORT hypothesis_num ASC
```

**Clinical plan for Misha** per [[Misha Compound-Het Therapy Stack Model]]: pursue #1 (PC for maternal E1659A) + #3 (AAV for paternal null) in parallel. h01 is the only monotherapy route to NORMAL (≤ 25 dB ABR); stack provides redundancy and lower drug/AAV burden. Phase 5 MD on #1 is the critical-path compute that produces the f_PC parameter.

## Kill list (D-tier, live)

```dataview
TABLE WITHOUT ID
  link(hypothesis_title) AS "Hypothesis",
  next_step AS "Reason killed",
  file.link AS "Hub"
FROM "research/strc/hypotheses"
WHERE type = "hypothesis-hub" AND tier = "D"
SORT hypothesis_num ASC
```

## Decision gates

- **→ D**: mechanism falsified OR Misha-fit=1 OR delivery impossible 10 yr.
- **→ C**: no near-term computational step; log external catalyst.
- **→ B**: evidence neutral, shallow engagement.
- **→ A**: real progress possible, not highest-leverage today.
- **→ S**: top 5 max, advances Misha directly, actionable next step.

## Literature-audit status (live)

Per [[feedback_literature_first]] rule + [[AGENTS]] §0c: before any computational proof, parameter-provenance must be verified against primary literature. Status is `lit_audit` frontmatter field on each hub. Detailed audit notes live in per-hub proof notes (see `[[STRC h{N} Parameter Provenance Audit YYYY-MM-DD]]`).

```dataview
TABLE WITHOUT ID
  hypothesis_num AS "#",
  link(hypothesis_title) AS "Hypothesis",
  tier AS "Tier",
  lit_audit AS "Audit",
  lit_audit_date AS "Date",
  file.link AS "Hub"
FROM "research/strc/hypotheses"
WHERE type = "hypothesis-hub"
SORT hypothesis_num ASC
```

## Connections

- `[part-of]` [[STRC Research Portal]]
- `[see-also]` [[STRC Hypothesis Ranking Log]]
- `[see-also]` [[STRC Mini-STRC Single-Vector Hypothesis]]
- `[see-also]` [[STRC Pharmacochaperone Virtual Screen E1659A]]
- `[see-also]` [[STRC Piezoelectric TM Bioelectronic Amplifier]]
- `[see-also]` [[STRC Calcium Oscillation Acoustic Therapy]]
- `[see-also]` [[STRC Synthetic Peptide Hydrogel HTC]]
- `[about]` [[Misha]]
