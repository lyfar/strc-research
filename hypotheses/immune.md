---
date: 2026-04-22
type: synthesis
title: Immune Response — Anti-AAV Antibodies and Re-dosing Strategy
tags:
- immune
- AAV
- neutralizing-antibodies
- seroprevalence
- re-dosing
- strc
- synthesis
- hypothesis
- legacy-hypothesis-sheet
status: active
priority: supporting
stage: computational-modeling
---

> Legacy STRC synthesis sheet. Canonical hypothesis state now lives in `research/strc/hypotheses/hXX-*/hub.md`; durable claims live as atomic notes in `notes/`.

## Core claim

AAV gene therapy is a one-shot treatment. After first administration, neutralizing antibodies (NAbs) form and block re-dosing. This is a major clinical problem for STRC therapy because (1) expression may decline over time, (2) children grow and may need repeat dosing, (3) first dose may fail (wrong serotype, subtherapeutic titer).

## Seroprevalence data

Pre-existing NAbs in general population (block gene therapy before first dose):
- AAV1: 67%
- AAV2: 72%
- AAV5: 40%
- AAV8: 38%
- AAV9: 47%

## Key questions

- What is seroprevalence in target age group (3–5 year old children)?
- Which AAV serotype has lowest pre-existing immunity for cochlear delivery?
- Can plasmapheresis or immunosuppression enable re-dosing?
- Does LNP-based delivery (no viral capsid) solve the immune memory problem?

## Strategic implications

LNP-based delivery (sonoporation hypothesis) has no AAV capsid → no immune memory → repeatable dosing. This is a key argument for the hybrid approach: AAV for initial dose, LNP for top-ups.

Anc80L65 (used by Iranfar 2026) shows ~80% OHC transduction but carries standard AAV immune constraints.

## Open questions

- NAb kinetics after cochlear AAV administration — does intracochlear route reduce systemic immunity?
- Minimum interval between doses with immunosuppression?
- Cross-reactivity between serotypes (Anc80L65 vs AAV9 etc.)?

## Recent Papers

### 2026-04-16 — [Gene Therapy for Infants and Children With Otoferlin-Related Auditory Neuropathy Spectrum Disorder](https://pubmed.ncbi.nlm.nih.gov/41979424/)
*Ear and Hearing · relevance: 7/10 · fda_path*

OTOF pediatric trials will generate critical real-world data on immune responses to cochlear AAV administration — including NAb formation kinetics, systemic vs. local immune responses, and whether intracochlear delivery reduces immunogenicity. This directly addresses our open question about NAb kinetics after cochlear AAV administration. Any immunosuppression protocols used in OTOF trials will establish precedent for managing immune barriers in subsequent hearing loss gene therapies.

## Connections

- `[part-of]` [[STRC]]
- [[STRC Anti-AAV Immune Response Model]]
- `[see-also]` [[STRC Hypothesis Ranking]]
