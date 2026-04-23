---
title: "Actin-bound structures of Wiskott–Aldrich syndrome protein (WASP)-homology domain 2 and the implications for filament assembly"
authors: ["Chereau D", "Kerff F", "Graceffa P", "Grabarek Z", "Langsetmo K", "Dominguez R"]
journal: "Proceedings of the National Academy of Sciences USA"
year: 2005
volume: 102
pages: "16644–16649"
pubmed_id: "16275905"
pmc: "PMC1283820"
doi: "10.1073/pnas.0507021102"
date: 2026-04-23
type: source
tags: [strc, paper, actin, wh2, wasp, biophysics, h09]
hypotheses: [h09]
---

## What they found

Crystal structures of WH2 domains from WASP, WAVE2, and WIP in complex with G-actin (ternary
complex with DNase I as crystallization cofactor). Parallel ITC biochemistry comparing WH2 vs
thymosin-β (Tβ) binding to actin. First structural evidence that short WH2 domains (~17 aa)
can coexist with intersubunit contacts in F-actin; the canonical G-actin binding site (barbed-end
groove between subdomains 1 and 3) is partly accessible on the filament surface.

## Numbers that matter

**WH2 G-actin binding affinities — ITC, 25°C, G-buffer (2 mM Tris pH 7.5, 0.2 mM CaCl₂,
0.2 mM ATP):**
- Exact Kd values are in Table 2 (SI, not open-access). Narrative in text only gives relatives:
  - WAVE WH2: strongest binder (5× higher affinity than WASP WH2)
  - WASP WH2: weakest binder in series
  - All WH2 domains: ~10× higher affinity than Tβ domains
- Since Tβ4 Kd = 1 μM (Husson 2010, Xue 2014), this implies WH2 Kd range ~50–200 nM
- WASP WA domain (WH2 + C region together): Kd = 0.6 μM from Rohatgi 2000 (cited)
- Isolated WASP WH2 (V domain alone, from Padrick/Kim 2011 citing this paper): Kd ≈ 3.1 μM

**CRITICAL NOTE on "100 nM" in model:** The ~10× boost vs Tβ4 is for construct-dependent
WH2 variants. The isolated WASP WH2 V-domain alone is ~3 μM, not 100 nM. The 100 nM figure
likely applies to long WH2 constructs (WIP-type) or in-context domains with favorable
electrostatics (WAVE with Arg-Arg in LKKT gives extra salt bridges). Model parameter
WH2_KD_GACTIN_M = 200 nM is near the favorable end; range is 50 nM–3 μM depending on construct.

**WH2 × F-actin side-binding — NOT MEASURED in this paper:**
- The paper shows WH2 binding site (barbed-end groove, subdomains 1/3) is ALSO the site of
  longitudinal actin–actin contacts in the filament long-pitch helix
- Short WH2 domains (~17 aa) can "coexist with intersubunit contacts in F-actin" structurally,
  but this means nucleation (WH2 templating monomers along a strand), not side-binding of
  preformed filaments
- No co-sedimentation or F-actin binding Kd measured

**Tβ4 × F-actin:** weak cooperative binding, Kd = 5–10 mM (cited from Husson 2010)
This is 3 orders of magnitude weaker than G-actin binding — the C-terminal α-helix in Tβ4
interferes with filament contacts, but the same region being absent in WH2 doesn't mean
WH2 side-binds better; it just means WH2 doesn't lock monomer as completely.

## Why WH2 × F-actin side-binding is structurally disfavored

The N-terminal amphipathic helix of WH2 binds in the cleft between actin subdomains 1 and 3.
In the filament (Holmes model), this same cleft is occupied by longitudinal actin–actin contacts
along the long-pitch helix — it is partially buried. Residual surface accessibility does not
constitute a high-affinity independent binding site. No evidence that WH2 can bind to the side
of a preformed filament with affinity < 1 mM. The model's WH2_KD_FACTIN_M = 5 μM is optimistic;
the real value is likely 1–100 mM or unmeasurable.

## Connections

- `[primary source for]` [[STRC H09 WH2 F-actin Bundling Hypothesis]] — provides structural
  rationale and Kd reference range for G-actin binding; establishes F-actin side-binding as
  unmeasured and structurally disfavored
- `[informs]` [[Hydrogel Phase 4d F-actin Bundling Model]] — WH2_KD_GACTIN_M defensible at
  50–200 nM for favorable WH2 constructs; WH2_KD_FACTIN_M = 5 μM is speculative/optimistic
- `[structural basis for]` [[WH2 Actin Binding Mechanism]] — barbed-end groove binding
