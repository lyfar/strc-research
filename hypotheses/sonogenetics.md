---
title: Sonogenetics — Sound-Activated Gene Expression
status: active
stage: computational-modeling
priority: secondary
tags: [sonogenetics, mechanosensitive, NFAT, calcium, promoter, ODE]
---

## Core claim

Outer hair cells (OHCs) are natural mechanosensors. Sound activates MET channels → Ca²⁺ influx → calcineurin/NFAT pathway → NFAT-responsive promoter drives therapeutic gene expression. Therapeutic dose only when sound is present — built-in safety switch.

## Cascade

Sound → stereocilia deflection → MET channels open → Ca²⁺ enters (134 channels/cell, 15% Ca²⁺ permeability) → [Ca²⁺] rises from 70 nM → 500–900 nM in 0.05 pL apical volume → calcineurin activated (Hill n=4) → NFAT dephosphorylated → nuclear translocation → NFAT-responsive promoter → transgene expression.

## Current evidence

- ODE model built and calibrated: Ca²⁺ dynamics, calcineurin kinetics, NFAT translocation, promoter activation
- Bifurcation analysis: system shows bistable behavior — threshold response to sound amplitude
- Sensitivity analysis: NFAT translocation rate is limiting step, not Ca²⁺ influx
- B8 synthetic enhancer (ARBITER panel) found to drive OHC-specific expression with zero ectopic

## Key parameters

- MET channels per OHC: 134
- Ca²⁺ permeability: 15%
- Apical volume: 0.05 pL
- Resting [Ca²⁺]: 70 nM
- Peak [Ca²⁺] during sound: 500–900 nM
- Hill coefficient for calcineurin activation: n=4

## Open questions

- Has NFAT pathway been validated as functional in mature OHCs in vivo?
- What sound intensity threshold triggers meaningful NFAT activation?
- Does repeated sound stimulation cause NFAT fatigue or desensitization?
- Can NFAT-responsive promoter drive therapeutic protein levels (not just reporter)?

## Relationship to other hypotheses

Sonogenetics addresses both delivery (sound-triggered uptake) and expression control. If validated, could be combined with mini-STRC (payload) and LNP delivery (vector) for a non-surgical, repeatable therapeutic system.

## Recent Papers

### 2026-04-20 — [Multicolored Mechanoluminescent Nanoparticles for Functional Sono-Optogenetics](https://pubmed.ncbi.nlm.nih.gov/41974592/)
*J Am Chem Soc · relevance: 7/10 · lateral*

Developed nanoparticles that convert focused ultrasound into tunable light emission (461-592 nm) via mechanoluminescence, successfully activating ChR2, eOPN3, and ChRmine in deep tissue. Directly enables a hybrid sono-optogenetic approach for the cochlea: AAV delivers a channelrhodopsin to hair cells, mechanoluminescent nanoparticles in perilymph convert sound vibrations or external FUS into local light, creating a sound-responsive prosthetic system. The red-shifted emission option (592 nm) would minimize cochlear phototoxicity. Key question: do cochlear acoustic pressures reach the ~MPa thresholds needed for mechanoluminescence activation?
