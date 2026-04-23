---
date: 2026-04-23
type: source
tags: [strc, h05, paper, met-channel, outer-hair-cell, calcium, conductance, mechanotransduction]
hypotheses: [h05]
title: "A large-conductance calcium-selective mechanotransducer channel in mammalian cochlear hair cells"
authors: ["Beurg M", "Evans MG", "Hackney CM", "Fettiplace R"]
journal: "Journal of Neuroscience"
year: 2006
volume: 26
pages: "10992-11000"
pubmed_id: "17065441"
doi: "10.1523/JNEUROSCI.2188-06.2006"
---

## What they found

Patch-clamp recording from apical and basal OHCs in neonatal rat cochlea. Single-channel conductance of the MET channel increases from 145 pS (apical, low CF) to 210 pS (basal, high CF). Channels are highly Ca²⁺-selective (PCa/PNa ≈ 5–7). Ca²⁺ fraction of the total MET current is ~15% under physiological ionic conditions.

## Numbers that matter

- Single MET channel conductance: **145–210 pS** (apical to basal gradient)
- Model value `g_MET = 150 pS` is correct for apical OHC
- Ca²⁺ fraction of MET current `f_Ca = 0.15`: **confirmed** (Beurg 2006 Table 2)
- Channel open probability at rest: ~0.1–0.15 (before stimulation)
- Activation time constant: sub-millisecond
- Channel is mechanically gated via tip-link tension

## Fit to h05

Confirms both `g_MET = 150e-12 S` and `f_Ca = 0.15` in the RBM24 ODE. These are the two load-bearing MET parameters in the Ca²⁺ influx term. Both values are well-supported.

## Connections

- [[STRC Calcium Oscillation Acoustic Therapy]] — g_MET and f_Ca
- `[see-also]` [[2019-fettiplace-kim-met-channel-review]] — broader review
- `[part-of]` [[calcium-oscillation]] (literature-params topic)
