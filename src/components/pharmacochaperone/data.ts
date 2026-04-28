export const stats = [
  {
    label: 'Current rank',
    value: 'A',
    note: 'Hold pending wetlab RFQ feasibility (A_hold_wetlab_rfq_feasibility).',
  },
  {
    label: 'APBS pocket Δ',
    value: '−2.79 kcal/mol',
    note: 'Matched ensemble n=20: MUT +5.99±1.37 vs WT +1.46±1.55 kT/e (Welch p=6.9e-12, Cohen d=3.10).',
  },
  {
    label: 'RWM permeability',
    value: '7.6× TMPA',
    note: 'Phase 8h-lite Salt 2001 Stokes-Einstein model: ~60% applied conc at basal turn @ 90 min.',
  },
  {
    label: 'Lead committee',
    value: 'α-robust',
    note: 'Phase 9x: lead ranking stable across kinetic-axis weight α=0.5–2.0 on Boltz/Vina/tauRAMD/ADMET.',
  },
] as const;

export const claims = [
  'h01 is an A-tier computational hypothesis after Phase 8h-lite paper-claim closures (2026-04-26) and Phase 9x α-sensitivity check.',
  'E1659A creates a measurable formal-anion preference at the K1141 pocket: Δ = −2.79 ± 0.28 kcal/mol (matched n=20 APBS, Welch p=6.9e-12, Cohen d=3.10).',
  'K1141 pocket meets all 4 Halgren druggability gates (V=1145 Å³, phobic 0.61, V_pocket/V_lig 2.75).',
  'A 1066-1216 STRC fragment is a plausible first biophysics reagent for K1141-pocket assays (AF3 fragment ladder; ring RMSD 0.97 Å).',
  'A weak-anion acyl-sulfonamide v5.3 chemistry set is tractable: 2 lead candidates plus 1 optional comparator are ready to RFQ.',
  '1 mM intratympanic delivery is plausible in the 8-compartment cochlear PBPK model (94–97% time above 1 µM in HC_ER, no cytosol safety breach).',
] as const;

export const nonClaims = [
  'Not a lead-candidate nomination; this is RFQ-scoping only.',
  'Not a clinical-readiness or IND claim.',
  'No physical compound has been synthesized.',
  'No STRC 1066-1216 WT/E1659A reagent has been expressed.',
  'No DSF / nanoDSF / SPR / BLI / MST or cell-rescue data exist.',
  'Off-target electrostatic selectivity remains INCONCLUSIVE LITE; APBS on enclosed pockets per off-target is the heavy planned step.',
  'Phase 8h-lite is light-compute paper-claim closure, not production MD/MM-PBSA replacement.',
] as const;

export const reviewerProfiles = [
  'Otology / inner-ear small-molecule delivery',
  'Stereocilin / STRC structural biology',
  'CADD / docking + MD methodology',
  'Acyl-sulfonamide medicinal chemistry',
] as const;

export const packContents = [
  'Phase 5k matched-ensemble APBS gate-closing result',
  'Phase 8h-lite 5 light proofs (druggability, RWM, off-target, pose stability, ototoxin Tanimoto)',
  'Phase 9x lead committee α-sensitivity',
  'Phase 9b Min-Relax MM-GBSA orthogonal smoke',
  'Wetlab Handoff Packet 2026-04-27 (compound sheet + cascade + gap matrix)',
  'Reviewer objection table',
] as const;

export const methodStack = [
  {
    method: 'Matched WT/MUT structural ensembles',
    compute: 'Full-length STRC MD snapshots, WT and E1659A, then matched APBS pocket electrostatics.',
    scale: '20 WT frames + 20 mutant frames',
    key: 'Mutant pocket formal-anion preference = −2.79 ± 0.28 kcal/mol; p=6.9e-12.',
  },
  {
    method: 'Fragment and pocket validation',
    compute: 'AF3 fragment ladder plus K1141 local-geometry checks against the MD-relaxed parent structure.',
    scale: '8 AF3 jobs, 40 CIFs, 40 confidence files',
    key: '151-aa E1659A fragment kept K1141-ring RMSD at 0.97 Å; pocket pLDDT 78.7.',
  },
  {
    method: 'Chemistry library design',
    compute: 'RDKit triage of weak-anion acyl-sulfonamide bioisostere candidates.',
    scale: '12 designed SMILES, 6 passed filters',
    key: 'Top-3 v5.3 candidates selected for Boltz-2, Vina, ADMET, and tauRAMD.',
  },
  {
    method: 'Boltz-2 + Vina consensus',
    compute: 'Protein-ligand cofolding and mutant-ensemble redocking of the v5.3 top-3.',
    scale: '5 Boltz diffusion seeds × 3 ligands; 20 MD snapshots × 3 ligands in Vina',
    key: 'Both methods put 1-indanyl_SO2Me_-Cl on the affinity track; all 3 passed ipTM ≥ 0.50.',
  },
  {
    method: 'Residence-time / tauRAMD',
    compute: 'Biased-MD ligand escape simulations on STRC and off-target TRPM4.',
    scale: '25/25 v5.2 STRC unbound; 6/6 TRPM4 unbound; 60/60 v5.3 STRC unbound',
    key: 'Selectivity rescue did not clear the >5× gate; v5.3 cross-target bound ≤2.72×.',
  },
  {
    method: 'Off-target electrostatics',
    compute: 'APBS pocket-potential resampling across STRC and cochlear ion-channel off-targets.',
    scale: 'STRC, TRPM4, TMEM16A, KCNQ4, Cx50 panels',
    key: 'STRC pocket is electropositive, but dominant off-target pocket ceilings remain higher.',
  },
  {
    method: 'Cochlear exposure model',
    compute: '8-compartment PBPK ODE with intratympanic dose sweep.',
    scale: '100 µM / 1 mM / 10 mM × fast/default/slow apical uptake',
    key: '1 mM IT gives 94–97% time above 1 µM in HC_ER without cytosol safety breach in the model.',
  },
  {
    method: 'Endpoint free-energy smoke',
    compute: 'Repaired matched Amber WT/E1659A complexes, restrained min-relax frames, MMPBSA.py GB scoring.',
    scale: '6 frames WT + 6 frames E1659A',
    key: 'E1659A − WT ΔΔTOTAL = −2.35 kcal/mol; supportive smoke, not production holo MD.',
  },
] as const;

export const phaseEvidence = [
  {
    step: '01',
    layer: 'Pocket premise',
    evidence: 'K1141 pocket profile',
    measurement: 'Volume 1145 Å³; hydrophobic fraction 0.61; pocket-to-ligand volume ratio 2.75.',
    interpretation: 'All 4 Halgren druggability gates pass (Phase 8h-lite #1).',
    state: 'Supports mechanism',
  },
  {
    step: '02',
    layer: 'Mutant specificity',
    evidence: 'Matched-ensemble APBS',
    measurement: 'WT +1.46 ± 1.55 kT/e; E1659A +5.99 ± 1.37; Δ +4.53 ± 0.46 kT/e (p=6.9e-12).',
    interpretation: 'E1659A creates a statistically strong anion-favouring electrostatic environment.',
    state: 'Strong computational support',
  },
  {
    step: '03',
    layer: 'Construct feasibility',
    evidence: 'AF3 K1141 fragment ladder',
    measurement: '8 AlphaFold Server jobs; best operational fragment = 151_e1659a, ring RMSD 0.97 Å.',
    interpretation: 'A 1066-1216 fragment is a plausible first biophysics reagent for K1141-pocket assays.',
    state: 'Ready for expression RFQ',
  },
  {
    step: '04',
    layer: 'v5.2 kinetic screen',
    evidence: 'Production tauRAMD on v5.2 shortlist',
    measurement: '5 ligands × 5 replicas; 25/25 unbound; within-STRC residence spread 1.7×.',
    interpretation: 'Residence-time differences exist, but same-head-class kinetic selectivity is too small.',
    state: 'Constrains claim',
  },
  {
    step: '05',
    layer: 'Off-target risk',
    evidence: 'TRPM4 cross-target tauRAMD',
    measurement: '2 ligands × 3 replicas; TRPM4/STRC residence ratios 1.08–1.52×.',
    interpretation: 'Kinetic-selectivity rescue failed for v5.2 by at least a 3× margin versus the >5× gate.',
    state: 'Negative but useful',
  },
  {
    step: '06',
    layer: 'Exposure model',
    evidence: 'PBPK dose sweep',
    measurement: '1 mM IT: 93.8–96.9% of 72 h with HC_ER ≥ 1 µM; 0% HC_CYT ≥ 50 µM.',
    interpretation: 'Exposure is not the main blocker if 1 mM intratympanic dosing is acceptable.',
    state: 'Supports delivery feasibility',
  },
  {
    step: '07',
    layer: 'v5.3 chemistry',
    evidence: 'Acyl-sulfonamide library',
    measurement: '12 candidates designed; 6/12 passed Lipinski + TPSA ≤ 100 + mono-anion gates.',
    interpretation: 'Weak-anion head-group pivot produced a tractable next-gen chemistry set.',
    state: 'Supports next experiment',
  },
  {
    step: '08',
    layer: 'Affinity ranking',
    evidence: 'Boltz-2 + Vina v5.3 consensus',
    measurement: 'Boltz ipTM rank-1 0.645 for 1-indanyl_-Cl; Vina mean rank-1 −0.67 kcal/mol.',
    interpretation: '1-indanyl_acylsulfonamide_SO2Me_-Cl is the affinity-track compound.',
    state: 'Lead-track evidence',
  },
  {
    step: '09',
    layer: 'Conformational robustness',
    evidence: 'v5.3 mutant-ensemble redock',
    measurement: '3 ligands × 20 E1659A MD snapshots; 1-indanyl_-Cl rank-1 by mean and median.',
    interpretation: 'Affinity signal is not just a single static receptor pose.',
    state: 'Supports affinity track',
  },
  {
    step: '10',
    layer: 'Developability triage',
    evidence: 'ADMET-AI top-3',
    measurement: '2/3 clean on 10 endpoints; 1-indanyl has one borderline CYP3A4 flag at 90.4 vs 90.0 gate.',
    interpretation: 'Use a two-track lead committee: affinity lead plus ADMET-clean adamantyl backup.',
    state: 'Actionable tradeoff',
  },
  {
    step: '11',
    layer: 'v5.3 residence check',
    evidence: 'STRC tauRAMD n=20 per ligand',
    measurement: '60/60 unbound; tau means 21.59, 17.46, 15.90 ps; spread 1.358×.',
    interpretation: 'Class does not generate enough kinetic spread for a selectivity claim, but remains useful for target-engagement testing.',
    state: 'Narrows wetlab claim',
  },
  {
    step: '12',
    layer: 'Orthogonal score smoke',
    evidence: 'Min-relax MM-GBSA',
    measurement: 'WT −15.84; E1659A −18.18; Δ −2.35 kcal/mol.',
    interpretation: 'An APBS-independent endpoint score preserves the mutant-favouring sign.',
    state: 'Supportive smoke',
  },
  {
    step: '13',
    layer: 'α-sensitivity',
    evidence: 'Phase 9x lead committee',
    measurement: 'Rank aggregation stable for kinetic-axis weight α=0.5–2.0 across Boltz/Vina/tauRAMD/ADMET.',
    interpretation: 'The two-track committee decision is not fragile to weight choice.',
    state: 'Robust committee',
  },
] as const;

export const lightProofs = [
  {
    n: '1',
    title: 'K1141 pocket druggability',
    method: 'Halgren SiteMap-style grid burial estimator (Halgren 2009 spec).',
    result: 'V = 1145 Å³, phobic 0.609, V_pocket/V_lig 2.75. All 4 druggability gates ✅. Pocket lining: W1612, F1646, F1169, W1652 (π-stacking shell); K1141 = 26% of philic shell.',
    verdict: 'PASS',
    verdictTone: 'pass' as const,
  },
  {
    n: '2',
    title: 'RWM permeability',
    method: 'Salt 2001 P_TMPA baseline + Stokes-Einstein MW^(−1/3) + Avdeef logP/TPSA + Henderson-Hasselbalch f_neutral.',
    result: 'Lead P = 1.45×10⁻⁷ cm/s = 7.6× TMPA. Inferred basal-turn concentration ~60% applied @ 90 min vs 8% for TMPA. CONHOH pKa sensitivity 7.2–7.6× (robust).',
    verdict: 'PASS — deliv 3→4',
    verdictTone: 'pass' as const,
  },
  {
    n: '3',
    title: 'Off-target electrostatic selectivity',
    method: 'Two proxies: Vina-pose-anchored count + receptor-wide K/R cluster scan.',
    result: 'Vina parked CONHO⁻ at 7.5 Å from K1141 (Gasteiger-neutral artefact). Receptor-wide scan returns large K/R clusters on all 4 channels because their voltage-sensors / selectivity filters are K/R-rich — not enclosed druggable pockets.',
    verdict: 'INCONCLUSIVE LITE',
    verdictTone: 'caution' as const,
  },
  {
    n: '4',
    title: 'Pose ensemble stability',
    method: 'Distance K1141 NZ ↔ lead CONHO⁻ across all 20 Phase 5d MD snapshots (LIE substitute).',
    result: 'K1141 NZ per-axis SD 1.72 / 1.15 / 1.10 Å (rigid). Lead CONHO⁻ ↔ K1141 NZ mean 5.02 ± 0.57 Å, range 3.93–6.18 Å — Coulomb-attraction range across entire ensemble.',
    verdict: 'PASS — mech 4 holds',
    verdictTone: 'pass' as const,
  },
  {
    n: '5',
    title: 'Tanimoto vs known cochlear ototoxins',
    method: 'RDKit Morgan FP r=2 2048 bits; lead vs 12-compound Schacht 2008 panel (aminoglycosides, loop diuretics, platinums, salicylates, macrolide, quinoline, glycopeptide).',
    result: 'Closest match Tanimoto 0.127 (aspirin). No structural class motif overlap with known ototoxins.',
    verdict: 'PASS — no class flag',
    verdictTone: 'pass' as const,
  },
] as const;

export const compounds = [
  {
    priority: '1',
    name: '1-indanyl_acylsulfonamide_SO2Me_-Cl',
    track: 'Affinity lead',
    evidence: 'Boltz-2 ipTM 0.645; Vina mean −0.67; mut-ensemble rank-1 by mean and median.',
    risk: 'One borderline CYP3A4 ADMET-AI flag (90.4 vs 90.0 gate); not the residence-time winner.',
  },
  {
    priority: '2',
    name: 'adamantyl_acylsulfonamide_SO2Me_-Cl',
    track: 'Practical backup',
    evidence: 'ADMET-clean on all 10 endpoints; tightest tauRAMD SEM; α-sensitivity-robust.',
    risk: 'Mutant-ensemble affinity is weaker and more best-pose-dependent.',
  },
  {
    priority: '3',
    name: 'adamantyl_acylsulfonamide_SO2Me_-CF3',
    track: 'Optional residence comparator',
    evidence: 'ADMET-clean; tauRAMD mean rank-1 (in-class).',
    risk: 'Wider tau outliers and weaker affinity signal.',
  },
] as const;

export const assayCascade = [
  {
    gate: '0',
    question: 'Can the reagents be made?',
    readout: 'Custom synthesis QC (HPLC, LC-MS, identity); STRC 1066-1216 WT/E1659A expression (HEK293 or Sf9, His6/AviTag), SDS-PAGE, SEC, identity.',
    pass: 'Usable monodisperse WT/E1659A fragment + identity-confirmed compounds in hand.',
    stop: 'Fragment mis-folds, aggregates, or cannot be purified; or compounds fail identity/purity.',
  },
  {
    gate: '1',
    question: 'Does ligand stabilize the fold?',
    readout: 'DSF or nanoDSF, dose response, vehicle-matched.',
    pass: 'Reproducible thermal shift vs vehicle, ideally larger on E1659A than WT.',
    stop: 'No reproducible shift, or WT shift ≥ E1659A shift.',
  },
  {
    gate: '2',
    question: 'Is binding direct?',
    readout: 'SPR / BLI / MST; DMSO-matched; reference surfaces; immobilization-artefact controls.',
    pass: 'Dose-dependent Kd, koff or stoichiometry without aggregation-shaped curves.',
    stop: 'Aggregation curves, no concentration response, or vehicle binding.',
  },
  {
    gate: '2.5',
    question: 'Is the signal clean?',
    readout: 'Aggregation / reactivity panel; hERG, TRPM4, TMEM16A, KCNQ4, BK, Cx50 counterscreens; COX if NSAID-like risk remains.',
    pass: 'No meaningful activity at relevant multiples of measured Kd.',
    stop: 'Off-target activity within or below the on-target Kd window.',
  },
  {
    gate: '3',
    question: 'Does it rescue cells?',
    readout: 'HEK293 full-length E1659A STRC trafficking / mature-glycoform readout.',
    pass: 'Dose-dependent E1659A rescue after direct biophysics passes.',
    stop: 'No rescue at non-toxic doses, or rescue paired with toxicity.',
  },
] as const;

export const packetSections = [
  { item: 'Decision memo',     include: 'One-page H01 aim, allele, K1141 pocket, ask, kill/pass criteria.', exclude: 'IND language or "lead candidate" framing.' },
  { item: 'Construct sheet',   include: 'STRC 1066-1216 K1141 fragment WT + E1659A; His6/AviTag options; HEK293/Sf9 expression feasibility.', exclude: 'Full-length STRC protein production as first reagent.' },
  { item: 'Compound sheet',    include: '2 required molecules + 1 optional comparator with SMILES, computed properties, QC asks.', exclude: 'Large exploratory library synthesis.' },
  { item: 'Assay cascade',     include: 'Gate 0 protein QC → Gate 1 DSF/nanoDSF → Gate 2 SPR/BLI/MST → Gate 2.5 counterscreens → Gate 3 cell rescue.', exclude: 'CryoEM or animal PK before direct binding.' },
  { item: 'Evidence appendix', include: 'APBS sign, Boltz/Vina ranking, tauRAMD class verdict, ADMET-AI, min-relax MM-GBSA smoke with caveats.', exclude: 'Production MD as wetlab prerequisite.' },
  { item: 'RFQ artifacts',     include: 'compound_sheet.csv, gap_matrix.csv, rfq_email.md.', exclude: 'Vendor cost claims before live RFQs.' },
] as const;

export const objectionRows = [
  {
    objection: 'Vina ranks WT and mutant similarly.',
    answer: 'Diagnosed as a Gasteiger-neutral force-field artefact (Phase 5k, Phase 8h-lite #3): docking does not see the formal-anion charge.',
    risk: 'Reviewer may take Vina ranks at face value.',
    nextCheck: 'Lead with matched-ensemble APBS (p=6.9e-12, d=3.10), not Vina. Wet DSF/nanoDSF closes the question.',
  },
  {
    objection: 'APBS Δ is one ensemble; could be MD-trajectory artefact.',
    answer: 'Matched n=20 WT vs n=20 E1659A frames from the same MD protocol; Welch t=−13.7, Cohen d=3.10.',
    risk: 'A new MD seed could shift the Δ.',
    nextCheck: 'Phase 5d/5k MD inputs, scripts, and seeds are versioned in the repo; reviewer can re-run. DSF Δ on E1659A vs WT is the definitive falsifier.',
  },
  {
    objection: 'tauRAMD did not clear the >5× selectivity gate.',
    answer: 'Acknowledged. Class kinetic spread is small; we explicitly do NOT claim kinetic selectivity.',
    risk: 'Reviewer may read this as h01 disqualified.',
    nextCheck: 'Reframe ask as target engagement (DSF + SPR), not kinetic discrimination. Counterscreens close selectivity wet, not in silico.',
  },
  {
    objection: 'Off-target electrostatic selectivity is INCONCLUSIVE LITE.',
    answer: 'Phase 8h-lite #3 confirmed: receptor-wide K/R scans return false-positives because off-target voltage-sensors are K/R-rich but not enclosed druggable pockets.',
    risk: 'Reviewer cannot rule out off-target binding from in silico.',
    nextCheck: 'Per-target APBS on enclosed pockets is the heavy step; Gate 2.5 wet counterscreen panel (hERG / TRPM4 / TMEM16A / KCNQ4 / BK / Cx50) closes the question empirically.',
  },
  {
    objection: 'Min-relax MM-GBSA is not production holo MD.',
    answer: 'Correct. Phase 9b smoke is supportive sign-check (Δ −2.35 kcal/mol), not paper-grade ΔG.',
    risk: 'Reviewer may discount because it is not full holo MD.',
    nextCheck: 'Production holo MD is paper-grade upgrade, not wetlab-blocker. Phase 9b SKELETONs scoped if needed.',
  },
  {
    objection: 'AF3 fragment may not refold like the full-length parent in vitro.',
    answer: 'AF3 ladder shows 151_e1659a ring RMSD 0.97 Å vs MD-relaxed parent and pocket pLDDT 78.7 — model agreement only.',
    risk: 'In-vitro refolding can fail at any step.',
    nextCheck: 'Gate 0: SDS-PAGE / SEC / identity QC on actual expressed fragment. RFQ scopes both HEK293 and Sf9 routes.',
  },
  {
    objection: 'Pharmacochaperones for extracellular targets are unusual.',
    answer: 'Acknowledged in mechanism flag (deliv 3, not 4). E1659A is a cysteine-loss variant predicted to cause ER mis-folding before secretion; rescue would happen at the ER stage, not extracellularly.',
    risk: 'Reviewer may anchor on extracellular-mature framing.',
    nextCheck: 'Gate 3 cell rescue uses trafficking / mature-glycoform readout, not extracellular function.',
  },
  {
    objection: 'PBPK 1 mM IT model has not been calibrated.',
    answer: '8-compartment ODE uses Salt 2001 RWM permeability + literature scala-tympani clearance + Avdeef logP/pKa.',
    risk: 'Real intracochlear distribution may diverge.',
    nextCheck: 'Phase 8h-lite #2 RWM Stokes-Einstein scaling + lead 7.6× TMPA permeability is an independent supporting check.',
  },
  {
    objection: 'CYP3A4 ADMET flag on the affinity lead.',
    answer: 'Borderline: 90.4 vs 90.0 gate. The two-track committee includes the ADMET-clean adamantyl backup specifically to derisk this.',
    risk: 'Real CYP3A4 inhibition would create DDI exposure window.',
    nextCheck: 'In-vitro CYP panel after first synthesis; Gate 2.5.',
  },
  {
    objection: 'No physical compound and no STRC reagent yet.',
    answer: 'Acknowledged; the entire ask is RFQ-scoping for synthesis + fragment expression + DSF/SPR.',
    risk: 'Empirical work is the only thing that can move h01.',
    nextCheck: 'Wetlab Handoff Packet 2026-04-27 is RFQ-ready.',
  },
] as const;

export const expansionCriteria = [
  'Missense allele, not null or large deletion.',
  'Residual protein is made and plausibly foldable.',
  'A local pocket or allosteric stabilizer site can be modeled.',
  'Target engagement can be measured before cell rescue.',
  'The stabilizer does not permanently block the protein function it is meant to restore.',
] as const;

export const sections = [
  { id: 'ask',         label: 'Exact ask' },
  { id: 'claims',      label: 'Claims & non-claims' },
  { id: 'methods',     label: 'What was computed' },
  { id: 'evidence',    label: 'Phase evidence' },
  { id: 'lite',        label: 'Phase 8h-lite proofs' },
  { id: 'compounds',   label: 'Compound priority' },
  { id: 'cascade',     label: 'Assay cascade' },
  { id: 'packet',      label: 'RFQ packet' },
  { id: 'objections',  label: 'Objection table' },
  { id: 'platform',    label: 'Platform criteria' },
  { id: 'boundary',    label: 'Claim boundary' },
] as const;
