export const stats = [
  {
    label: 'Current rank',
    value: 'S',
    note: 'Promoted after Phase 4m closed the conservative avidity gate.',
  },
  {
    label: 'Gate-closing result',
    value: '1.48 mM',
    note: 'Phase 4m F-actin-adjusted WH2 Kd soft floor used in the conservative avidity model.',
  },
  {
    label: 'Pass probability',
    value: '80.3%',
    note: 'Phase 4n ideal Monte Carlo pass probability; 59.0% under a 30x entropic haircut.',
  },
  {
    label: 'Dose window',
    value: '0.13-1.32 mg',
    note: 'Phase 4e modeled lower-to-upper ototopical dosing window before the toxic floor.',
  },
] as const;

export const claims = [
  'h09 is an S-tier computational hypothesis after Phase 4m Kd Mono Soft Floor Result (2026-04-27).',
  'The central computational gate closed is Kd_mono soft-floor / avidity plausibility: predicted h09 Kd_F = 1.48 mM after Tbeta4 G-to-F scaling; conservative avidity rerun gives Kd_eff = 1.30 uM at C_eff = 50 mM, N = 3.',
  'The therapy concept is local peptide-material rescue, not gene therapy.',
  'The intended population is STRC-loss hearing loss with residual outer hair cell function (incl. Misha).',
  'A Boltz-2 control panel is supportive (not definitive): h09 stays above scrambled / alanine / poly-GS negatives, but Boltz does not reproduce the full ITC rank order.',
] as const;

export const nonClaims = [
  'Not a cure claim.',
  'Not a clinical-readiness claim.',
  'Not a proof of cochlear delivery.',
  'Not a proof of absolute WH2 x F-actin Kd.',
  'Not a universal hearing-loss treatment.',
  'Not a claim that AF3 or Boltz predicts Kd directly.',
  'Not a claim that Phase 4k passed; Phase 4k is an explicit geometry caution.',
] as const;

export const reviewerProfiles = [
  'Inner-ear mechanobiology / stereocilia',
  'Biomaterials peptide hydrogel',
  'Cochlear delivery / translational otology',
  'Actin-binding motif biophysics',
] as const;

export const packContents = [
  'Phase 4m gate-closing result',
  'Phase 4k geometry caution',
  'Phase 4p Boltz-2 control panel',
  'Precedent table with literature anchors',
  'Reviewer objection table',
  '5-step minimum feasibility package',
] as const;

export const mechanismStack = [
  {
    step: '01',
    layer: 'Therapeutic shape',
    evidence: 'Synthetic horizontal top connector',
    measurement: 'A self-assembling peptide hydrogel is meant to replace STRC mechanical coupling, not STRC sequence.',
    interpretation: 'The hypothesis is mutation-agnostic for STRC-null patients if bundle-scale coupling can be restored.',
    state: 'Mechanism',
  },
  {
    step: '02',
    layer: 'Construct',
    evidence: 'PEPTIDE_TAIL91',
    measurement: 'WH2 + GSGSG + RADA16 + GSGSG + STRC aa 1620-1710; 134 aa total.',
    interpretation: 'One peptide combines actin engagement, self-assembly, and TMEM145-facing STRC tail mimicry.',
    state: 'Design lock',
  },
  {
    step: '03',
    layer: 'TMEM145 axis',
    evidence: 'Tail91 epitope retool',
    measurement: 'TMEM145 GOLD-zone complex ipTM 0.68 after switching to the correct STRC tail epitope.',
    interpretation: 'The TMEM145 side moved from a generic tail premise to a specific modeled contact surface.',
    state: 'Strong support',
  },
  {
    step: '04',
    layer: 'Actin axis',
    evidence: 'WH2 x G-actin and tail91 actin contact',
    measurement: 'Actin-axis ipTM 0.51 in Phase 3 / 4b family runs; Phase 4m WH2 calibration predicts Kd_F 1.48 mM.',
    interpretation: 'Single-site F-actin binding is weak, but multivalent presentation can make it useful.',
    state: 'Avidity-gated',
  },
  {
    step: '05',
    layer: 'Bundle mechanics',
    evidence: 'Coupling fraction model',
    measurement: 'Original mechanics model estimated about 22 dB rescue at f=0.3 and about 30 dB at f=0.6.',
    interpretation: 'The biological target is partial mechanical coupling rescue, not full restoration of native STRC architecture.',
    state: 'Model target',
  },
] as const;

export const evidenceStack = [
  {
    phase: 'Phase 1',
    result: 'Self-assembly screen',
    evidence: '5/6 candidate scaffolds passed the initial assembly filter; RADA16-WH2-native became the top track.',
    decision: 'Advance self-assembling peptide scaffold rather than soluble protein replacement.',
  },
  {
    phase: 'Phase 3',
    result: 'Tail91 retool',
    evidence: 'Correct-epitope STRC tail91 reached TMEM145 ipTM 0.68 and actin ipTM 0.51.',
    decision: 'Lock the 134-aa PEPTIDE_TAIL91 design for Phase 4 stress testing.',
  },
  {
    phase: 'Phase 4b v2',
    result: 'Developability scrub',
    evidence: 'Cys11/52->Ser candidate held TMEM145 median ipTM 0.55 and actin median ipTM 0.51 over n=5 each.',
    decision: 'Production candidate can remove oxidation liabilities with modest TMEM145 cost.',
  },
  {
    phase: 'Phase 4j',
    result: 'Avidity model',
    evidence: 'Original conservative Kd_mono 10 mM, C_eff 50 mM, N=4 gave Kd_eff 80 uM, failing by 1.6x; nominal case passed at 31 nM.',
    decision: 'The hypothesis needed a defensible Kd_mono floor, C_eff proof, or higher-valency geometry.',
  },
  {
    phase: 'Phase 4n',
    result: 'Monte Carlo stress test',
    evidence: 'P(PASS ideal)=0.803; P(PASS 30x entropic haircut)=0.590; Kd_mono was the dominant uncertainty.',
    decision: 'Prioritize WH2/F-actin affinity calibration over extra structural storytelling.',
  },
  {
    phase: 'Phase 4m',
    result: 'Kd soft-floor calibration',
    evidence: 'WH2 calibration fit R2=0.690; h09 predicted Kd_F=1.48 mM after 10,000x G-to-F scaling.',
    decision: 'Re-run conservative avidity gate: N=3, C_eff 50 mM gives Kd_eff 1.30 uM, passing the therapeutic threshold.',
  },
  {
    phase: 'Phase 4k',
    result: 'Fibril geometry caution',
    evidence: 'Tandem proxies reached four geometric WH2-actin contacts, but zero confidence-supported contacts across five AF3 jobs.',
    decision: 'Keep geometry as a caution, not as the S-tier promotion evidence.',
  },
  {
    phase: 'Phase 4p',
    result: 'Orthogonal WH2 controls',
    evidence: 'Boltz-2 placed h09 WH2 at ipTM 0.706; all negative controls were lower, while known positives remained higher or nearby.',
    decision: 'Supportive orthogonal check; not a replacement for wet WH2/F-actin binding data.',
  },
] as const;

export const constructRows = [
  {
    construct: 'tail91 v0',
    composition: 'WH2-RADA16-STRC aa 1620-1710',
    evidence: 'TMEM145 ipTM 0.68; actin ipTM 0.51 in the Phase 3 handoff.',
    role: 'Reference computational design.',
  },
  {
    construct: 'tail91 v2',
    composition: 'Cys11/52->Ser production candidate',
    evidence: 'TMEM145 median ipTM 0.55, sigma 0.010; actin median ipTM 0.51, sigma 0.011.',
    role: 'Developability-safe candidate for synthesis discussion.',
  },
  {
    construct: 'RADA16-WH2 scaffold',
    composition: 'Self-assembling SAP backbone with actin-binding handle',
    evidence: 'Top Phase 1 scaffold; later carried into PEPTIDE_TAIL91.',
    role: 'Assembly and multivalency engine.',
  },
] as const;

export const avidityRows = [
  {
    run: 'Phase 4j conservative',
    inputs: 'Kd_mono 10 mM, C_eff 50 mM, N=4',
    output: 'Kd_eff 80 uM',
    call: 'Fail by 1.6x; forced the next proof.',
  },
  {
    run: 'Phase 4j nominal',
    inputs: 'Kd_mono 5 mM, C_eff 100 mM, N=5',
    output: 'Kd_eff 31 nM',
    call: 'Pass, but too assumption-heavy to carry the claim alone.',
  },
  {
    run: 'Phase 4n Monte Carlo',
    inputs: 'Kd_mono 1-20 mM, C_eff 13.4-500 mM, N=3-6',
    output: '80.3% ideal pass; 59.0% stressed pass',
    call: 'Identified Kd_mono as the bottleneck.',
  },
  {
    run: 'Phase 4m conservative rerun',
    inputs: 'Kd_mono 1.48 mM, C_eff 50 mM, N=3',
    output: 'Kd_eff 1.30 uM',
    call: 'Pass; promotion from tentative-S to S.',
  },
  {
    run: 'Phase 4m higher valency',
    inputs: 'Kd_mono 1.48 mM, C_eff 50 mM, N=4',
    output: 'Kd_eff 38.39 nM',
    call: 'Pass with margin if four contacts are physically available.',
  },
] as const;

export const precedentRows = [
  {
    established: 'STRC-related hearing loss is a human DFNB16 phenotype with mostly mild/moderate onset and a long residual-cell window.',
    citation: 'GeneReviews STRC / DFNB16 (2024)',
    extrapolation: 'h09 scopes only STRC-loss ears with residual OHC function, not broad hearing loss.',
    risk: 'Late-stage disease or absent OHCs may be outside the response window.',
  },
  {
    established: 'Horizontal top connectors are mechanically load-bearing in mature OHC bundles. Strc-null/top-connector loss reduces bundle stiffness ~60% and damping ~74%.',
    citation: 'Verpy 2008; Verpy 2011; Rivera/Dulon 2019',
    extrapolation: 'A synthetic connector-like material could be useful if it restores even partial inter-stereocilia cohesion.',
    risk: 'Native STRC geometry and timed maturation may matter more than generic mechanical bridging.',
  },
  {
    established: 'WH2 domains bind G-actin with measured ITC Kd: WAVE2 52 nM, WIP 160 nM, MIM 230 nM, WASP 250 nM, Tbeta4 760 nM.',
    citation: 'Chereau 2005 PNAS, Fig. 1C',
    extrapolation: 'h09 uses a WH2-class actin-binding motif; Phase 4m places h09 in the WH2-family soft-floor band.',
    risk: 'WH2 x F-actin side-binding Kd is not directly measured; Phase 4m is a soft floor, not a measurement.',
  },
  {
    established: 'Polyvalency/avidity can convert weak monovalent binding into useful effective binding when geometry is matched (Mammen / Jencks / Karpen formalism).',
    citation: 'Mammen, Choi, Whitesides 1998 Angew Chem',
    extrapolation: "h09's RADA16 display is meant to turn weak mono WH2 contacts into cooperative effective binding.",
    risk: 'Flexible or misregistered display can erase avidity; Phase 4k remained confidence-inconclusive on geometry.',
  },
  {
    established: 'RADA16/SAP scaffolds self-assemble into nanofibers and hydrogels; characterized fiber length 615 +/- 104 nm; 0.6-3 mM hydrogel formation range.',
    citation: 'Yokoi/Kinoshita/Zhang 2005 PNAS',
    extrapolation: 'h09 uses a low-fraction tail-bearing RADA16 blend, not pure 134-aa tail91 self-assembly.',
    risk: 'Tail91 is far beyond characterized appendage length; rheology/self-assembly must be measured first.',
  },
  {
    established: 'Round-window entry into cochlear perilymph is experimentally quantifiable; primary RWM permeability and scala-tympani clearance parameters exist.',
    citation: 'Salt & Ma 2001 Hear Res',
    extrapolation: 'h09 delivery is a local exposure feasibility problem, not systemic delivery.',
    risk: '14 kDa peptide/hydrogel permeation is extrapolated; target-zone localization remains unproved.',
  },
  {
    established: 'Independent Boltz-2 WH2 controls do not place h09 below sequence-destroyed negatives.',
    citation: 'Phase 4p, this work (2026-04-27)',
    extrapolation: 'The Phase 4m WH2-family claim survives one matched alternative-model sanity check.',
    risk: 'Boltz does not recover the ITC rank order cleanly; supportive only, not affinity evidence.',
  },
  {
    established: "h09's novelty is the combination: local peptide-material display + WH2/avidity + synthetic HTC-style mechanical rescue.",
    citation: 'h09 hub, this work',
    extrapolation: 'The reviewer question is whether this combination merits the minimum feasibility package.',
    risk: 'The combination has no precedent; it must earn feasibility through localization, toxicity, and rescue assays.',
  },
] as const;

export const objectionRows = [
  {
    objection: 'AF3 is not a Kd predictor.',
    answer: 'Agree. Phase 4m uses matched WH2-family rank calibration only, not absolute AF3 thermodynamics.',
    risk: 'The rank relation could be a model artifact.',
    nextCheck: 'Run an orthogonal WH2 panel (positive references + h09 + scrambled + alanine) under one protocol. Phase 4p Boltz-2 already does the first pass.',
  },
  {
    objection: 'WH2 x G-actin does not prove WH2 x F-actin.',
    answer: 'Agree. Phase 4m uses Tbeta4 G-to-F scaling as an explicit soft floor.',
    risk: 'Scaling may not transfer to h09 WH2 or to filament-side binding.',
    nextCheck: 'Label Kd_F = 1.48 mM as a soft floor, not a measurement. Reserve a filament assay for paper-grade C_eff / F-actin support.',
  },
  {
    objection: 'Short WH2 motifs are weakly represented by structure models.',
    answer: 'Phase 4k is filed as INCONCLUSIVE and h09 S does not rely on it.',
    risk: 'A reviewer may discount all motif-model outputs.',
    nextCheck: 'Lead with primary WH2 Kd literature and Phase 4m calibration; use Phase 4k as method discipline, not support.',
  },
  {
    objection: 'Hydrogel may not reach the hair-bundle target zone.',
    answer: 'Delivery remains a feasibility premise, not a closed proof.',
    risk: 'Ototopical / RWM / LIFU delivery could fail before biology matters.',
    nextCheck: 'Run the precedent-separated delivery panel: RWM permeability, LIFU permeation, peptide localization, hair-bundle imaging.',
  },
  {
    objection: 'Self-assembly could cause local toxicity or occlusion.',
    answer: 'Current dose window is computational and class-based only.',
    risk: 'Local ototoxicity, vestibular exposure, or gel occlusion could kill the route.',
    nextCheck: 'Step 4 of the feasibility package: explant localization, dose titration, washout, toxicity before any rescue claim.',
  },
  {
    objection: 'A synthetic bridge may not rescue STRC-null OHC mechanics.',
    answer: 'h09 is a synthetic substitute concept, not STRC restoration.',
    risk: 'Native STRC geometry, timing, or protein interactions may be required.',
    nextCheck: 'Define success as partial mechanical rescue. Ask reviewers which assay can falsify mechanical benefit fastest.',
  },
  {
    objection: 'Patient window may be narrower than assumed.',
    answer: 'Best fit is STRC-loss with residual outer hair cells and mild/moderate phenotype.',
    risk: 'Late disease or advanced OHC loss may be unrecoverable.',
    nextCheck: 'Specify genotype, residual-cell inclusion criteria; avoid broad STRC-cure language.',
  },
  {
    objection: 'Phase 4k did not independently pass.',
    answer: 'Correct. Phase 4k is a caution: tandem-proxy geometric contacts yes, confidence-supported contacts no.',
    risk: 'Reviewer may treat this as contradiction rather than model limitation.',
    nextCheck: 'State the hierarchy plainly: Phase 4m is gate-closing; Phase 4k is a negative-control-style limitation.',
  },
  {
    objection: 'The hypothesis is a device/material concept, not a gene therapy.',
    answer: 'Correct. h09 is local peptide-material therapy.',
    risk: 'Wrong reviewer set will reject for the wrong reasons.',
    nextCheck: 'Send to biomaterials and inner-ear delivery reviewers before gene-therapy-only review.',
  },
  {
    objection: 'Addressable population may be overstated.',
    answer: 'Scope is STRC-loss with viable residual OHCs.',
    risk: 'Impact claims could look inflated.',
    nextCheck: 'Use conservative prevalence and response-window language; do not imply universal STRC-null rescue.',
  },
] as const;

export const boltzPanel = [
  { construct: 'MIM',                role: 'positive',         iptm: '0.925', read: 'high positive' },
  { construct: 'Tbeta4',             role: 'weak reference',   iptm: '0.914', read: 'false-high vs ITC' },
  { construct: 'WASP',               role: 'positive',         iptm: '0.873', read: 'positive' },
  { construct: 'h09 WH2',            role: 'query',            iptm: '0.706', read: 'in positive-reference band' },
  { construct: 'WAVE2',              role: 'positive',         iptm: '0.699', read: 'positive, not rank-best' },
  { construct: 'scrambled h09 WH2',  role: 'negative',         iptm: '0.633', read: 'below h09, close' },
  { construct: 'alanine h09 WH2',    role: 'negative',         iptm: '0.588', read: 'below h09' },
  { construct: 'WIP',                role: 'positive',         iptm: '0.340', read: 'false-low positive' },
  { construct: 'poly-GS spacer',     role: 'negative',         iptm: '0.319', read: 'low negative' },
] as const;

export const feasibilitySteps = [
  {
    n: '1',
    title: 'Peptide synthesis & QC',
    question: 'Can the h09 material be made cleanly enough to test?',
    materials: 'tail91_v2, plain RADA16, scrambled control, alanine-WH2 control, fluorescent tracer aliquot.',
    readouts: 'HPLC purity (>=90% discovery, >=95% explant), LC-MS intact mass, analytical HPLC trace, endotoxin for explant lots, lyophilized solubility.',
    pass: 'Identity / purity / solubility / endotoxin all clear; tail91_v2 + plain RADA16 in hand.',
    stop: 'Peptide cannot be made cleanly or remains insoluble/aggregated under assay buffer.',
  },
  {
    n: '2',
    title: 'Self-assembly',
    question: 'Do low-fraction tail91_v2 / RADA16 blends preserve beta-sheet fiber and gel behavior?',
    materials: '0%, 1%, 5%, 10%, 100% tail91_v2 in plain RADA16; matched buffer / pH / salt / temp.',
    readouts: 'CD spectroscopy (216-218 nm beta-sheet), AFM/TEM fiber morphology, DLS aggregate distribution, rheology G′ > G″, visual handling.',
    pass: '5% tail91_v2 blend forms RADA16-like beta-sheet/fiber/gel material without large-aggregate failure.',
    stop: 'Long-tail blend destroys assembly or forms micron-scale aggregates.',
  },
  {
    n: '3',
    title: 'Binding assays',
    question: 'Does h09 separate from sequence-negative controls on actin and TMEM145 axes?',
    materials: 'h09 tail91_v2 blend (Step 2 pass condition), plain RADA16, scrambled control, alanine-WH2, poly-GS, optional WH2 positive reference.',
    readouts: 'A) G-actin: fluorescence anisotropy / MST / BLI / pull-down. B) F-actin: co-sedimentation, TIRF bundle imaging. C) TMEM145: best available extracellular/GOLD proxy. Dose curves with all controls on the same axes.',
    pass: 'Dose-dependent h09 signal separates from scrambled / alanine / poly-GS controls on each axis.',
    stop: 'Negatives match h09, or no F-actin / TMEM145 signal is separable.',
  },
  {
    n: '4',
    title: 'Cochlear explant: localization & toxicity',
    question: 'Does the material reach the hair-bundle region without obvious local harm?',
    materials: 'Organ of Corti / cochlear explants. Exposures 0.3, 1, 3, 10 uM peptide-equivalent + vehicle + plain RADA16 controls.',
    readouts: 'Hair-bundle localization (fluorescent tracer aliquot + unlabeled toxicity controls), washout/reversibility, OHC survival (myo7a / prestin / phalloidin), bundle morphology, MET / FM dye / Ca proxy.',
    pass: 'Bundle-region enrichment at non-toxic doses, no irreversible occlusion, no OHC loss vs controls.',
    stop: 'No bundle-zone localization, OHC toxicity, irreversible occlusion, or vehicle-independent clumping.',
  },
  {
    n: '5',
    title: 'Mechanical rescue (gated)',
    question: 'Does h09 produce a measurable mechanical / functional benefit?',
    materials: 'Vehicle, plain RADA16, scrambled/alanine control, h09 at the best Step 4 non-toxic/localizing dose, optional positive rescue if available.',
    readouts: 'Hair-bundle stiffness / deflection, MET-current or FM-dye proxy, high-speed bundle motion coherence, Strc-null vs WT comparison if available. Predefine the primary endpoint before starting.',
    pass: 'h09 improves the prespecified primary endpoint vs vehicle and sequence-negative controls without toxicity.',
    stop: 'No improvement, or apparent improvement is paired with toxicity/occlusion.',
  },
] as const;

export const sections = [
  { id: 'ask',         label: 'Exact ask' },
  { id: 'claims',      label: 'Claims & non-claims' },
  { id: 'mechanism',   label: 'Mechanism' },
  { id: 'evidence',    label: 'Phase evidence' },
  { id: 'constructs',  label: 'Constructs' },
  { id: 'avidity',     label: 'Avidity gate' },
  { id: 'precedents',  label: 'Precedent table' },
  { id: 'objections',  label: 'Objection table' },
  { id: 'controls',    label: 'Boltz-2 controls' },
  { id: 'feasibility', label: 'Feasibility package' },
  { id: 'boundary',    label: 'Claim boundary' },
] as const;
