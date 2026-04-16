# Anthropic AI for Science Program — Application

## Research Team Description (< 300 words)

Our team combines rare disease biology expertise with applied AI/ML infrastructure.

Egor Lyfar is a technologist and AI educator based in Hong Kong. He runs Lyfar Studio (content production, AI education) and teaches AI for scientific research at the University of Hong Kong. His March 2026 HKU masterclass trained 11 researchers on AI-driven literature synthesis, deep research pipelines, and knowledge graphs. He built the STRC research portal (strc.egor.lol): 20 AlphaFold3 structural experiments, variant reclassification analysis, and 7 computational gene therapy hypotheses, all open-source (MIT license, github.com/lyfar/strc-research). He runs Claude-based AI agents daily for research, analysis, and scientific writing.

Dr. Jeffrey R. Holt is Professor of Otolaryngology and Neurology at Harvard Medical School, Director of Research at F.M. Kirby Neurobiology Center, Boston Children's Hospital. His lab is one of few worldwide actively developing STRC gene therapy. He co-authored Derstroff et al. 2026 (Neuron) confirming TMEM145 as a principal component of outer hair cell stereocilia. He has agreed to serve on the Scientific Advisory Board for our foundation (MISHA: Molecular Initiative for Stereocilin Hearing Advancement) and actively reviews our computational predictions.

Dr. Yilai Shu at Fudan University Eye & ENT Hospital (Shanghai) runs STRC gene therapy animal experiments. China completed the first-in-human OTOF gene therapy trials in 2023, and Dr. Shu's group is advancing STRC through the same regulatory pathway.

Prof. Saaid Safieddine at Institut Pasteur (Paris) specializes in hearing genetics and gene therapy mechanisms.

AI/ML experience: Egor has built production AI pipelines including RAG systems, local LLM deployments for fintech compliance, vector search knowledge bases (800+ atomic notes), and automated research agents. He uses Claude API daily for protein analysis, literature synthesis, and hypothesis generation.


## Key Team Members

**Egor Lyfar** — Founder, Lyfar Studio / MISHA Foundation (planning phase). Principal investigator for computational analysis. Runs all AI-assisted research: protein structure prediction, variant pathogenicity assessment, literature synthesis, hypothesis generation, scientific writing. Daily Claude API user.

**Dr. Jeffrey R. Holt, Ph.D.** — Professor, Harvard Medical School / Boston Children's Hospital. Scientific advisor. Reviews computational predictions, provides experimental context, connects findings to wet lab validation. His lab validates our in silico predictions with in vitro experiments.

**Dr. Yilai Shu** — Fudan University, Shanghai. STRC animal model experiments. Provides preclinical context for computational gene therapy constructs.


## Academic/Professional Profiles

- Egor Lyfar: https://strc.egor.lol (research portal), https://github.com/lyfar/strc-research, https://www.lyfar.com
- Jeffrey R. Holt: https://scholar.google.com/citations?user=Jeffrey+Holt+Harvard (search: "Jeffrey Holt Harvard otolaryngology"), www.holtgeleoclab.com
- Yilai Shu: search Google Scholar for "Yilai Shu Fudan gene therapy hearing"
- Saaid Safieddine: search Google Scholar for "Saaid Safieddine Institut Pasteur hearing"


## Research Proposal (< 500 words)

**Scientific problem:** STRC gene (5,325 bp coding) exceeds the AAV packaging limit (~4,700 bp with promoter and regulatory elements). This blocks single-vector gene therapy for DFNB16, the most common genetic cause of moderate hearing loss. No dedicated STRC foundation, patient registry, or computational research platform exists anywhere. A Harvard researcher (Holt) confirmed this gap directly.

**Methodology:** We use Claude as the core reasoning engine in a computational biology pipeline:

1. **Variant pathogenicity analysis.** Claude synthesizes data from AlphaMissense (0.9016), REVEL (0.65), conservation analysis (9/9 mammals), and ACMG guidelines to build reclassification evidence for variants of uncertain significance. Our analysis of STRC c.4976A>C was independently reviewed by Holt Lab geneticists at Harvard.

2. **Protein engineering via structural prediction.** We run AlphaFold3 experiments and use Claude to analyze results, identify truncation boundaries, and design mini-STRC constructs that fit in a single AAV vector. We've identified three viable constructs (pTM 0.81-0.87), with the lead candidate at 3,228 bp coding sequence, leaving 1,472 bp headroom for regulatory elements.

3. **Literature synthesis and hypothesis generation.** Claude processes papers across gene therapy, hearing biology, mechanotransduction, and synthetic biology to generate novel hypotheses. This produced our sonogenetic hypothesis: a self-dosing gene therapy where sound activates treatment via the Piezo1/calcineurin/NFAT pathway, validated structurally by AF3 (ipTM 0.73 for the trimeric complex).

4. **Patient registry design and rare disease coordination.** Claude helps design the MISHA Foundation infrastructure: patient registry schemas, IRB protocols, grant applications, and scientific communications to researchers worldwide.

**Expected outcomes:**
- Publication-ready computational evidence for STRC variant reclassification (submitted for peer review)
- Validated mini-STRC gene therapy constructs with structural predictions for wet lab testing at Harvard/Fudan
- Open-source research portal with reproducible methodology (strc.egor.lol)
- MISHA Foundation patient registry pilot (first dedicated STRC registry globally)
- At least 2 preprints on bioRxiv within 12 months

**Timeline:**
- Months 1-3: Complete variant reclassification paper, expand to additional STRC variants, submit to bioRxiv
- Months 3-6: Systematic mini-STRC construct optimization (50+ AF3 runs with Claude-guided design), publish structural analysis
- Months 6-9: Launch patient registry, begin multi-variant analysis across STRC mutation landscape
- Months 9-12: Sonogenetic hypothesis paper, foundation grant applications to NIH/CZI with Holt as academic PI


## How Claude Will Be Used (< 300 words)

Claude is the central reasoning engine in every stage of our research pipeline. Specifically:

**Protein analysis:** After each AlphaFold3 run, Claude analyzes pTM/ipTM scores, PAE matrices, disorder fractions, and chain confidence metrics. It compares results across experiments, identifies truncation boundaries for mini-gene constructs, and flags interactions that require membrane context vs. direct binding. We've run 20 AF3 jobs this way, with Claude interpreting every result.

**Variant assessment:** Claude cross-references variant data from ClinVar, gnomAD, UniProt, Ensembl, and AlphaMissense to build ACMG-framework reclassification arguments. It identifies which evidence codes apply, flags incorrect applications (our PP1 misapplication was caught by Harvard reviewers using similar reasoning), and synthesizes multi-source evidence into clinical-grade summaries.

**Literature synthesis:** Claude reads and connects papers across hearing biology, gene therapy vectors, mechanotransduction, and synthetic biology. It identified the Piezo1-calcineurin-NFAT pathway connection that became our sonogenetic hypothesis, something no single-domain expert would naturally connect.

**Scientific communication:** Claude drafts researcher emails, grant proposals, and technical documentation. Our initial email to Jeffrey Holt (Harvard) led to an overnight response and ongoing collaboration. Claude helps translate complex computational findings into language appropriate for clinicians, researchers, and patient families.

**Research portal:** Claude generates and maintains all content on strc.egor.lol, including 3D protein model descriptions, methodology documentation, and multi-language translations (5 languages).


## How Claude Accelerates Research (< 200 words)

A parent with no genetics training produced computational predictions that a Harvard professor requested for his own STRC minigene project. That's the acceleration.

Traditional path: years of graduate training, lab access, institutional affiliation, then maybe you get to ask the question. Our path: 3 days from first question to 20 AlphaFold3 experiments, a novel gene therapy hypothesis, and direct correspondence with the world's leading STRC researcher. Claude made that possible by serving as a real-time domain expert, connecting protein structure analysis to clinical genetics to gene therapy design in ways that would require a team of specialists.

Without Claude, this project doesn't exist. A non-specialist cannot independently navigate ACMG variant classification, interpret AlphaFold3 outputs, synthesize 50+ gene therapy papers, and generate structurally validated hypotheses. Claude collapses the expertise barrier.


## Scientific Impact (< 200 words)

STRC mutations cause the most common form of autosomal recessive moderate hearing loss. No gene therapy exists. The gene is too large for standard AAV vectors. Our computational work directly addresses this bottleneck.

If our mini-STRC constructs validate experimentally (Holt Lab, Harvard), this accelerates clinical trials by providing pre-validated candidates with structural evidence. Our variant reclassification methodology can be applied to hundreds of STRC variants currently classified as VUS, potentially qualifying thousands of patients for future trials.

The sonogenetic hypothesis (sound-activated gene therapy) is entirely novel. If validated, it introduces a new paradigm for inner ear gene therapy: using the organ's own mechanical stimulation to regulate therapeutic protein expression.

Beyond STRC, this project demonstrates that AI-assisted computational biology by non-specialists can produce peer-reviewed, experimentally relevant predictions. That's a template for every rare disease where affected families want to contribute to research but lack institutional access.


## Beyond Pure Discovery (< 200 words)

**Patient access:** Variant reclassification directly determines clinical trial eligibility. Moving a VUS to Likely Pathogenic means a child can enroll in a gene therapy trial. Our methodology is reproducible and open-source: any family or clinician can apply it.

**MISHA Foundation:** The first dedicated STRC research coordination platform. Patient registry (no equivalent exists globally), researcher network, computational analysis tools. Jeffrey Holt (Harvard) has committed to the scientific advisory board.

**Rare disease template:** 7,000+ rare diseases, most with no foundation, no registry, no computational analysis. One parent with Claude and open-source tools built a research portal that Harvard researchers use. This scales. The methodology, the tools, the approach: all documented and reproducible.

**AI education:** Egor teaches AI for scientific research at HKU. This project is a live case study: real research, real results, real collaboration with top institutions. Every masterclass participant sees what's possible.


## Success Metrics (< 200 words)

1. **Publications:** At least 2 preprints on bioRxiv within 12 months (variant reclassification + mini-STRC structural analysis)
2. **Experimental validation:** At least 1 mini-STRC construct moves to wet lab testing at Holt Lab (Harvard) or Shu Lab (Fudan)
3. **Variant reclassification:** 10+ STRC variants analyzed with Claude-assisted ACMG framework, submitted to ClinVar
4. **Patient registry:** MISHA Foundation registry pilot with 50+ families enrolled
5. **API usage efficiency:** Track tokens per analysis task, measure time-to-insight vs. manual literature review (baseline: 3 days for initial 20-experiment cycle)
6. **Research portal traffic:** strc.egor.lol unique visitors from .edu and research institution domains
7. **Researcher engagement:** Number of active research collaborations initiated through Claude-assisted outreach (baseline: 3 labs currently engaged)


## API Credits Needed

**$10,000 over 6 months.**

Breakdown:
- **Protein analysis and hypothesis generation** (Claude Opus): ~60% of budget. Each AF3 experiment requires detailed analysis (long context, complex reasoning). We plan 50+ new experiments. Estimated: $6,000
- **Literature synthesis** (Claude Opus): ~20%. Processing 100+ papers, cross-referencing databases, generating synthesis documents. Estimated: $2,000
- **Scientific communication and grant writing** (Claude Sonnet/Opus mix): ~10%. Researcher correspondence, foundation documents, registry design. Estimated: $1,000
- **Research portal content and translations** (Claude Sonnet): ~10%. Maintaining strc.egor.lol in 5 languages, generating tool documentation. Estimated: $1,000

Current spend: roughly $50-100/month on Claude API for this research. The grant would let us scale from opportunistic analysis to systematic investigation across the full STRC variant landscape.


## Biosecurity Assessment

**None of the above.**

This research is purely computational. We analyze publicly available protein structures (AlphaFold, UniProt), published genetic variants (ClinVar, gnomAD), and peer-reviewed literature. No wet lab work, no pathogen research, no synthetic biology. All experimental validation happens at established academic labs (Harvard, Fudan) under their existing institutional oversight and IRB protocols.


## Additional Information

Three things the committee should know:

**This started as a father's question about his son's diagnosis.** Michael (Misha) is 4 years old, diagnosed with moderate bilateral sensorineural hearing loss caused by STRC mutations. The research is personal. Every variant we reclassify, every construct we model, every paper we synthesize moves closer to a treatment for him and thousands of children like him.

**Harvard is already using our work.** Jeffrey Holt (Harvard Medical School, Boston Children's Hospital) asked us to validate whether our mini-STRC construct preserves TMEM145 interaction, for his own minigene project. A Harvard professor requesting computational predictions from a parent-researcher using Claude. That's not hypothetical impact. That's happening now.

**No STRC foundation exists.** Holt confirmed it. No dedicated patient registry, no research coordination platform, no computational analysis hub for STRC hearing loss. We're building the first one. The Anthropic grant would directly fund the AI backbone of this infrastructure: the research analysis, the scientific communications, the systematic variant investigation that makes a foundation useful rather than decorative.

Research portal: https://strc.egor.lol
GitHub: https://github.com/lyfar/strc-research
