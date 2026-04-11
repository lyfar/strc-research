#!/usr/bin/env python3
"""Populate how_to_use, api_url, example_query for all tools."""
import sqlite3
from pathlib import Path

DB = Path(__file__).parent.parent / "data/tools.db"

USAGE = {
    # === VARIANT EFFECT PREDICTION ===
    "alphamissense": (
        "Web: AlphaFold DB heatmap. API: download pre-computed TSV from GCS. Python: pandas on TSV.",
        "https://alphafold.ebi.ac.uk/entry/Q7RTU9",
        'curl "https://alphafold.ebi.ac.uk/api/prediction/Q7RTU9" | jq .[0].amAnnotationUrl',
    ),
    "revel": (
        "Pre-computed scores. Download from sites.google.com or extract via dbNSFP/myvariant.info.",
        "https://myvariant.info/v1/variant/chr15:g.43600551T%3EG?fields=dbnsfp.REVEL",
        'curl "https://myvariant.info/v1/variant/chr15:g.43600551T%3EG?fields=dbnsfp.REVEL"',
    ),
    "cadd": (
        "Web: cadd.gs.washington.edu/snv. API: REST. Best via Ensembl VEP with CADD plugin.",
        "https://cadd.gs.washington.edu/api/v1.0/GRCh38/15:43600551",
        "curl 'https://rest.ensembl.org/vep/human/region/15:43600551-43600551:1/C?CADD=1' -H 'Content-Type:application/json'",
    ),
    "dann": (
        "Pre-computed scores via dbNSFP. No standalone API.",
        "https://myvariant.info/v1/variant/chr15:g.43600551T%3EG?fields=dbnsfp.DANN",
        'curl "https://myvariant.info/v1/variant/chr15:g.43600551T%3EG?fields=dbnsfp.DANN"',
    ),
    "clinpred": (
        "Download pre-computed scores from Google Sites or via dbNSFP.",
        "https://myvariant.info/v1/variant/chr15:g.43600551T%3EG?fields=dbnsfp.ClinPred",
        'curl "https://myvariant.info/v1/variant/chr15:g.43600551T%3EG?fields=dbnsfp.ClinPred"',
    ),
    "fathmm-mkl": (
        "Web form at fathmm.biocompute.org.uk. Also in dbNSFP.",
        "http://fathmm.biocompute.org.uk/",
        "Input: chr15 43600551 T G (web form). Or via dbNSFP.",
    ),
    "metarnn": (
        "Site down. Scores available via dbNSFP/myvariant.info.",
        "https://myvariant.info/v1/variant/chr15:g.43600551T%3EG?fields=dbnsfp.MetaRNN",
        'curl "https://myvariant.info/v1/variant/chr15:g.43600551T%3EG?fields=dbnsfp.MetaRNN"',
    ),
    "bayesdel": (
        "Pre-computed via dbNSFP. No standalone web interface.",
        "https://myvariant.info/v1/variant/chr15:g.43600551T%3EG?fields=dbnsfp.BayesDel",
        'curl "https://myvariant.info/v1/variant/chr15:g.43600551T%3EG?fields=dbnsfp.BayesDel"',
    ),
    "vest4": (
        "Site down. Scores in dbNSFP.",
        "https://myvariant.info/v1/variant/chr15:g.43600551T%3EG?fields=dbnsfp.VEST4",
        'curl "https://myvariant.info/v1/variant/chr15:g.43600551T%3EG?fields=dbnsfp.VEST4"',
    ),
    "varsome": (
        "Web: paste variant in search bar. API: REST (free account for rate limits).",
        "https://api.varsome.com/lookup/15-43600551-A-C?add-ACMG-annotation=1",
        "https://varsome.com/variant/hg38/NM_153700.2:c.4976A>C",
    ),
    "intervar": (
        "Web form: paste chr/pos/ref/alt. Automated ACMG classification.",
        "https://wintervar.wglab.org/",
        "Input: Chr15 43600551 43600551 A C (web form)",
    ),
    "franklin": (
        "Web: paste HGVS or VCF. Requires free account.",
        "https://franklin.genoox.com",
        "Search: NM_153700.2:c.4976A>C",
    ),
    "dbnsfp": (
        "Download full DB (35GB) or query via myvariant.info API. Contains 40+ precomputed scores.",
        "https://myvariant.info/v1/variant/chr15:g.43600551T%3EG?fields=dbnsfp",
        'curl "https://myvariant.info/v1/variant/chr15:g.43600551T%3EG?fields=dbnsfp" | python3 -m json.tool',
    ),
    "mpc": (
        "Pre-computed via dbNSFP. Missense constraint score.",
        "https://myvariant.info/v1/variant/chr15:g.43600551T%3EG?fields=dbnsfp.MPC",
        'curl "https://myvariant.info/v1/variant/chr15:g.43600551T%3EG?fields=dbnsfp.MPC"',
    ),
    "linsight": (
        "Pre-computed via dbNSFP. Non-coding focused but has coding scores.",
        "https://myvariant.info/v1/variant/chr15:g.43600551T%3EG?fields=dbnsfp.LINSIGHT",
        'curl "https://myvariant.info/v1/variant/chr15:g.43600551T%3EG?fields=dbnsfp.LINSIGHT"',
    ),
    "fitcons": (
        "Pre-computed via dbNSFP. Fitness consequences of functional annotation.",
        "https://myvariant.info/v1/variant/chr15:g.43600551T%3EG?fields=dbnsfp.fitCons",
        'curl "https://myvariant.info/v1/variant/chr15:g.43600551T%3EG?fields=dbnsfp.fitCons"',
    ),
    "eigen-eigen-pc": (
        "Pre-computed via dbNSFP. Unsupervised functional score.",
        "https://myvariant.info/v1/variant/chr15:g.43600551T%3EG?fields=dbnsfp.Eigen",
        'curl "https://myvariant.info/v1/variant/chr15:g.43600551T%3EG?fields=dbnsfp.Eigen"',
    ),

    # === STRUCTURAL BIOLOGY ===
    "alphafold-db": (
        "Web: search UniProt ID. API: REST for PDB/PAE/pLDDT.",
        "https://alphafold.ebi.ac.uk/api/prediction/Q7RTU9",
        'curl "https://alphafold.ebi.ac.uk/api/prediction/Q7RTU9" | jq .[0].pdbUrl',
    ),
    "alphafold-3": (
        "Web: alphafoldserver.com. Submit sequences for complex prediction. 10 jobs/day free.",
        "https://alphafoldserver.com",
        "Submit STRC sequence (Q7RTU9) + interaction partner. Download PDB result.",
    ),
    "dynamut": (
        "Web: upload PDB + specify mutation. Returns ddG prediction.",
        "https://biosig.lab.uq.edu.au/dynamut/",
        "Input: PDB AF-Q7RTU9-F1, Chain A, Mutation E1659A. Result: -0.913 kcal/mol",
    ),
    "mcsm": (
        "Web: upload PDB + specify mutation. Graph-based stability prediction.",
        "https://biosig.lab.uq.edu.au/mcsm/",
        "Input: PDB AF-Q7RTU9-F1, Chain A, Mutation E1659A",
    ),
    "duet": (
        "Web: upload PDB + mutation. Combines mCSM + SDM predictions.",
        "https://biosig.lab.uq.edu.au/duet/",
        "Input: PDB AF-Q7RTU9-F1, Chain A, Mutation E1659A",
    ),
    "swiss-model": (
        "Web: paste sequence, auto-find template. API: REST for automated modeling.",
        "https://swissmodel.expasy.org/automodel",
        "Input: STRC protein sequence (Q7RTU9). Template: AlphaFold structure.",
    ),
    "esmfold": (
        "API: POST protein sequence, returns PDB. No GPU needed.",
        "https://esmatlas.com/resources?action=fold",
        'curl -X POST --data "STRC_sequence" https://api.esmatlas.com/foldSequence/v1/pdb/',
    ),
    "iupred3": (
        "Web: paste sequence or UniProt ID. Returns disorder prediction per residue.",
        "https://iupred3.elte.hu",
        "Input: Q7RTU9. Check position 1659 disorder score.",
    ),
    "signalp": (
        "Web: paste protein sequence. Predicts signal peptide and cleavage site.",
        "https://services.healthtech.dtu.dk/services/SignalP-6.0/",
        "Input: STRC full sequence. Result: signal peptide at N-terminus.",
    ),
    "netgpi": (
        "Web: paste protein sequence. Predicts GPI-anchor signal.",
        "https://services.healthtech.dtu.dk/services/NetGPI-1.1/",
        "Input: STRC full sequence. Result: GPI anchor at C-terminus.",
    ),
    "netnglyc": (
        "Web: paste protein sequence. Predicts N-glycosylation sites.",
        "https://services.healthtech.dtu.dk/services/NetNGlyc-1.0/",
        "Input: STRC full sequence. Check glycosylation near position 1659.",
    ),

    # === SPLICING ===
    "spliceai": (
        "Web: spliceailookup.broadinstitute.org. Also via Ensembl VEP plugin.",
        "https://spliceailookup.broadinstitute.org/",
        "curl 'https://rest.ensembl.org/vep/human/region/15:43600551-43600551:1/C?SpliceAI=1' -H 'Content-Type:application/json'",
    ),
    "pangolin": (
        "CLI: pip install pangolin-splice. Run on VCF input.",
        "https://github.com/bw2/Pangolin",
        "pangolin input.vcf output.vcf --genome hg38",
    ),

    # === POPULATION DATABASES ===
    "gnomad": (
        "Web: search variant. API: GraphQL.",
        "https://gnomad.broadinstitute.org/api",
        'POST {"query":"{ variant(variantId: \\"15-43600551-A-C\\", dataset: gnomad_r4) { variant_id exome { ac af } } }"}',
    ),
    "clinvar": (
        "Web: search gene or variant. API: NCBI E-utilities.",
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar",
        'curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=STRC[gene]&retmode=json"',
    ),
    "dbsnp": (
        "Web: search rsID or coordinates. API: NCBI SPDI.",
        "https://api.ncbi.nlm.nih.gov/variation/v0/",
        'curl "https://api.ncbi.nlm.nih.gov/variation/v0/spdi/NC_000015.10:43600550:T:G/rsids"',
    ),
    "lovd": (
        "Web: search gene. API: REST.",
        "https://databases.lovd.nl/shared/api/rest.php/variants/STRC",
        'curl "https://databases.lovd.nl/shared/api/rest.php/variants/STRC?format=application/json"',
    ),
    "bravo-topmed": (
        "Web: search variant by coordinates.",
        "https://bravo.sph.umich.edu/freeze8/hg38/",
        "https://bravo.sph.umich.edu/freeze8/hg38/variant/snv/15-43600551-A-C",
    ),
    "clingen": (
        "Web: search gene/variant. Allele Registry API for ClinGen IDs.",
        "https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_caid",
        "https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_caid?caid=CA392159910",
    ),

    # === NOMENCLATURE ===
    "variantvalidator": (
        "API: REST. Validates and normalizes HGVS descriptions.",
        "https://rest.variantvalidator.org/VariantValidator/variantvalidator/",
        'curl "https://rest.variantvalidator.org/VariantValidator/variantvalidator/GRCh38/NM_153700.2:c.4976A>C/all?content-type=application/json"',
    ),
    "mutalyzer": (
        "API: REST. Normalizes HGVS notation, checks syntax.",
        "https://mutalyzer.nl/api/normalize/",
        'curl "https://mutalyzer.nl/api/normalize/NM_153700.2:c.4976A%3EC"',
    ),

    # === GENE-LEVEL ===
    "uniprot": (
        "Web: search protein. API: REST for sequence, features, domains.",
        "https://rest.uniprot.org/uniprotkb/Q7RTU9",
        'curl "https://rest.uniprot.org/uniprotkb/Q7RTU9.json"',
    ),
    "ensembl": (
        "API: REST. VEP, variant consequences, gene info.",
        "https://rest.ensembl.org",
        "curl 'https://rest.ensembl.org/vep/human/region/15:43600551-43600551:1/C?content-type=application/json'",
    ),
    "ucsc": (
        "Web: genome browser visualization. API: REST for track data.",
        "https://api.genome.ucsc.edu",
        "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr15:43600000-43601000",
    ),

    # === LITERATURE ===
    "pubmed": (
        "Web: search. API: NCBI E-utilities.",
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
        'curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=STRC+hearing+loss&retmax=5&retmode=json"',
    ),
    "litvar-20": (
        "API: search variant mentions in literature.",
        "https://www.ncbi.nlm.nih.gov/research/litvar2-api/",
        'curl "https://www.ncbi.nlm.nih.gov/research/litvar2-api/variant/search?query=STRC%20E1659A"',
    ),
    "semantic-scholar": (
        "API: search papers by keyword.",
        "https://api.semanticscholar.org/graph/v1/paper/search",
        'curl "https://api.semanticscholar.org/graph/v1/paper/search?query=STRC+stereocilin&limit=5&fields=title,year"',
    ),
    "pubtator": (
        "API: search for variant/gene mentions with NER annotations.",
        "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/",
        'curl "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/search/?text=STRC%20E1659A"',
    ),

    # === REGULATORY ===
    "regulomedb": (
        "Web: search by coordinates. REST API available.",
        "https://regulomedb.org/regulome-search/",
        "https://regulomedb.org/regulome-search/?regions=chr15:43600551-43600551&genome=GRCh38",
    ),
    "encode": (
        "Web: search experiments. API: REST for metadata and files.",
        "https://www.encodeproject.org/search/",
        'curl "https://www.encodeproject.org/search/?type=Experiment&target.genes.symbol=STRC&format=json"',
    ),
    "gtex-portal": (
        "Web: search gene expression by tissue. API: REST.",
        "https://gtexportal.org/api/v2/",
        'curl "https://gtexportal.org/api/v2/expression/geneExpression?gencodeId=ENSG00000242866"',
    ),

    # === CLINICAL ===
    "omim": (
        "Web: search gene/disease. API requires key.",
        "https://omim.org",
        "https://omim.org/entry/606440 (STRC) / https://omim.org/entry/603720 (DFNB16)",
    ),
    "open-targets": (
        "API: GraphQL for gene-disease associations.",
        "https://api.platform.opentargets.org/api/v4/graphql",
        'POST {"query":"{ target(ensemblId: \\"ENSG00000242866\\") { approvedSymbol associatedDiseases { count } } }"}',
    ),
    "decipher": (
        "Web: search gene for phenotype-genotype correlations.",
        "https://www.deciphergenomics.org",
        "https://www.deciphergenomics.org/gene/STRC/overview",
    ),
    "human-protein-atlas": (
        "Web: search gene for tissue expression. No JSON API.",
        "https://www.proteinatlas.org",
        "https://www.proteinatlas.org/ENSG00000242866-STRC",
    ),
    "clinicaltrials": (
        "Web: search trials. API: REST.",
        "https://clinicaltrials.gov/api/v2/studies",
        'curl "https://clinicaltrials.gov/api/v2/studies?query.cond=STRC+hearing+loss&pageSize=5"',
    ),

    # === HEARING LOSS ===
    "dvd": (
        "Web: search gene for hearing loss variant classifications.",
        "https://deafnessvariationdatabase.org",
        "https://deafnessvariationdatabase.org/gene/STRC",
    ),
    "hhl": (
        "Web: browse gene list for hearing loss loci.",
        "https://hereditaryhearingloss.org",
        "https://hereditaryhearingloss.org (search STRC / DFNB16)",
    ),
    "audiogene": (
        "Web: input audiogram shape to predict causative gene.",
        "https://audiogene.eng.uiowa.edu/",
        "Input: bilateral moderate SNHL audiogram. Expected: STRC as candidate.",
    ),
    "gear": (
        "Web: search gene for inner ear single-cell expression.",
        "https://umgear.org/",
        "Search: STRC. View expression in cochlear cell types.",
    ),

    # === CONSERVATION ===
    "consurf": (
        "Web: submit UniProt ID or PDB. Returns per-residue conservation scores.",
        "https://consurf.tau.ac.il/",
        "Input: UniProt Q7RTU9. Check conservation grade at position 1659.",
    ),
    "orthodb": (
        "Web: search gene for orthologs across species.",
        "https://www.orthodb.org/v11/",
        "https://www.orthodb.org/v11/search?query=STRC&species=9606",
    ),
    "clustal-omega": (
        "Web: paste multiple sequences. API: EMBL-EBI Job Dispatcher.",
        "https://www.ebi.ac.uk/jdispatcher/msa/clustalo",
        "Input: STRC orthologs from Human, Mouse, Rat, Cow, Dog. Check position 1659.",
    ),

    # === AI / FOUNDATION MODELS ===
    "alphagenome": (
        "Python API: predict variant effects on expression, splicing, chromatin.",
        "https://www.alphagenomedocs.com",
        'alphagenome.predict(chrom="chr15", pos=43600551, ref="T", alt="G", genome="hg38")',
    ),
}


def main():
    db = sqlite3.connect(str(DB))
    for slug, (how, api, example) in USAGE.items():
        db.execute(
            "UPDATE tools SET how_to_use = ?, api_url = ?, example_query = ? WHERE slug = ?",
            (how, api, example, slug),
        )
    db.commit()

    filled = db.execute("SELECT COUNT(*) FROM tools WHERE how_to_use != ''").fetchone()[0]
    total = db.execute("SELECT COUNT(*) FROM tools").fetchone()[0]
    print(f"Usage info populated: {filled}/{total} tools")
    db.close()


if __name__ == "__main__":
    main()
