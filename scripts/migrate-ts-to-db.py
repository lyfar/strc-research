#!/usr/bin/env python3
"""Migrate tools.ts data into SQLite database.

Stores: slug, brain_note, name, category, short description, url, status,
tags, is_free, has_api, acmg_relevance, e1659a results.
Short descriptions extracted from Brain vault notes.
"""
import re, sqlite3, json
from pathlib import Path

ROOT = Path(__file__).parent.parent
TS_FILE = ROOT / "src/data/tools.ts"
DB_FILE = ROOT / "data/tools.db"
SCHEMA = ROOT / "scripts/init-db.sql"
BRAIN = Path.home() / "Library/Mobile Documents/iCloud~md~obsidian/Documents/Brain/notes"


def get_brain_data(note_name):
    """Extract short description and tags from Brain note."""
    path = BRAIN / f"{note_name}.md"
    if not path.exists():
        parts = note_name.split(' / ')
        if len(parts) == 2:
            path = BRAIN / f"{parts[0]} " / f"{parts[1]}.md"
    if not path.exists():
        return '', []

    content = path.read_text()

    # Tags from frontmatter
    tags = []
    fm = re.match(r'^---\s*\n(.*?)\n---', content, re.DOTALL)
    if fm:
        tm = re.search(r'tags:\s*\[([^\]]*)\]', fm.group(1))
        if tm:
            tags = [t.strip() for t in tm.group(1).split(',') if t.strip()]
    tags = [t for t in tags if t not in ('tool', 'strc-research')]

    # Short description: first line of "What It Does" or first paragraph
    body = re.sub(r'^---.*?---\s*', '', content, flags=re.DOTALL)
    body = re.sub(r'^# .+\n+', '', body)
    desc = ''
    wm = re.search(r'## What It Does\s*\n+(.+?)(?:\n\n|\n##)', body, re.DOTALL)
    if wm:
        lines = [l.strip().lstrip('- ') for l in wm.group(1).strip().split('\n') if l.strip()]
        if lines:
            desc = lines[0]
    if not desc:
        for line in body.split('\n'):
            line = line.strip()
            if line and not line.startswith(('#', '-', '|', '*', '`')):
                desc = line
                break
    desc = re.sub(r'\[([^\]]+)\]\([^\)]+\)', r'\1', desc)
    desc = re.sub(r'\*\*([^*]+)\*\*', r'\1', desc)
    return desc[:200], tags

# slug → Brain note filename (without .md)
BRAIN_MAP = {
    "alphafold-db": "AlphaFold Database",
    "alphafold-3": "AlphaFold 3 Server",
    "alphagenome": "AlphaGenome",
    "swiss-model": "SWISS-MODEL",
    "molstar": "Mol Star Viewer",
    "pymol": "PyMOL",
    "esmfold": "ESMFold",
    "iupred3": "IUPred3",
    "alphamissense": "AlphaMissense",
    "revel": "REVEL",
    "clinvar": "ClinVar",
    "spliceai": "SpliceAI",
    "cadd": "CADD",
    "varsome": "VarSome",
    "intervar": "InterVar",
    "franklin": "Franklin",
    "dbnsfp": "dbNSFP",
    "gnomad": "gnomAD",
    "seqr": "seqr",
    "decipher": "DECIPHER",
    "clingen": "ClinGen",
    "uniprot": "UniProt",
    "ucsc": "UCSC Genome Browser",
    "ensembl": "Ensembl REST API",
    "orthodb": "OrthoDB",
    "clustal-omega": "Clustal Omega",
    "omim": "OMIM",
    "dvd": "Deafness Variation Database",
    "hhl": "Hereditary Hearing Loss Homepage",
    "clinicaltrials": "ClinicalTrials.gov",
    "pubmed": "PubMed",
    "aav-capsid": "AAV Capsid Database",
    "asgct": "ASGCT Gene Therapy Database",
    "signalp": "SignalP",
    "netgpi": "NetGPI",
    "netnglyc": "NetNGlyc",
    "esm1v": "ESM1v",
    "eve-evolutionary-model-of-variant-effect": "EVE (Evolutionary model of Variant Effect)",
    "primateai-3d": "PrimateAI-3D",
    "mutscore": "MutScore",
    "gmvp": "gMVP",
    "metarnn": "MetaRNN",
    "bayesdel": "BayesDel",
    "mpc": "MPC",
    "dann": "DANN",
    "fathmm-mkl": "FATHMM-MKL",
    "linsight": "LINSIGHT",
    "eigen-eigen-pc": "Eigen",
    "fitcons": "fitCons",
    "clinpred": "ClinPred",
    "vest4": "VEST4",
    "snpeff": "SnpEff",
    "sparkinferno": "SparkINFERNO",
    "favor": "FAVOR",
    "rosettafold": "RoseTTAFold",
    "ifum": "IFUM",
    "foldx": "FoldX",
    "rosettaddg": "RosettaDDG",
    "ddgun-ddgun3d": "DDGun",
    "dynamut": "DynaMut",
    "mcsm": "mCSM",
    "sdm": "SDM",
    "duet": "DUET",
    "popmusic": "PoPMuSiC",
    "strum": "STRUM",
    "elaspic": "ELASPIC",
    "pangolin": "Pangolin",
    "absplice2": "AbSplice2",
    "maxentscan": "MaxEntScan",
    "genesplicer": "GeneSplicer",
    "nnsplice": "NNSplice",
    "spanr": "SPANR",
    "hal": "HAL",
    "spip": "SPiP",
    "splicerover": "SpliceRover",
    "mmsplice": "MMSplice",
    "fraser2": "FRASER2",
    "gtex-portal": "GTEx Portal",
    "encode": "ENCODE",
    "roadmap-epigenomics": "Roadmap Epigenomics",
    "regulomedb": "RegulomeDB",
    "deepsea": "DeepSEA",
    "basenji": "Basenji",
    "enformer": "Enformer",
    "tland": "TLand",
    "haploreg": "HaploReg",
    "piq": "PIQ",
    "bravo-topmed": "BRAVO / TOPMed",
    "uk-biobank": "UK Biobank",
    "dbsnp": "dbSNP",
    "exac-legacy": "ExAC",
    "alfa": "ALFA",
    "lovd": "LOVD",
    "biobank-japan-bbj-pheweb": "Biobank Japan",
    "cmdb-chinese-millionome": "CMDB",
    "kob-kdna": "KoB KDNA",
    "sage": "SAGE",
    "hgmd-professional": "HGMD Professional",
    "mastermind": "Mastermind",
    "civic": "CIViC",
    "pharmgkb": "PharmGKB",
    "dgidb": "DGIdb",
    "gwas-catalog": "GWAS Catalog",
    "open-targets": "Open Targets",
    "cgar": "CGAR",
    "human-protein-atlas": "Human Protein Atlas",
    "expression-atlas": "Expression Atlas",
    "allen-brain-atlas": "Allen Brain Atlas",
    "archs4": "ARCHS4",
    "audiogene": "AudioGene",
    "hiedra": "HIEDRA",
    "gear": "gEAR",
    "sensorion-aav-platform": "Sensorion AAV Platform",
    "db-oto-trial-resources": "DB-OTO",
    "consurf": "ConSurf",
    "evcouplings": "EVcouplings",
    "rate4site": "Rate4Site",
    "classifycnv": "ClassifyCNV",
    "annotsv": "AnnotSV",
    "cnvnator": "CNVnator",
    "gatk-gcnv": "GATK gCNV",
    "manta": "Manta",
    "delly": "Delly",
    "cnest": "CNest",
    "variantvalidator": "VariantValidator",
    "mutalyzer": "Mutalyzer",
    "crossmap": "CrossMap",
    "litvar-20": "LitVar 2.0",
    "pubtator": "PubTator",
    "tmvar": "tmVar",
    "europe-pmc": "Europe PMC",
    "semantic-scholar": "Semantic Scholar",
    "litsuggest": "LitSuggest",
    "litsense": "LitSense",
    "galaxy-project": "Galaxy Project",
    "terrabio": "Terra.bio",
    "illumina-connected-insights": "Illumina Connected Insights",
    "variant-tools-vtools": "variant tools (vtools)",
    "evo-evo-2": "Evo / Evo 2",
    "nucleotide-transformer-ntv3": "Nucleotide Transformer (NTv3)",
    "caduceus": "Caduceus",
    "hyenadna": "HyenaDNA",
    "dnabert-2": "DNABERT-2",
}


def parse_tools_ts():
    """Extract tool entries from tools.ts."""
    content = TS_FILE.read_text()
    tools = []

    blocks = re.split(r'\n  \},\n', content)

    for block in blocks:
        sm = re.search(r"slug: '([^']+)'", block)
        if not sm:
            continue
        slug = sm.group(1)

        def extract(field):
            m = re.search(rf"{field}: '((?:[^'\\]|\\.)*)'", block)
            return m.group(1).replace("\\'", "'") if m else ''

        name = extract('name')
        category = extract('category')
        url = extract('url')
        status = extract('status')

        e1659a = None
        em = re.search(r"e1659a: \{ score: '([^']*)', interpretation: '([^']*)' \}", block)
        if em:
            e1659a = {'score': em.group(1), 'interpretation': em.group(2)}

        acmg = extract('acmgRelevance') or None
        is_free = 'isFree: true' in block
        has_api = 'hasApi: true' in block

        tools.append({
            'slug': slug, 'name': name, 'category': category,
            'url': url, 'status': status, 'e1659a': e1659a,
            'acmg': acmg, 'is_free': is_free, 'has_api': has_api,
        })

    return tools


def main():
    DB_FILE.parent.mkdir(exist_ok=True)
    if DB_FILE.exists():
        DB_FILE.unlink()

    db = sqlite3.connect(str(DB_FILE))
    db.executescript(SCHEMA.read_text())

    tools = parse_tools_ts()
    print(f"Parsed {len(tools)} tools from tools.ts")

    brain_found = 0
    for t in tools:
        brain_note = BRAIN_MAP.get(t['slug'], t['name'])
        desc, tags = get_brain_data(brain_note)
        if desc:
            brain_found += 1

        db.execute("""
            INSERT INTO tools (slug, brain_note, name, category, description, url, status, is_free, has_api, acmg_relevance)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (t['slug'], brain_note, t['name'], t['category'], desc,
              t['url'], t['status'], int(t['is_free']), int(t['has_api']), t['acmg']))

        for tag in tags:
            db.execute("INSERT OR IGNORE INTO tags (tool_slug, tag) VALUES (?, ?)", (t['slug'], tag))

        if t['e1659a']:
            db.execute("""
                INSERT INTO e1659a_results (tool_slug, score, interpretation)
                VALUES (?, ?, ?)
            """, (t['slug'], t['e1659a']['score'], t['e1659a']['interpretation']))

    db.commit()

    count = db.execute("SELECT COUNT(*) FROM tools").fetchone()[0]
    e1659a_count = db.execute("SELECT COUNT(*) FROM e1659a_results").fetchone()[0]
    tag_count = db.execute("SELECT COUNT(DISTINCT tag) FROM tags").fetchone()[0]
    by_status = db.execute("SELECT status, COUNT(*) FROM tools GROUP BY status ORDER BY COUNT(*) DESC").fetchall()

    print(f"\nDatabase: {DB_FILE}")
    print(f"Tools: {count}, E1659A: {e1659a_count}, Tags: {tag_count}")
    print(f"Descriptions from Brain: {brain_found}/{count}")
    print(f"\nBy status:")
    for status, cnt in by_status:
        print(f"  {status}: {cnt}")

    db.close()


if __name__ == '__main__':
    main()
