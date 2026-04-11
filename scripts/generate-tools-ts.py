#!/usr/bin/env python3
"""Generate src/data/tools.ts from Brain vault notes."""
import re, json
from pathlib import Path

BRAIN = Path.home() / "Library/Mobile Documents/iCloud~md~obsidian/Documents/Brain/notes"
OUT = Path(__file__).parent.parent / "src/data/tools.ts"

# slug → Brain note name (from sync script)
TOOLS = {
    "alphafold-db": "AlphaFold Database", "alphafold-3": "AlphaFold 3 Server",
    "alphagenome": "AlphaGenome", "swiss-model": "SWISS-MODEL",
    "molstar": "Mol Star Viewer", "pymol": "PyMOL", "esmfold": "ESMFold",
    "iupred3": "IUPred3", "alphamissense": "AlphaMissense", "revel": "REVEL",
    "clinvar": "ClinVar", "spliceai": "SpliceAI", "cadd": "CADD",
    "varsome": "VarSome", "intervar": "InterVar", "franklin": "Franklin",
    "dbnsfp": "dbNSFP", "gnomad": "gnomAD", "seqr": "seqr",
    "decipher": "DECIPHER", "clingen": "ClinGen", "uniprot": "UniProt",
    "ucsc": "UCSC Genome Browser", "ensembl": "Ensembl REST API",
    "orthodb": "OrthoDB", "clustal-omega": "Clustal Omega", "omim": "OMIM",
    "dvd": "Deafness Variation Database", "hhl": "Hereditary Hearing Loss Homepage",
    "clinicaltrials": "ClinicalTrials.gov", "pubmed": "PubMed",
    "aav-capsid": "AAV Capsid Database", "asgct": "ASGCT Gene Therapy Database",
    "signalp": "SignalP", "netgpi": "NetGPI", "netnglyc": "NetNGlyc",
    "esm1v": "ESM1v",
    "eve-evolutionary-model-of-variant-effect": "EVE (Evolutionary model of Variant Effect)",
    "primateai-3d": "PrimateAI-3D", "mutscore": "MutScore", "gmvp": "gMVP",
    "metarnn": "MetaRNN", "bayesdel": "BayesDel", "mpc": "MPC", "dann": "DANN",
    "fathmm-mkl": "FATHMM-MKL", "linsight": "LINSIGHT",
    "eigen-eigen-pc": "Eigen", "fitcons": "fitCons",
    "clinpred": "ClinPred", "vest4": "VEST4", "snpeff": "SnpEff",
    "sparkinferno": "SparkINFERNO", "favor": "FAVOR",
    "rosettafold": "RoseTTAFold", "ifum": "IFUM", "foldx": "FoldX",
    "rosettaddg": "RosettaDDG", "ddgun-ddgun3d": "DDGun",
    "dynamut": "DynaMut", "mcsm": "mCSM", "sdm": "SDM", "duet": "DUET",
    "popmusic": "PoPMuSiC", "strum": "STRUM", "elaspic": "ELASPIC",
    "pangolin": "Pangolin", "absplice2": "AbSplice2", "maxentscan": "MaxEntScan",
    "genesplicer": "GeneSplicer", "nnsplice": "NNSplice", "spanr": "SPANR",
    "hal": "HAL", "spip": "SPiP", "splicerover": "SpliceRover",
    "mmsplice": "MMSplice", "fraser2": "FRASER2", "gtex-portal": "GTEx Portal",
    "encode": "ENCODE", "roadmap-epigenomics": "Roadmap Epigenomics",
    "regulomedb": "RegulomeDB", "deepsea": "DeepSEA", "basenji": "Basenji",
    "enformer": "Enformer", "tland": "TLand", "haploreg": "HaploReg",
    "piq": "PIQ", "bravo-topmed": "BRAVO / TOPMed", "uk-biobank": "UK Biobank",
    "dbsnp": "dbSNP", "exac-legacy": "ExAC", "alfa": "ALFA", "lovd": "LOVD",
    "biobank-japan-bbj-pheweb": "Biobank Japan",
    "cmdb-chinese-millionome": "CMDB", "kob-kdna": "KoB KDNA", "sage": "SAGE",
    "hgmd-professional": "HGMD Professional", "mastermind": "Mastermind",
    "civic": "CIViC", "pharmgkb": "PharmGKB", "dgidb": "DGIdb",
    "gwas-catalog": "GWAS Catalog", "open-targets": "Open Targets", "cgar": "CGAR",
    "human-protein-atlas": "Human Protein Atlas",
    "expression-atlas": "Expression Atlas", "allen-brain-atlas": "Allen Brain Atlas",
    "archs4": "ARCHS4", "audiogene": "AudioGene", "hiedra": "HIEDRA",
    "gear": "gEAR", "sensorion-aav-platform": "Sensorion AAV Platform",
    "db-oto-trial-resources": "DB-OTO",
    "consurf": "ConSurf", "evcouplings": "EVcouplings", "rate4site": "Rate4Site",
    "classifycnv": "ClassifyCNV", "annotsv": "AnnotSV", "cnvnator": "CNVnator",
    "gatk-gcnv": "GATK gCNV", "manta": "Manta", "delly": "Delly",
    "cnest": "CNest", "variantvalidator": "VariantValidator",
    "mutalyzer": "Mutalyzer", "crossmap": "CrossMap", "litvar-20": "LitVar 2.0",
    "pubtator": "PubTator", "tmvar": "tmVar", "europe-pmc": "Europe PMC",
    "semantic-scholar": "Semantic Scholar", "litsuggest": "LitSuggest",
    "litsense": "LitSense", "galaxy-project": "Galaxy Project",
    "terrabio": "Terra.bio",
    "illumina-connected-insights": "Illumina Connected Insights",
    "variant-tools-vtools": "variant tools (vtools)",
    "evo-evo-2": "Evo / Evo 2",
    "nucleotide-transformer-ntv3": "Nucleotide Transformer (NTv3)",
    "caduceus": "Caduceus", "hyenadna": "HyenaDNA", "dnabert-2": "DNABERT-2",
}

# Category mapping from tags
TAG_TO_CAT = {
    'variant-effect-prediction': 'Variant Effect Prediction',
    'variant-classification': 'Variant Effect Prediction',
    'pathogenicity': 'Variant Effect Prediction',
    'structural-biology': 'Structural Biology',
    'protein-structure': 'Structural Biology',
    'splicing-prediction': 'Splicing Prediction',
    'regulatory-genomics': 'Regulatory & Non-Coding',
    'population-genetics': 'Population Databases',
    'clinical-genomics': 'Clinical Databases',
    'clinical': 'Clinical Databases',
    'gene-level': 'Gene-Level Resources',
    'hearing-loss': 'Hearing Loss & Inner Ear',
    'conservation': 'Conservation & Evolution',
    'cnv-analysis': 'Structural Variants & CNV',
    'nomenclature': 'Nomenclature & Validation',
    'literature-mining': 'Literature Mining',
    'workflow': 'Workflow Platforms',
    'emerging-ai': 'AI & Foundation Models',
    'foundation-model': 'AI & Foundation Models',
}

# Manual category overrides for pre-research tools
MANUAL_CAT = {
    'alphafold-db': 'Structural Biology', 'alphafold-3': 'Structural Biology',
    'swiss-model': 'Structural Biology', 'molstar': 'Structural Biology',
    'pymol': 'Structural Biology', 'esmfold': 'Structural Biology',
    'iupred3': 'Structural Biology', 'signalp': 'Structural Biology',
    'netgpi': 'Structural Biology', 'netnglyc': 'Structural Biology',
    'alphamissense': 'Variant Effect Prediction', 'revel': 'Variant Effect Prediction',
    'cadd': 'Variant Effect Prediction', 'spliceai': 'Splicing Prediction',
    'dbnsfp': 'Variant Effect Prediction',
    'varsome': 'Variant Effect Prediction', 'intervar': 'Variant Effect Prediction',
    'franklin': 'Variant Effect Prediction',
    'clinvar': 'Clinical Databases', 'omim': 'Clinical Databases',
    'clinicaltrials': 'Clinical Databases', 'pubmed': 'Literature Mining',
    'dvd': 'Hearing Loss & Inner Ear', 'hhl': 'Hearing Loss & Inner Ear',
    'aav-capsid': 'Hearing Loss & Inner Ear', 'asgct': 'Hearing Loss & Inner Ear',
    'gnomad': 'Population Databases', 'seqr': 'Population Databases',
    'decipher': 'Clinical Databases', 'clingen': 'Population Databases',
    'uniprot': 'Gene-Level Resources', 'ucsc': 'Gene-Level Resources',
    'ensembl': 'Gene-Level Resources', 'dbsnp': 'Population Databases',
    'orthodb': 'Conservation & Evolution', 'clustal-omega': 'Conservation & Evolution',
    'alphagenome': 'AI & Foundation Models',
}

# E1659A results
E1659A = {
    'alphamissense': ('0.9016', 'Likely Pathogenic'),
    'revel': ('0.789', 'Pathogenic range'),
    'cadd': ('PHRED 25.5', 'Top 0.3% deleterious'),
    'dann': ('0.9946', 'Highly deleterious'),
    'gnomad': ('Not found', 'Absent from 251K controls (PM2)'),
    'clinvar': ('N/A', 'Not yet submitted'),
    'bayesdel': ('0.2255', 'Damaging'),
    'metarnn': ('0.8552', 'Damaging'),
    'clinpred': ('0.9869', 'Very strong pathogenic signal'),
    'vest4': ('0.5900', 'Moderate pathogenic signal'),
    'fathmm-mkl': ('0.9748', 'Strong damaging'),
    'dynamut': ('-0.913 kcal/mol', 'Destabilizing'),
    'variantvalidator': ('Valid', 'HGVS confirmed, MANE select'),
    'mutalyzer': ('Valid', 'NM_153700.2:c.4976A>C confirmed'),
    'regulomedb': ('N/A', 'No regulatory variants (coding region)'),
    'open-targets': ('0.731', '73 STRC disease associations'),
    'litvar-20': ('N/A', 'E1659A not found in literature'),
    'pubtator': ('N/A', 'E1659A not found'),
    'semantic-scholar': ('13 papers', 'None mention E1659A'),
    'dbsnp': ('N/A', 'No rsID assigned'),
    'lovd': ('N/A', 'STRC variants present, E1659A absent'),
    'bravo-topmed': ('Not found', 'Absent from TOPMed'),
    'alphafold-db': ('pLDDT 68.75', 'Moderate confidence at E1659'),
    'alphagenome': ('54,276 scores', 'Splice quantile 0.997+'),
    'spliceai': ('Low', 'Missense, not splice-disrupting'),
    'varsome': ('VUS', 'PM2_Supporting + BP1'),
    'dbnsfp': ('40 scores', 'All extracted via myvariant.info'),
}

STATUS_MAP = {
    'verified': 'verified', 'partial': 'tested', 'documented': 'available',
    'stub': 'stub', 'site-down': 'site-down', 'link-broken': 'link-broken',
    'link-fixed': 'available', 'upcoming': 'stub', 'unknown': 'stub',
}

def parse_note(path):
    """Extract frontmatter and first paragraph."""
    try:
        content = path.read_text()
    except:
        return None
    # Frontmatter
    fm = {}
    m = re.match(r'^---\s*\n(.*?)\n---', content, re.DOTALL)
    if m:
        for line in m.group(1).split('\n'):
            if ':' in line:
                k, v = line.split(':', 1)
                fm[k.strip()] = v.strip()
    # Tags
    tags_raw = fm.get('tags', '')
    tags = [t.strip() for t in tags_raw.strip('[]').split(',') if t.strip()]
    # Status
    status = fm.get('status', 'unknown')
    # Description from first paragraph or "Why It Matters"
    body = re.sub(r'^---.*?---\s*', '', content, re.DOTALL)
    body = re.sub(r'^# .+\n+', '', body)
    desc = ''
    # Try "Why It Matters" section
    wm = re.search(r'## Why It Matters\s*\n+(.+?)(?:\n\n|\n##)', body, re.DOTALL)
    if wm:
        desc = wm.group(1).strip().split('\n')[0]
    if not desc:
        # First non-empty paragraph
        for line in body.split('\n'):
            line = line.strip()
            if line and not line.startswith('#') and not line.startswith('-') and not line.startswith('|') and not line.startswith('*'):
                desc = line
                break
    # URL from Key Info
    url = ''
    um = re.search(r'\*\*URL\*\*[:\s]*(?:~~[^~]+~~\s*)?(?:BROKEN\s*[—-]\s*)?(.+?)(?:\s*\(|$)', body, re.MULTILINE)
    if um:
        url = um.group(1).strip()
        # Clean markdown links
        url = re.sub(r'\[([^\]]+)\]\([^\)]+\)', r'\1', url)
    if not url:
        um2 = re.search(r'https?://[^\s<>\)\]]+', body)
        if um2:
            url = um2.group(0)
    # API
    has_api = bool(re.search(r'\*\*API\*\*[:\s]*Yes', body, re.IGNORECASE))
    is_free = bool(re.search(r'\*\*Free\*\*[:\s]*Free', body, re.IGNORECASE))
    # ACMG
    acmg = ''
    am = re.search(r'\*\*ACMG[^*]*\*\*[:\s]*(.+)', body)
    if am:
        acmg = am.group(1).strip()
    return {
        'status': status, 'tags': tags, 'description': desc[:200],
        'url': url, 'has_api': has_api, 'is_free': is_free, 'acmg': acmg,
    }

def get_category(slug, tags):
    if slug in MANUAL_CAT:
        return MANUAL_CAT[slug]
    for tag in tags:
        if tag in TAG_TO_CAT:
            return TAG_TO_CAT[tag]
    return 'Variant Effect Prediction'  # fallback

def esc(s):
    return s.replace('\\', '\\\\').replace("'", "\\'").replace('\n', ' ')

def main():
    entries = []
    for slug, note_name in TOOLS.items():
        # Try multiple note name variants
        path = BRAIN / f"{note_name}.md"
        if not path.exists():
            # Try alternate names
            for alt in [note_name.replace(' / ', ' '), note_name.split(' (')[0]]:
                p2 = BRAIN / f"{alt}.md"
                if p2.exists():
                    path = p2
                    break
        data = parse_note(path) if path.exists() else None
        if not data:
            data = {'status': 'stub', 'tags': [], 'description': note_name,
                    'url': '', 'has_api': False, 'is_free': True, 'acmg': ''}

        status = STATUS_MAP.get(data['status'], 'stub')
        cat = get_category(slug, data['tags'])
        e1659a = E1659A.get(slug)

        entry = {
            'slug': slug, 'name': note_name, 'category': cat,
            'description': data['description'], 'url': data['url'],
            'status': status, 'tags': [t for t in data['tags'] if t not in ('tool', 'strc-research')],
            'acmg': data['acmg'], 'isFree': data['is_free'], 'hasApi': data['has_api'],
            'e1659a': e1659a,
        }
        entries.append(entry)

    # Generate TypeScript
    lines = []
    lines.append("// Auto-generated from Brain vault. Do not edit manually.")
    lines.append("// Run: python3 scripts/generate-tools-ts.py")
    lines.append("")
    lines.append("export type ToolStatus = 'verified' | 'tested' | 'available' | 'stub' | 'site-down' | 'link-broken';")
    lines.append("")
    lines.append("export type ToolCategory =")
    cats = ['AI & Foundation Models', 'Structural Biology', 'Variant Effect Prediction',
            'Splicing Prediction', 'Regulatory & Non-Coding', 'Population Databases',
            'Clinical Databases', 'Gene-Level Resources', 'Hearing Loss & Inner Ear',
            'Conservation & Evolution', 'Structural Variants & CNV',
            'Nomenclature & Validation', 'Literature Mining', 'Workflow Platforms']
    for i, c in enumerate(cats):
        sep = ';' if i == len(cats) - 1 else ''
        lines.append(f"  | '{c}'{sep}")
    lines.append("")
    lines.append("export interface Tool {")
    lines.append("  slug: string;")
    lines.append("  name: string;")
    lines.append("  category: ToolCategory;")
    lines.append("  description: string;")
    lines.append("  url: string;")
    lines.append("  status: ToolStatus;")
    lines.append("  e1659a?: { score: string; interpretation: string };")
    lines.append("  tags: string[];")
    lines.append("  acmgRelevance?: string;")
    lines.append("  isFree: boolean;")
    lines.append("  hasApi: boolean;")
    lines.append("}")
    lines.append("")
    lines.append("export const STATUS_CONFIG: Record<ToolStatus, { label: string; bg: string; text: string }> = {")
    lines.append("  verified:      { label: 'VERIFIED',  bg: 'bg-green-500/20',   text: 'text-green-400' },")
    lines.append("  tested:        { label: 'TESTED',    bg: 'bg-green-500/20',   text: 'text-green-400' },")
    lines.append("  available:     { label: 'AVAILABLE', bg: 'bg-blue-500/20',    text: 'text-blue-400' },")
    lines.append("  stub:          { label: 'STUB',      bg: 'bg-neutral-500/20', text: 'text-neutral-400' },")
    lines.append("  'site-down':   { label: 'DOWN',      bg: 'bg-red-500/20',     text: 'text-red-400' },")
    lines.append("  'link-broken': { label: 'BROKEN',    bg: 'bg-red-500/20',     text: 'text-red-400' },")
    lines.append("};")
    lines.append("")
    lines.append(f"export const CATEGORIES: ToolCategory[] = {json.dumps(cats)};")
    lines.append("")
    lines.append("export const tools: Tool[] = [")

    for e in entries:
        lines.append("  {")
        lines.append(f"    slug: '{e['slug']}',")
        lines.append(f"    name: '{esc(e['name'])}',")
        lines.append(f"    category: '{e['category']}',")
        lines.append(f"    description: '{esc(e['description'])}',")
        lines.append(f"    url: '{esc(e['url'])}',")
        lines.append(f"    status: '{e['status']}',")
        if e['e1659a']:
            s, i = e['e1659a']
            lines.append(f"    e1659a: {{ score: '{esc(s)}', interpretation: '{esc(i)}' }},")
        tags_str = json.dumps(e['tags'])
        lines.append(f"    tags: {tags_str},")
        if e['acmg'] and e['acmg'] != 'N/A':
            lines.append(f"    acmgRelevance: '{esc(e['acmg'])}',")
        lines.append(f"    isFree: {'true' if e['isFree'] else 'false'},")
        lines.append(f"    hasApi: {'true' if e['hasApi'] else 'false'},")
        lines.append("  },")

    lines.append("];")
    lines.append("")
    lines.append("export const verifiedCount = tools.filter(t => t.status === 'verified').length;")
    lines.append("export const testedCount = tools.filter(t => t.e1659a !== undefined).length;")
    lines.append("")
    lines.append("export function getToolsByCategory(): Map<ToolCategory, Tool[]> {")
    lines.append("  const map = new Map<ToolCategory, Tool[]>();")
    lines.append("  for (const cat of CATEGORIES) map.set(cat, []);")
    lines.append("  for (const tool of tools) map.get(tool.category)?.push(tool);")
    lines.append("  return map;")
    lines.append("}")
    lines.append("")
    lines.append("export function getToolBySlug(slug: string): Tool | undefined {")
    lines.append("  return tools.find(t => t.slug === slug);")
    lines.append("}")

    OUT.write_text('\n'.join(lines))
    print(f"Generated {OUT} with {len(entries)} tools")
    # Stats
    by_status = {}
    by_cat = {}
    for e in entries:
        by_status[e['status']] = by_status.get(e['status'], 0) + 1
        by_cat[e['category']] = by_cat.get(e['category'], 0) + 1
    print(f"  By status: {dict(sorted(by_status.items()))}")
    print(f"  By category: {dict(sorted(by_cat.items()))}")
    print(f"  With E1659A data: {sum(1 for e in entries if e['e1659a'])}")

if __name__ == '__main__':
    main()
