#!/usr/bin/env python3
"""
Sync tool pages from Brain vault to site.
Brain = single source of truth.

Usage: python3 scripts/sync-tools-from-brain.py
Run before: npm run build && wrangler pages deploy dist ...
"""

import re
import os
from pathlib import Path

BRAIN = Path.home() / "Library/Mobile Documents/iCloud~md~obsidian/Documents/Brain"
SITE = Path(__file__).parent.parent
TOOLS_DIR = SITE / "src/pages/tools/tool"

# Tool slug → Brain note filename (without .md)
TOOLS = {
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
    # === Deep Research additions (107 tools) ===
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
    "eigen-eigen-pc": "Eigen / Eigen-PC",
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
    "ddgun-ddgun3d": "DDGun / DDGun3D",
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
    "exac-legacy": "ExAC (Legacy)",
    "alfa": "ALFA",
    "lovd": "LOVD",
    "biobank-japan-bbj-pheweb": "Biobank Japan (BBJ) / PheWeb",
    "cmdb-chinese-millionome": "CMDB (Chinese Millionome)",
    "kob-kdna": "KoB / KDNA",
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
    "db-oto-trial-resources": "DB-OTO Trial Resources",
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


def escape(text):
    """Escape for Astro template."""
    return (text
            .replace('&', '&amp;')
            .replace('<', '&lt;')
            .replace('>', '&gt;')
            .replace('{', '&#123;')
            .replace('}', '&#125;'))


def parse_frontmatter(content):
    """Extract YAML frontmatter values."""
    m = re.match(r'^---\s*\n(.*?)\n---', content, re.DOTALL)
    if not m:
        return {}
    fm = {}
    for line in m.group(1).split('\n'):
        if ':' in line:
            key, val = line.split(':', 1)
            fm[key.strip()] = val.strip()
    return fm


def parse_sections(content):
    """Split markdown into (title, body) sections, stopping at Connections."""
    # Strip frontmatter
    content = re.sub(r'^---.*?---\s*', '', content, flags=re.DOTALL)
    # Strip H1
    content = re.sub(r'^# .+\n+', '', content)

    sections = []
    cur_title = None
    cur_lines = []

    for line in content.split('\n'):
        if line.startswith('## '):
            if cur_title:
                sections.append((cur_title, '\n'.join(cur_lines)))
            title = line[3:].strip()
            if title == 'Connections':
                break
            cur_title = title
            cur_lines = []
        else:
            cur_lines.append(line)

    if cur_title and cur_lines:
        sections.append((cur_title, '\n'.join(cur_lines)))
    return sections


def md_to_html(text):
    """Minimal markdown → HTML for a section body."""
    lines = text.strip().split('\n')
    out = []
    in_code = False
    in_list = False

    for line in lines:
        # Code fences
        if line.startswith('```'):
            if in_code:
                out.append('</code></pre>')
                in_code = False
            else:
                out.append(
                    '<pre class="bg-black/40 border border-white/10 '
                    'rounded p-3 text-xs overflow-x-auto my-2"><code>')
                in_code = True
            continue

        if in_code:
            out.append(escape(line))
            continue

        # List items
        if line.startswith('- '):
            if not in_list:
                out.append('<ul class="space-y-1 my-2">')
                in_list = True
            item = _inline(line[2:])
            out.append(f'<li class="text-xs text-neutral-400">{item}</li>')
            continue

        # Close list if needed
        if in_list:
            out.append('</ul>')
            in_list = False

        # Skip table rows
        if line.startswith('|'):
            continue

        # Regular paragraph
        if line.strip():
            out.append(
                f'<p class="text-xs text-neutral-400 my-1">{_inline(line)}</p>')

    if in_list:
        out.append('</ul>')
    if in_code:
        out.append('</code></pre>')
    return '\n        '.join(out)


def _inline(text):
    """Convert inline markdown (bold, code, wikilinks)."""
    text = re.sub(r'\*\*(.+?)\*\*',
                  r'<strong class="text-neutral-200">\1</strong>', text)
    text = re.sub(r'\[\[(.+?)\]\]', r'\1', text)
    text = re.sub(r'`(.+?)`',
                  r'<code class="text-blue-300 text-xs">\1</code>', text)
    return text


def section_color(title):
    """Pick heading color by section type."""
    if 'Result' in title:
        return 'text-green-400'
    if 'Next' in title:
        return 'text-yellow-400'
    if 'Verified' in title or 'Status' in title:
        return 'text-blue-400'
    if 'Critical' in title or 'Limitation' in title:
        return 'text-red-400'
    if 'How to Use' in title:
        return 'text-purple-400'
    return 'text-neutral-300'


def build_page(note_name, content):
    """Generate an Astro page from a Brain note."""
    fm = parse_frontmatter(content)
    status = fm.get('status', 'unknown')
    tags_raw = fm.get('tags', '')
    # Parse [a, b, c] format
    tags = [t.strip() for t in tags_raw.strip('[]').split(',')
            if t.strip() and t.strip() not in ('tool', 'strc-research')]

    # First paragraph after H1 = description
    body = re.sub(r'^---.*?---\s*', '', content, flags=re.DOTALL)
    desc_m = re.search(r'^# .+\n+(.+?)(?:\n\n|\n##)', body, re.DOTALL)
    description = desc_m.group(1).strip()[:300] if desc_m else note_name

    # Sections → HTML
    sections = parse_sections(content)
    sections_html = []
    for title, body in sections:
        if not body.strip():
            continue
        html = md_to_html(body)
        color = section_color(title)
        sections_html.append(f'''
      <div class="mt-6">
        <h3 class="text-sm font-semibold {color} mb-2">{escape(title)}</h3>
        {html}
      </div>''')

    # Badge
    badges = {
        'verified': ('bg-green-500/20 text-green-400', 'VERIFIED'),
        'partial':  ('bg-yellow-500/20 text-yellow-400', 'PARTIAL'),
    }
    cls, label = badges.get(status, ('bg-neutral-500/20 text-neutral-400', 'DOCUMENTED'))
    badge = f'<span class="px-2 py-0.5 text-[10px] rounded {cls} font-semibold">{label}</span>'

    # Tag pills
    tag_html = '\n            '.join(
        f'<span class="text-[10px] px-1.5 py-0.5 rounded bg-white/5 '
        f'text-neutral-500">{t}</span>'
        for t in tags
    )

    return f'''---
import Layout from '../../../layouts/Layout.astro';
const lang = 'en';
---

<Layout lang={{lang}}>
  <div class="px-6 pt-6 pb-2 max-w-3xl mx-auto">
    <a href="/tools" class="text-xs text-neutral-500 hover:text-white transition-colors no-underline">&larr; All Tools</a>
  </div>

  <section class="px-6 pt-4 pb-16">
    <div class="max-w-3xl mx-auto">
      <div class="flex items-center gap-3 mb-3">
        <h1 class="text-lg font-bold text-white">{escape(note_name)}</h1>
        {badge}
      </div>

      <p class="text-sm text-neutral-300 leading-relaxed mb-4">{escape(description)}</p>

      <div class="flex flex-wrap gap-1 mb-6">
        {tag_html}
      </div>

      <div class="bg-white/[0.02] border border-white/5 rounded-lg p-5">
        {''.join(sections_html)}
      </div>
    </div>
  </section>
</Layout>
'''


def main():
    TOOLS_DIR.mkdir(parents=True, exist_ok=True)

    synced = 0
    skipped = 0

    for slug, note_name in TOOLS.items():
        note_path = BRAIN / "notes" / f"{note_name}.md"
        if not note_path.exists():
            print(f"  SKIP  {note_name} (not found)")
            skipped += 1
            continue

        content = note_path.read_text()
        page = build_page(note_name, content)

        out_path = TOOLS_DIR / f"{slug}.astro"
        out_path.write_text(page)
        synced += 1

    print(f"\nSynced {synced} tool pages from Brain → site")
    if skipped:
        print(f"Skipped {skipped} (iCloud dataless? run: open <file>)")


if __name__ == '__main__':
    main()
