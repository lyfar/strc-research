#!/usr/bin/env python3
"""Generate src/data/tools.ts from SQLite database.

All data comes from DB: description, tags, status, url, e1659a.
DB is the single source of truth for the site.
Brain vault is upstream — sync via migrate-ts-to-db.py or tools-cli.py.
"""
import sqlite3, json
from pathlib import Path

ROOT = Path(__file__).parent.parent
DB_FILE = ROOT / "data/tools.db"
OUT = ROOT / "src/data/tools.ts"


def esc(s):
    return s.replace('\\', '\\\\').replace("'", "\\'").replace('\n', ' ')


def main():
    db = sqlite3.connect(str(DB_FILE))
    db.row_factory = sqlite3.Row

    categories = [r['name'] for r in db.execute(
        "SELECT name FROM categories ORDER BY sort_order"
    )]

    tools = db.execute("""
        SELECT t.*, c.sort_order
        FROM tools t
        JOIN categories c ON t.category = c.name
        ORDER BY c.sort_order, t.name
    """).fetchall()

    tags_map = {}
    for r in db.execute("SELECT tool_slug, tag FROM tags ORDER BY tool_slug, tag"):
        tags_map.setdefault(r['tool_slug'], []).append(r['tag'])

    e1659a_map = {}
    for r in db.execute("SELECT tool_slug, score, interpretation FROM e1659a_results"):
        e1659a_map[r['tool_slug']] = {'score': r['score'], 'interpretation': r['interpretation']}

    lines = []
    lines.append("// Auto-generated from tools.db. Do not edit manually.")
    lines.append("// Run: python3 scripts/db-to-ts.py")
    lines.append("")
    lines.append("export type ToolStatus = 'verified' | 'tested' | 'available' | 'stub' | 'site-down' | 'link-broken';")
    lines.append("")
    lines.append("export type ToolCategory =")
    for i, c in enumerate(categories):
        sep = ';' if i == len(categories) - 1 else ''
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
    lines.append(f"export const CATEGORIES: ToolCategory[] = {json.dumps(categories)};")
    lines.append("")
    lines.append("export const tools: Tool[] = [")

    for t in tools:
        tags = tags_map.get(t['slug'], [])
        e1659a = e1659a_map.get(t['slug'])
        desc = t['description'] if t['description'] else t['name']

        lines.append("  {")
        lines.append(f"    slug: '{t['slug']}',")
        lines.append(f"    name: '{esc(t['name'])}',")
        lines.append(f"    category: '{t['category']}',")
        lines.append(f"    description: '{esc(desc)}',")
        lines.append(f"    url: '{esc(t['url'])}',")
        lines.append(f"    status: '{t['status']}',")
        if e1659a:
            lines.append(f"    e1659a: {{ score: '{esc(e1659a['score'])}', interpretation: '{esc(e1659a['interpretation'])}' }},")
        lines.append(f"    tags: {json.dumps(tags)},")
        if t['acmg_relevance']:
            lines.append(f"    acmgRelevance: '{esc(t['acmg_relevance'])}',")
        lines.append(f"    isFree: {'true' if t['is_free'] else 'false'},")
        lines.append(f"    hasApi: {'true' if t['has_api'] else 'false'},")
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

    verified = sum(1 for t in tools if t['status'] == 'verified')
    tested = len(e1659a_map)
    print(f"Generated {OUT}")
    print(f"  {len(tools)} tools, {verified} verified, {tested} with E1659A data")

    db.close()


if __name__ == '__main__':
    main()
