import { readFileSync, readdirSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';

const vault = process.env.STRC_VAULT_PATH ?? '/Users/egorlyfar/STRC';
const hypothesesDir = join(vault, 'hypotheses');
const outputPath = new URL('../src/data/hypothesis-ranking.json', import.meta.url);

function unquote(value) {
  const trimmed = value.trim();
  if (trimmed.startsWith('"') && trimmed.endsWith('"')) {
    return trimmed.slice(1, -1).replace(/\\"/g, '"');
  }
  if (trimmed.startsWith("'") && trimmed.endsWith("'")) {
    return trimmed.slice(1, -1).replace(/\\'/g, "'");
  }
  return trimmed;
}

function parseFrontmatter(markdown) {
  const match = markdown.match(/^---\n([\s\S]*?)\n---/);
  if (!match) return {};

  const data = {};
  for (const line of match[1].split('\n')) {
    if (!line || /^\s/.test(line) || line.trim().startsWith('-')) continue;
    const sep = line.indexOf(':');
    if (sep === -1) continue;

    const key = line.slice(0, sep).trim();
    const value = line.slice(sep + 1).trim();
    if (value === '') continue;
    data[key] = unquote(value);
  }
  return data;
}

function clean(value) {
  return String(value ?? '')
    .replace(/\[\[([^\]|]+)\|([^\]]+)\]\]/g, '$2')
    .replace(/\[\[([^\]]+)\]\]/g, '$1')
    .replace(/`([^`]+)`/g, '$1')
    .replace(/\*\*/g, '')
    .replace(/\s+/g, ' ')
    .trim();
}

function shorten(value, max = 260) {
  const text = clean(value);
  if (text.length <= max) return text;
  return `${text.slice(0, max - 1).trim()}…`;
}

const hypotheses = readdirSync(hypothesesDir, { withFileTypes: true })
  .filter((entry) => entry.isDirectory() && /^h\d+/.test(entry.name))
  .map((entry) => {
    const markdown = readFileSync(join(hypothesesDir, entry.name, 'index.md'), 'utf8');
    const fm = parseFrontmatter(markdown);
    return {
      num: Number(fm.hypothesis_num),
      title: clean(fm.hypothesis_title || fm.title),
      mech: clean(fm.mech),
      deliv: clean(fm.deliv),
      misha: clean(fm.misha_fit),
      tier: clean(fm.tier),
      status: shorten(fm.excerpt, 220),
      nextStep: shorten(fm.next_step, 320),
    };
  })
  .filter((h) => Number.isFinite(h.num) && h.title)
  .sort((a, b) => a.num - b.num);

const generated = {
  updatedDate: new Date().toISOString().slice(0, 10),
  source: `${hypothesesDir}/*/index.md`,
  hypotheses,
};

writeFileSync(outputPath, `${JSON.stringify(generated, null, 2)}\n`);
console.log(`Wrote ${hypotheses.length} hypotheses to ${outputPath.pathname}`);
