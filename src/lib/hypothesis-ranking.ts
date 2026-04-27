import rankingJson from '../data/hypothesis-ranking.json';
import slugsJson from '../data/hypothesis-slugs.json';

export type Tier = 'S' | 'A' | 'B' | 'C' | 'D' | 'D/C' | 'reference';

export interface Hypothesis {
  num: number;
  title: string;
  slug: string | null;
  mech: string;
  deliv: string;
  misha: string;
  status: string;
  tier: Tier;
  nextStep: string;
}

const slugMap = (slugsJson as { slugs: Record<string, string | null> }).slugs;

function cleanCell(cell: string): string {
  return cell
    .replace(/\*\*/g, '')
    .replace(/\p{Extended_Pictographic}/gu, '')
    .trim();
}

function normalizeTier(raw: string): Tier {
  const t = cleanCell(raw);
  if (t === '—' || t === '-' || t === '') return 'reference';
  if (/\bD\s*\/\s*C\b|\bC\s*\/\s*D\b/.test(t)) return 'D/C';
  const m = t.match(/\b([SABCD])\b/);
  if (m) return m[1] as Tier;
  return 'reference';
}

export function parseRanking(): {
  hypotheses: Hypothesis[];
  updatedDate: string | null;
} {
  const data = rankingJson as {
    updatedDate: string | null;
    hypotheses: Array<Omit<Hypothesis, 'slug' | 'tier'> & { tier: string }>;
  };

  return {
    updatedDate: data.updatedDate,
    hypotheses: data.hypotheses.map((h) => ({
      ...h,
      slug: slugMap[h.title] ?? null,
      tier: normalizeTier(h.tier),
    })),
  };
}

export const TIER_ORDER: Tier[] = ['S', 'A', 'B', 'C', 'D/C', 'D', 'reference'];

export const TIER_META: Record<Tier, { label: string; caption: string; color: string }> = {
  S: {
    label: 'S-tier',
    caption: 'Primary — active compute now',
    color: 'bg-red-500/15 text-red-300 border-red-500/30',
  },
  A: {
    label: 'A-tier',
    caption: 'Active backburner — advances when S blocked',
    color: 'bg-amber-500/15 text-amber-300 border-amber-500/30',
  },
  B: {
    label: 'B-tier',
    caption: 'Watch — incremental engagement',
    color: 'bg-blue-500/15 text-blue-300 border-blue-500/30',
  },
  C: {
    label: 'C-tier',
    caption: 'Paused — awaiting external catalyst',
    color: 'bg-slate-500/15 text-slate-300 border-slate-500/30',
  },
  'D/C': {
    label: 'D/C',
    caption: 'Partially killed, partially paused',
    color: 'bg-neutral-500/15 text-neutral-400 border-neutral-500/30',
  },
  D: {
    label: 'D-tier',
    caption: 'Killed — mechanism falsified or wrong-patient',
    color: 'bg-neutral-700/30 text-neutral-500 border-neutral-700/40',
  },
  reference: {
    label: 'Reference',
    caption: 'Supporting models, engineering helpers, reference notes',
    color: 'bg-neutral-800/40 text-neutral-400 border-neutral-700/30',
  },
};
