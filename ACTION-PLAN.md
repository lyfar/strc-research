# STRC Site Tools Refactor — Action Plan

## Overview
Migrating 143 hardcoded tool pages to data-driven architecture.
Before: 143 individual .astro files, each hand-crafted or sync'd badly.
After: 1 data file + 1 template + 1 index page.

## Architecture

```
src/data/tools.ts          ← Single source of truth (all 143 tools)
src/pages/tools/tool/[slug].astro  ← Dynamic route (generates all tool pages)
src/pages/tools.astro      ← Index page with categories + search
src/components/ToolCard.astro      ← Reusable card component
src/components/StatusBadge.astro   ← Status badge component
scripts/sync-tools-from-brain.py   ← Updated: now generates tools.ts
scripts/cleanup-old-tools.sh       ← Removes old static pages
```

## Steps

### Phase 1: Data Layer ✅
- [x] Create `src/data/tools.ts` with all tool data from Brain vault
- [x] Define TypeScript interfaces (Tool, ToolStatus, ToolCategory)
- [x] Include E1659A test results for verified tools

### Phase 2: Components ✅
- [x] Create `src/components/StatusBadge.astro`
- [x] Create `src/components/ToolCard.astro`

### Phase 3: Pages ✅
- [x] Create `src/pages/tools/tool/[slug].astro` (dynamic route)
- [x] Rewrite `src/pages/tools.astro` (data-driven index)

### Phase 4: Test Build
- [ ] `cd /Users/egorlyfar/Sites/site-strc-egor-lol`
- [ ] `npm run build` — must succeed with 0 errors
- [ ] `npm run preview` — check at http://localhost:4321
- [ ] Verify: /tools page loads with all categories
- [ ] Verify: /tools/tool/alphamissense loads with correct data
- [ ] Verify: /tools/tool/dynamut shows DynaMut2 result (-0.913)
- [ ] Verify: search filters tools correctly
- [ ] Verify: broken tool pages show correct DOWN/BROKEN badges

### Phase 5: Cleanup
- [ ] Run `bash scripts/cleanup-old-tools.sh` (dry run first)
- [ ] Verify build still works after cleanup
- [ ] Run cleanup with --yes
- [ ] Final build + preview check

### Phase 6: Deploy
- [ ] `npm run build`
- [ ] `npx wrangler pages deploy dist --project-name=site-strc-egor-lol`
- [ ] Verify live at https://strc.egor.lol/tools
- [ ] Check 3-4 random tool pages on live site
- [ ] Verify no 404s on old URLs

### Phase 7: Localized Pages (later)
- [ ] Update zh/tools.astro, fr/tools.astro, ru/tools.astro, es/tools.astro, ja/tools.astro
- [ ] These should import from same tools.ts data
- [ ] Only translate UI chrome (headings, labels), not tool names

## Key URLs to test after deploy
- https://strc.egor.lol/tools — index with all 143 tools
- https://strc.egor.lol/tools/tool/alphamissense — flagship tool (verified, E1659A data)
- https://strc.egor.lol/tools/tool/dynamut — DynaMut2 result
- https://strc.egor.lol/tools/tool/gnomad — population database (verified)
- https://strc.egor.lol/tools/tool/mutscore — broken link example
- https://strc.egor.lol/tools/tool/sage — site-down example
- https://strc.egor.lol/tools/tool/evo-evo-2 — foundation model

## Sync from Brain (ongoing)
To update tools data from Brain vault:
```bash
cd /Users/egorlyfar/Sites/site-strc-egor-lol
python3 scripts/sync-tools-from-brain.py
npm run build
npx wrangler pages deploy dist --project-name=site-strc-egor-lol
```

## Rollback
If something breaks:
```bash
git checkout src/pages/tools.astro
git checkout src/pages/tools/tool/
git checkout src/data/tools.ts
```
