# STRC Research Monitor — System Architecture

This document describes how the research-monitor pipeline actually works end-to-end: what lives where, which symlinks point which direction, who writes what, and how changes propagate. Read this before making structural changes.

Living outside the code (in the agent's head) this was getting lossy. Writing it down so we don't lose it.

---

## What the system does

1. **Egor writes** hypothesis ranking + research context in Obsidian (Brain vault).
2. **A local file watcher** (`launchd` agent) detects changes and regenerates `PROBLEM.md` at the root of this repo — a single agent-facing context file — then commits and pushes.
3. **A scheduled cloud agent** (Anthropic "routine" / RemoteTrigger) wakes up daily. It clones this repo, reads `PROBLEM.md` to orient itself, scans PubMed + bioRxiv for new papers, deduplicates against existing ones, writes relevant ones into `papers/`, updates matching files in `hypotheses/` with synthesis notes, commits, and pushes.
4. **GitHub Actions** builds the Astro site on every push and deploys to Cloudflare Pages at `strc.egor.lol`.
5. **Egor reads** the result on the site (`/hypotheses`, `/journal`) and in Obsidian (papers and models are symlinked back into Brain).

The whole system is an autonomous loop. Human touches the vault, everything else updates without manual action.

---

## Data flow

```
Egor edits in Obsidian                      Cloud agent (daily 00:47 HKT)
        │                                             │
        ▼                                             ▼
~/Brain/notes/STRC Hypothesis Ranking.md      RemoteTrigger
~/Brain/research/strc/_context.md             (trig_01XfQV8QaaR5fjLPczpJJXVu)
~/Brain/research/strc/site-hypothesis-slugs.json     │
        │                                             │
        │ (launchd WatchPaths)                        │ clones repo
        ▼                                             ▼
~/bin/strc-sync-problem                       PROBLEM.md (read)
  ├─ assembles PROBLEM.md                     papers/ (read for dedup)
  ├─ git commit                               hypotheses/ (read for synthesis targets)
  └─ git push                                         │
        │                                             ▼
        ▼                                      PubMed + bioRxiv API
github.com/lyfar/strc-research (master)              │
        ▲                                             ▼
        │                               New/updated paper .md files
        └───────── git push ◄──────────  Synthesis .md files
        │
        ▼
GitHub Actions → Cloudflare Pages build
        │
        ▼
strc.egor.lol (site)
```

---

## The symlink topology (most important part)

This system has two git repos (Brain + site) that need to share content. The direction of symlinks matters because **cloud runners can't resolve symlinks that point to local paths**.

### Rule

> Content canonical in Brain, symlinked into site repo: for content Egor authors and cloud agents only *read*.
>
> Content canonical in site repo, symlinked into Brain: for data that cloud agents *write* (they need real files at clone time).

### Current state

| Path in Brain | Path in site repo | Canonical side | Why |
|---|---|---|---|
| `~/Brain/notes/STRC Hypothesis Ranking.md` | (not in repo) | Brain | Egor's priority register, authored in Obsidian. Read by sync script. |
| `~/Brain/research/strc/_context.md` | (not in repo) | Brain | Problem context, authored by Egor. Read by sync script. |
| `~/Brain/research/strc/site-hypothesis-slugs.json` | `src/data/hypothesis-slugs.json` (copy) | Brain | Edited rarely; re-copied into repo by sync for site build + by sync for `PROBLEM.md` generation. |
| `~/Brain/research/strc/papers` → symlink | `papers/` (real dir) | **Site repo** | Cloud agent writes here. Must be real files when cloud clones. |
| `~/Brain/research/strc/models` → symlink | `models/` (real dir) | **Site repo** | Same — cloud agents produce CIFs, PNGs, analysis outputs. |

Brain's `research/strc/papers` and `.../models` are symlinks into `~/Sites/site-strc-egor-lol/papers` and `.../models`. Obsidian follows symlinks transparently, so from Egor's side nothing changed.

### Memory / backup implications

- Brain lives in iCloud (`~/Brain` → `~/Library/Mobile Documents/iCloud~md~obsidian/Documents/Brain`). Brain files are synced to iCloud.
- Papers + models live in `~/Sites/site-strc-egor-lol/` which is NOT in iCloud but IS in git + pushed to GitHub. Net: papers+models backed up via git, not via iCloud. That's fine — git is a better backup than iCloud anyway.
- `.gitignore` excludes `models/**/__pycache__`, `*.png`, `*.jpg`, `.DS_Store`, and a few other build artifacts. Python scripts (`.py`), JSON results, PDBs, FASTAs, CIFs, and other computational outputs ARE tracked.

---

## PROBLEM.md — the agent's context file

**Location:** `PROBLEM.md` at repo root.

**Purpose:** single-file context the cloud agent reads on every run. Without it, every routine run would need its context hardcoded in the prompt, and the prompt would drift from Egor's priorities.

**Generator:** `~/bin/strc-sync-problem`. Bash script. Reads three Brain files, assembles PROBLEM.md, does `git commit` + `git push` if content changed.

Sections of PROBLEM.md:

1. **Core problem** — STRC biology, Misha's genotype, engineering constraints. Copied from `_context.md`.
2. **Active hypothesis register** — 25-row markdown table with M/D/F scoring, status, tier, evidence, next step. Copied verbatim from `STRC Hypothesis Ranking.md`.
3. **Hypothesis → synthesis file slug table** — mapping of hypothesis wikilinks to `hypotheses/<slug>.md`. Generated from `site-hypothesis-slugs.json`.
4. **Kill list** — D-tier hypotheses. Copied from ranking note.
5. **Agent directive** — how to use the file. Static text in the generator script.

**Regeneration trigger:** `launchd` agent (`~/Library/LaunchAgents/com.egor.strc-sync.plist`) watches:
- `STRC Hypothesis Ranking.md`
- `_context.md`
- `site-hypothesis-slugs.json`

On any change → sleeps 10s (to let Obsidian finish its atomic-write) → runs `strc-sync-problem` → commits + pushes automatically.

**Maintenance:** the system is built to require zero maintenance of `PROBLEM.md`. Egor never edits it directly. If he wants to change the agent's orientation, he edits the upstream Brain files.

---

## The cloud routine

**ID:** `trig_01XfQV8QaaR5fjLPczpJJXVu`
**Name:** "STRC Research Monitor"
**Manage:** https://claude.ai/code/scheduled/trig_01XfQV8QaaR5fjLPczpJJXVu
**Schedule:** `47 16 * * *` UTC → 00:47 Asia/Hong_Kong daily
**Model:** claude-opus-4-6
**Repo:** `github.com/lyfar/strc-research`

**Prompt flow (10 steps):**

- Step 0: `git pull`, `cat PROBLEM.md`, `ls hypotheses/`. Agent orients from PROBLEM.md.
- Step 1: Build four dedup lookup sets from existing `papers/*.md` — DOI, PMID, bioRxiv ID, and normalized-title hash (via inline Python).
- Step 2: Direct PubMed scan (last 2 days, STRC-adjacent queries + S-tier-specific queries).
- Step 3: Lateral PubMed scan (last 14 days, cross-disciplinary bridges).
- Step 4: bioRxiv scan (last 7 days, filtered by title/abstract keywords).
- Step 5: Batch-fetch abstracts for unique PMIDs; pre-filter against existing PMIDs/biorxiv set.
- Step 6: Score relevance, classify candidate, write. Three classes:
  - `skip-dup-*`: DOI/PMID/bioRxiv match → no-op.
  - `update-in-place`: title matches existing file but IDs differ → this is a preprint→published transition → update the existing file (preserving `date_added`), don't create a duplicate.
  - `write-new`: fresh paper → write new file.
  - Also detects `publication_status: ahead-of-print | published` by comparing `date` to today.
- Step 7: `git commit` + `git push` of `papers/`.
- Step 8: Match each new/updated paper to hypotheses via `PROBLEM.md` slug table.
- Step 9: Append "Recent Papers" entries to `hypotheses/<slug>.md` files. For hypotheses with `null` slug → skip (don't auto-create new files).
- Step 10: `git commit` + `git push` of `hypotheses/`.

**Management API:** `RemoteTrigger` tool (via `/schedule` skill). Actions: list, get, create, update, run, delete (delete must be done in the web UI).

---

## CF Pages deployment

**Workflow:** `.github/workflows/deploy.yml`. Triggered on push to `master`. Runs `npm ci && npm run build && wrangler pages deploy dist --project-name=strc-egor-lol --branch=main`.

**Why `--branch=main` when the git branch is `master`:** CF Pages was set up with a `main` branch, we kept the git branch as `master`. This mismatch is fine but confusing. Don't "fix" it without coordinating.

**Build time:** ~40-60 seconds from push to live on strc.egor.lol.

**Schema:** `src/content.config.ts` defines the `papers` content collection schema. If you add a frontmatter field to papers written by the agent, also add it here (optional) or build will fail.

---

## The site pages

- **`/hypotheses`** (`src/pages/hypotheses/index.astro`) — live priority register. Parses `src/data/hypothesis-ranking.md` (copy of Brain's ranking note) via `src/lib/hypothesis-ranking.ts`, groups by tier, joins with `hypothesis-slugs.json` to show links only where a detail page exists.
- **`/hypotheses/<slug>`** (individual files under `src/pages/hypotheses/`) — hand-authored long-form articles per hypothesis. Only 11 of the 25 register entries have these. Rest of the ranking shows as cards with "No detail page yet" on the index.
- **`/journal`** (`src/pages/journal.astro`) — research monitor. Sorts papers by `date_added` (agent run date) then by publication date. Shows S/A-tier priority context at the top so visitors see the scoring frame. Cards render "Ahead of print" badge for papers with future `date`.

---

## File inventory — what each file does

### In site repo
- `PROBLEM.md` — agent context, auto-generated by sync script.
- `papers/*.md` — paper summaries, written by cloud routine.
- `models/*` — computational outputs (CIFs, scripts, JSON). Egor produces these in other sessions; committed here. Accessible from Brain via symlink.
- `hypotheses/<slug>.md` — auto-maintained synthesis files (cloud routine appends "Recent Papers"). Can also be hand-edited.
- `src/pages/hypotheses/*.astro` — long-form hypothesis articles (hand-authored essays).
- `src/data/hypothesis-ranking.md` — copy of Brain's ranking note. Read by site build. Regenerated by sync script? **No — currently only `PROBLEM.md` is regenerated.** See "gap" below.
- `src/data/hypothesis-slugs.json` — copy of Brain's slug map. Same — currently hand-synced.
- `src/lib/hypothesis-ranking.ts` — parser for the ranking table. Produces typed objects for Astro pages.
- `.github/workflows/deploy.yml` — CF Pages deploy.

### Outside the repo
- `~/bin/strc-sync-problem` — the sync script.
- `~/Library/LaunchAgents/com.egor.strc-sync.plist` — launchd watcher config.
- `~/Library/Logs/strc-sync.log` — sync run log.
- `~/Brain/notes/STRC Hypothesis Ranking.md` — ranking (Brain-canonical).
- `~/Brain/research/strc/_context.md` — problem context (Brain-canonical).
- `~/Brain/research/strc/site-hypothesis-slugs.json` — slug map (Brain-canonical).
- `~/Brain/research/strc/papers` — symlink → `~/Sites/site-strc-egor-lol/papers`.
- `~/Brain/research/strc/models` — symlink → `~/Sites/site-strc-egor-lol/models`.

---

## Sync fan-out

`strc-sync-problem` writes three files into the site repo on every run:

1. `PROBLEM.md` — assembled from `_context.md` + ranking table + slug map + kill list.
2. `src/data/hypothesis-ranking.md` — direct `cp -L` of the Brain ranking note (site's `/hypotheses` page parses this).
3. `src/data/hypothesis-slugs.json` — direct `cp -L` of the slug map (same).

All three get committed in a single commit (`sync: refresh PROBLEM.md + ranking mirror from Brain`) and pushed together. So `/hypotheses` page and the cloud agent's context never diverge — both update on the same push.

---

## Maintenance runbook

### Adding a new hypothesis

1. Add a row to `~/Brain/notes/STRC Hypothesis Ranking.md` (new wikilink + tier + evidence + next step).
2. If you want papers to get synthesized into a detail file: add an entry to `~/Brain/research/strc/site-hypothesis-slugs.json` mapping the wikilink → slug (e.g., `"STRC New Hypothesis": "new-hyp"`). Use `null` if you just want the hypothesis visible in the register but no synthesis file.
3. Save. Watcher fires → `PROBLEM.md` regenerates → push → cloud agent sees the new hypothesis on next run.

### Changing a tier / killing a hypothesis

1. Edit the ranking row in `STRC Hypothesis Ranking.md`.
2. Save. Watcher syncs.
3. Next cloud run: agent weights relevance per new tier, logs any D-tier papers that would need flip-evidence.

### Swapping the deploy key / git auth

The watcher uses `osxkeychain` credential helper (inherited from user's git config). If git push fails for the watcher:
1. Check `~/Library/Logs/strc-sync.log` for `WARN: push failed`.
2. Run `gh auth status` — verify `lyfar` account is active.
3. If account switched: `gh auth switch` or `git config --global credential.helper osxkeychain` and re-auth.
4. Alternative: use an SSH deploy key scoped only to this repo. Set `GIT_SSH_COMMAND` in the sync script: `GIT_SSH_COMMAND='ssh -i ~/.ssh/strc-research_ed25519' git push`.

### Inspecting / running the cloud routine

- List/manage: `https://claude.ai/code/scheduled/trig_01XfQV8QaaR5fjLPczpJJXVu`
- Trigger manually: use `/schedule` skill → "run trigger"
- Change schedule: edit `cron_expression` via RemoteTrigger update
- Enable/disable: toggle `enabled: true|false`

### Changing the agent's orientation

Edit the prompt via `RemoteTrigger update` with a new `events[0].data.message.content`. The current prompt is large (~12 KB) — keep the 10-step structure and only touch specific steps.

Or: for structural guidance, edit `PROBLEM.md`'s generator (`~/bin/strc-sync-problem`) to change what gets assembled. That propagates on next watcher fire.

### Disabling / uninstalling the watcher

```bash
launchctl bootout gui/$(id -u)/com.egor.strc-sync
rm ~/Library/LaunchAgents/com.egor.strc-sync.plist
```

Re-install:

```bash
launchctl bootstrap gui/$(id -u) ~/Library/LaunchAgents/com.egor.strc-sync.plist
```

---

## Troubleshooting

**Cloud routine push fails with non-fast-forward:** means the routine's local clone fell behind. Happens if the watcher committed between the routine's clone and push. Routine doesn't handle this gracefully — currently exits. Fix is to add a `git pull --rebase` retry in the routine prompt. Hasn't happened enough to prioritize.

**`/hypotheses` page shows stale ranking:** see "Known gap" above. Manually `cp` the ranking file until auto-sync is extended.

**Build fails on CF Pages with `InvalidContentEntryDataError`:** schema validation. Either the routine wrote a paper with a new field not in `content.config.ts`, or a field the schema expects is missing. Fix: add the field as `.optional()` in schema, push.

**Watcher doesn't fire:** `launchctl list | grep strc-sync` — if not listed, reload with `launchctl bootstrap`. iCloud occasionally delays mtime events on synced files; if this becomes a problem, switch the watcher to `StartInterval` polling (every 5 minutes) instead of `WatchPaths`.

**Routine runs but makes 0 new papers:** normal if PubMed has no new STRC-adjacent work in the last 2 days. Check the commit message — it says `papers: YYYY-MM-DD scan — 0 direct, 0 lateral`. Nothing to worry about.

**Duplicate papers appearing:** dedup check failed somehow. Inspect `/tmp/existing_*.txt` files from the most recent routine run. Most likely cause is frontmatter formatting inconsistency (stray whitespace, quoting). The dedup grep is defensive but not exhaustive.

---

## Evolution history

Sessions where this system took shape, for future context recovery:

- **2026-04-16**: Research monitor page shipped. Papers journal wired with Astro content collections.
- **2026-04-17**: Routine prompt strengthened with PMID dedup. `date_added` field added for "New" badge. Schema divergence: `pubmed_id` made optional, `biorxiv_id` added.
- **2026-04-21**: Large architectural refactor — `/hypotheses` became live register from Brain ranking; `PROBLEM.md` + hypothesis-to-slug system introduced; symlink direction flipped for `papers/` and `models/`; sync script + launchd watcher added; four-key dedup + ahead-of-print handling added to routine prompt; this architecture doc written.

---

## Related documents

- `~/Brain/AGENTS.md` — Brain vault rules (note conventions, frontmatter requirements, wikilink standards).
- `~/Brain/research/strc/_context.md` — STRC-specific research context.
- `~/Brain/notes/STRC Hypothesis Ranking.md` — priority register.
- `~/.claude/CLAUDE.md` — behavioral guidelines for Claude when working with Brain.
- `~/.claude/projects/-Users-egorlyfar/CLAUDE.md` — identity / voice guidelines.

If this document falls out of sync with reality, treat the working system as source of truth, fix the doc, and commit.
