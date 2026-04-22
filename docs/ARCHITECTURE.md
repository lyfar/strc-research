# STRC Research Monitor — System Architecture

How the research-monitor pipeline works end-to-end after the 2026-04-22 Brain-canonical migration: what lives where, who writes what, how changes propagate. Read this before making structural changes.

---

## What the system does

1. **Egor writes** hypothesis ranking + research context in Obsidian (Brain vault).
2. **A local file watcher** (`launchd` agent) regenerates `research/strc/PROBLEM.md` inside Brain. `brain-keeper` picks it up on its next poll and pushes Brain to Forgejo.
3. **A scheduled cloud agent** (Anthropic "routine" / RemoteTrigger) wakes up daily. It clones **Brain from Forgejo** (`git.egor.lol/git-egor/brain.git`), `cd`s into `research/strc/`, reads `PROBLEM.md`, scans PubMed + bioRxiv for new papers, dedups against existing `papers/*.md`, writes new ones, updates matching `hypotheses/<slug>.md` synthesis files, commits + pushes back to Forgejo.
4. **A second local file watcher** (`com.egor.strc-content-sync`) sees Brain's `research/strc/` changed and runs `strc-sync-brain-to-site`: pulls Brain from Forgejo, rsyncs papers/models/hypotheses into this site repo, pushes to GitHub.
5. **GitHub Actions** builds the Astro site on every push and deploys to Cloudflare Pages at `strc.egor.lol`.
6. **Egor reads** the result on the site and (transparently) in Obsidian, since the papers now live physically in Brain.

The whole system is an autonomous loop. Brain is the single source of truth. Site is a derived artifact.

---

## Data flow

```
                     ┌─ Egor edits in Obsidian
                     ▼
       ~/Brain/notes/STRC Hypothesis Ranking.md
       ~/Brain/research/strc/_context.md
       ~/Brain/research/strc/site-hypothesis-slugs.json
                     │
                     │ (launchd com.egor.strc-sync WatchPaths)
                     ▼
              ~/bin/strc-sync-problem
               assembles PROBLEM.md → Brain/research/strc/PROBLEM.md
                     │
                     │ (brain-keeper auto-commits + pushes)
                     ▼
              git.egor.lol/git-egor/brain ◄─────┐
                     │                          │
              ┌──────┴──────┐                   │ push papers/ + hypotheses/
              │             │                   │
              │             │             Cloud routine (daily 00:47 HKT)
              │             │             clones Brain, cd research/strc,
              │             │             scans PubMed + bioRxiv,
              │             │             writes papers/ + hypotheses/
              │             │                   ▲
              │             ▼                   │
              │     brain-keeper pulls ─────────┘
              │     (inside strc-sync-brain-to-site)
              ▼
 (launchd com.egor.strc-content-sync WatchPaths)
              │
              ▼
  ~/bin/strc-sync-brain-to-site
   pull Brain ← Forgejo
   rsync Brain/research/strc/{papers,models,hypotheses}/ → site/
   cp   Brain/notes/STRC Hypothesis Ranking.md         → site/src/data/hypothesis-ranking.md
   cp   Brain/research/strc/site-hypothesis-slugs.json → site/src/data/hypothesis-slugs.json
   git commit + push to github.com/lyfar/strc-research
              │
              ▼
     GitHub Actions → Cloudflare Pages build
              │
              ▼
        strc.egor.lol (site)
```

---

## Canonicity

All STRC research content is canonical in Brain. The site is a projection. Anything you edit directly in the site repo under `papers/`, `models/`, `hypotheses/`, `src/data/hypothesis-ranking.md`, or `src/data/hypothesis-slugs.json` will be **overwritten on next sync**.

| Path | Canonical location | Writer |
|---|---|---|
| `STRC Hypothesis Ranking.md` | `~/Brain/notes/` | Egor (Obsidian) |
| `_context.md` | `~/Brain/research/strc/` | Egor (Obsidian) |
| `site-hypothesis-slugs.json` | `~/Brain/research/strc/` | Egor (rarely) |
| `PROBLEM.md` | `~/Brain/research/strc/` (regen) | `strc-sync-problem` script |
| `papers/*.md` | `~/Brain/research/strc/papers/` | cloud routine |
| `models/*` | `~/Brain/research/strc/models/` | Egor's compute scripts |
| `hypotheses/*.md` | `~/Brain/research/strc/hypotheses/` | cloud routine |
| `src/pages/**`, `astro.config.mjs`, etc. | site repo | Egor |

Pre-migration (before 2026-04-22), papers and models lived in the site repo and Brain symlinked into it. That was required because the cloud routine ran against GitHub, which can't resolve local-path symlinks on the runner. Flipping the routine to clone Brain directly from Forgejo removed that constraint.

### iCloud + git

- Brain lives in iCloud (`~/Brain` → `~/Library/Mobile Documents/iCloud~md~obsidian/Documents/Brain`) and is simultaneously a git repo pushed to Forgejo at `git.egor.lol/git-egor/brain`.
- Site repo is NOT in iCloud; lives at `~/Sites/site-strc-egor-lol/`, tracked by git, pushed to GitHub.
- Brain's `.gitignore` excludes `.obsidian/`, `.trash/`, `.DS_Store`, `*.tmp`.
- Site's `.gitignore` excludes `models/**/__pycache__`, `*.png`, `*.jpg`, `.DS_Store` + build artifacts.

---

## Local scripts

### `~/bin/strc-sync-problem`

Regenerates `Brain/research/strc/PROBLEM.md` from three Brain sources. Triggered by `com.egor.strc-sync.plist` WatchPaths on the three source files.

```
~/Brain/notes/STRC Hypothesis Ranking.md  ─┐
~/Brain/research/strc/_context.md          ├─► assemble PROBLEM.md ─► Brain/research/strc/PROBLEM.md
~/Brain/research/strc/site-hypothesis-slugs.json ─┘
```

Doesn't commit — `brain-keeper` handles pushing Brain to Forgejo.

### `~/bin/strc-sync-brain-to-site`

Mirrors Brain's STRC research content into the site repo. Triggered by `com.egor.strc-content-sync.plist` WatchPaths on:
- `Brain/research/strc/papers/`
- `Brain/research/strc/models/`
- `Brain/research/strc/hypotheses/`
- `Brain/notes/STRC Hypothesis Ranking.md`
- `Brain/research/strc/site-hypothesis-slugs.json`

Plus `StartInterval=900` (15 min fallback).

Steps:
1. `git pull --rebase --autostash forgejo main` in Brain (pick up cloud-routine writes)
2. `rsync --delete` the three content dirs Brain → site
3. `cp` the ranking + slugs files
4. `git pull --rebase --autostash origin master` in site
5. Re-rsync after pull (Brain state wins)
6. `git add` + commit + push site to GitHub

Idempotent: if no diff after rsync, exits early without committing.

### `brain-keeper` (external, pre-existing)

Compiled binary at `/usr/local/bin/brain-keeper`, launched by `com.egor.brain-keeper.plist` with `KeepAlive=true`. Polls Brain for local changes, commits with `auto: YYYY-MM-DD (N files changed)` messages, pushes to Forgejo. We don't have its source; black-box.

---

## PROBLEM.md — the agent's context file

**Location:** `research/strc/PROBLEM.md` inside Brain (was: site repo root pre-migration).

**Purpose:** single compiled context file the cloud routine reads on every run. Without it, every routine run would need its context hardcoded in the prompt, and the prompt would drift from Egor's priorities.

Sections:
1. **Core problem** — STRC biology, Misha's genotype, engineering constraints. From `_context.md`.
2. **Active hypothesis register** — 25-row markdown table. From `STRC Hypothesis Ranking.md`.
3. **Hypothesis → synthesis file slug table** — from `site-hypothesis-slugs.json`.
4. **Kill list** — D-tier. From ranking note.
5. **Agent directive** — how to use the file.

**Regeneration trigger:** `com.egor.strc-sync.plist` watches the three source files in Brain, sleeps 10s after any change (lets Obsidian finish atomic-write), runs `strc-sync-problem`.

**Maintenance:** zero. Egor never edits PROBLEM.md directly. To change agent orientation, edit the upstream Brain files.

---

## The cloud routine

**ID:** `trig_01XfQV8QaaR5fjLPczpJJXVu`
**Name:** "STRC Research Monitor"
**Manage:** https://claude.ai/code/scheduled/trig_01XfQV8QaaR5fjLPczpJJXVu
**Schedule:** `47 16 * * *` UTC → 00:47 Asia/Hong_Kong daily
**Model:** claude-opus-4-6
**Source:** `https://git.egor.lol/git-egor/brain.git` (Forgejo, private; auth via embedded PAT in URL)
**Working dir:** `research/strc/` inside cloned Brain

**Prompt flow (10 steps):** unchanged from pre-migration except:
- Step 0 clones Brain (not GitHub), `cd research/strc`
- Step 7/10 push to Forgejo (origin main), not GitHub

Full prompt lives in the routine config, readable via `RemoteTrigger get`.

**Collision handling:** routine starts with `git pull --rebase --autostash origin main`. `brain-keeper` pushing local Brain between routine steps is harmless because the push is the last step. If routine's final push gets rejected (brain-keeper pushed concurrently), retry-on-failure logic would help — currently not present; will surface if it happens.

---

## CF Pages deployment

**Workflow:** `.github/workflows/deploy.yml`. Triggered on push to `master`. Runs `npm ci && npm run build && wrangler pages deploy dist --project-name=strc-egor-lol --branch=main`.

**Why `--branch=main` when the git branch is `master`:** CF Pages was set up with a `main` branch, we kept the git branch as `master`. This mismatch is fine but confusing. Don't "fix" it without coordinating.

**Build time:** ~40-60 seconds from push to live on strc.egor.lol.

**Schema:** `src/content.config.ts` defines the `papers` content collection schema. If you add a frontmatter field to papers written by the routine, also add it here (optional) or build will fail.

---

## The site pages

- **`/hypotheses`** (`src/pages/hypotheses/index.astro`) — live priority register. Parses `src/data/hypothesis-ranking.md` (copy of Brain's ranking note, maintained by `strc-sync-brain-to-site`).
- **`/hypotheses/<slug>`** (individual files under `src/pages/hypotheses/`) — hand-authored long-form articles per hypothesis.
- **`/journal`** (`src/pages/journal.astro`) — research monitor. Sorts `papers/*.md` by `date_added` then by publication date.

---

## File inventory — what each file does

### In site repo (downstream from Brain)
- `papers/*.md` — mirror of `~/Brain/research/strc/papers/`. Written by cloud routine (to Brain), synced here.
- `models/*` — mirror of `~/Brain/research/strc/models/`. Egor's compute outputs (produced in Brain), synced here.
- `hypotheses/<slug>.md` — mirror of `~/Brain/research/strc/hypotheses/`. Auto-maintained synthesis by cloud routine.
- `src/data/hypothesis-ranking.md` — copy of `~/Brain/notes/STRC Hypothesis Ranking.md`.
- `src/data/hypothesis-slugs.json` — copy of `~/Brain/research/strc/site-hypothesis-slugs.json`.
- `src/pages/**`, `src/lib/**`, `src/content.config.ts` — Astro source (authored in site, not derived).
- `.github/workflows/deploy.yml` — CF Pages deploy.

### In Brain (canonical)
- `~/Brain/notes/STRC Hypothesis Ranking.md` — ranking (Egor-authored).
- `~/Brain/research/strc/_context.md` — problem context (Egor-authored).
- `~/Brain/research/strc/site-hypothesis-slugs.json` — slug map (Egor-maintained).
- `~/Brain/research/strc/PROBLEM.md` — compiled context for routine (script-generated).
- `~/Brain/research/strc/papers/` — paper index (routine-written).
- `~/Brain/research/strc/hypotheses/` — synthesis files (routine-written).
- `~/Brain/research/strc/models/` — compute outputs (Egor's scripts).

### Outside both repos
- `~/bin/strc-sync-problem` — PROBLEM.md generator.
- `~/bin/strc-sync-brain-to-site` — Brain→site mirror.
- `/usr/local/bin/brain-keeper` — Brain auto-committer.
- `~/Library/LaunchAgents/com.egor.strc-sync.plist` — watcher for PROBLEM.md sources.
- `~/Library/LaunchAgents/com.egor.strc-content-sync.plist` — watcher for Brain→site sync.
- `~/Library/LaunchAgents/com.egor.brain-keeper.plist` — brain-keeper daemon.
- `~/Library/Logs/strc-sync.log` — PROBLEM.md sync log.
- `~/Library/Logs/strc-content-sync.log` — Brain→site sync log.

---

## Failure modes + how to debug

### Site not updating despite Brain changes
1. Check `~/Library/Logs/strc-content-sync.log` — latest rsync run.
2. Run `~/bin/strc-sync-brain-to-site` manually. Watch output.
3. If "pull failed" — site repo's working tree may have conflicting changes. `cd ~/Sites/site-strc-egor-lol && git status`.
4. Launchd status: `launchctl list | grep strc-content-sync`.

### Routine failing to write
1. Check routine run history: `RemoteTrigger get trig_01XfQV8QaaR5fjLPczpJJXVu` — inspect last run.
2. Verify PROBLEM.md is current on Forgejo: `cd ~/Brain && git log --oneline forgejo/main -5`.
3. Verify Brain is pushed: `cd ~/Brain && git status` — should show nothing relevant behind forgejo.

### brain-keeper silent
1. Check `launchctl list | grep brain-keeper`.
2. `ps aux | grep brain-keeper`.
3. Binary location: `/usr/local/bin/brain-keeper` (Mach-O, source not in repo).
4. Worst case: manual `cd ~/Brain && git add -A && git commit -m "manual" && git push forgejo main`.

### Cloud routine + brain-keeper race
If both push within seconds, one loses with non-FF error.
- brain-keeper: we don't know its retry behavior; KeepAlive=true means at least the daemon stays up.
- Routine: no built-in retry. A failed push means that day's papers sit only in routine's ephemeral workspace — lost on next run.
- Workaround: schedule routine off-peak (00:47 HKT is already during Egor's sleep, low chance of brain-keeper activity from manual edits).
- Fix if it becomes a problem: add retry loop to routine's Step 7/10.
