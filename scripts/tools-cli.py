#!/usr/bin/env python3
"""CLI for managing STRC tools database.

Usage:
  python3 scripts/tools-cli.py status              — show summary
  python3 scripts/tools-cli.py check-urls           — check all URLs and record results
  python3 scripts/tools-cli.py list [--status=X]    — list tools, optionally filtered
  python3 scripts/tools-cli.py get <slug>           — show tool details
  python3 scripts/tools-cli.py set <slug> <field> <value> — update a field
  python3 scripts/tools-cli.py add-e1659a <slug> <score> <interpretation> — add E1659A result
  python3 scripts/tools-cli.py health               — run health checks
  python3 scripts/tools-cli.py generate             — regenerate tools.ts from db
"""
import sys, sqlite3, subprocess, time
from pathlib import Path
from urllib.request import urlopen, Request
from urllib.error import URLError, HTTPError

ROOT = Path(__file__).parent.parent
DB_FILE = ROOT / "data/tools.db"


def get_db():
    db = sqlite3.connect(str(DB_FILE))
    db.row_factory = sqlite3.Row
    db.execute("PRAGMA foreign_keys = ON")
    return db


def cmd_status():
    db = get_db()
    total = db.execute("SELECT COUNT(*) FROM tools").fetchone()[0]
    by_status = db.execute("SELECT status, COUNT(*) as cnt FROM tools GROUP BY status ORDER BY cnt DESC").fetchall()
    by_cat = db.execute("""
        SELECT t.category, COUNT(*) as cnt
        FROM tools t JOIN categories c ON t.category = c.name
        GROUP BY t.category ORDER BY c.sort_order
    """).fetchall()
    e1659a = db.execute("SELECT COUNT(*) FROM e1659a_results").fetchone()[0]
    tags = db.execute("SELECT COUNT(DISTINCT tag) FROM tags").fetchone()[0]
    print(f"STRC Tools Database")
    print(f"{'='*40}")
    print(f"Total tools:    {total}")
    print(f"E1659A tested:  {e1659a}")
    print(f"Unique tags:    {tags}")
    print(f"\nBy status:")
    icons = {'verified': '✅', 'tested': '🧪', 'available': '🔵', 'stub': '⬜', 'site-down': '🔴', 'link-broken': '💀'}
    for r in by_status:
        print(f"  {icons.get(r['status'], '?')} {r['status']:12s} {r['cnt']}")
    print(f"\nBy category:")
    for r in by_cat:
        print(f"  {r['category']:30s} {r['cnt']}")
    db.close()


def cmd_list(status_filter=None):
    db = get_db()
    if status_filter:
        rows = db.execute("SELECT slug, name, status, category FROM tools WHERE status = ? ORDER BY category, name", (status_filter,)).fetchall()
    else:
        rows = db.execute("""
            SELECT t.slug, t.name, t.status, t.category
            FROM tools t JOIN categories c ON t.category = c.name
            ORDER BY c.sort_order, t.name
        """).fetchall()

    icons = {'verified': '✅', 'tested': '🧪', 'available': '🔵', 'stub': '⬜', 'site-down': '🔴', 'link-broken': '💀'}
    current_cat = None
    for r in rows:
        if r['category'] != current_cat:
            current_cat = r['category']
            print(f"\n  {current_cat}")
        print(f"    {icons.get(r['status'], '?')} {r['name']:40s} [{r['slug']}]")
    print(f"\n  Total: {len(rows)}")
    db.close()


def cmd_get(slug):
    db = get_db()
    t = db.execute("SELECT * FROM tools WHERE slug = ?", (slug,)).fetchone()
    if not t:
        print(f"Tool '{slug}' not found")
        return
    e1659a = db.execute("SELECT * FROM e1659a_results WHERE tool_slug = ?", (slug,)).fetchone()
    last_check = db.execute("SELECT * FROM url_checks WHERE tool_slug = ? ORDER BY checked_at DESC LIMIT 1", (slug,)).fetchone()

    tags = [r['tag'] for r in db.execute("SELECT tag FROM tags WHERE tool_slug = ?", (slug,))]
    print(f"{'='*50}")
    print(f"  Name:        {t['name']}")
    print(f"  Slug:        {t['slug']}")
    print(f"  Brain note:  {t['brain_note']}")
    print(f"  Category:    {t['category']}")
    print(f"  Status:      {t['status']}")
    print(f"  URL:         {t['url']}")
    print(f"  Description: {t['description'][:100]}")
    print(f"  Tags:        {', '.join(tags) or 'none'}")
    print(f"  Free:        {'Yes' if t['is_free'] else 'No'}")
    print(f"  API:         {'Yes' if t['has_api'] else 'No'}")
    print(f"  ACMG:        {t['acmg_relevance'] or 'N/A'}")
    if e1659a:
        print(f"  E1659A:      {e1659a['score']} — {e1659a['interpretation']}")
    if last_check:
        print(f"  Last check:  HTTP {last_check['http_code']} ({last_check['response_ms']}ms) at {last_check['checked_at']}")
    print(f"  Updated:     {t['updated_at']}")
    db.close()


def cmd_set(slug, field, value):
    db = get_db()
    valid_fields = ['name', 'brain_note', 'category', 'url', 'status', 'is_free', 'has_api', 'acmg_relevance']
    if field not in valid_fields:
        print(f"Invalid field. Valid: {', '.join(valid_fields)}")
        return
    if field in ('is_free', 'has_api'):
        value = 1 if value.lower() in ('true', '1', 'yes') else 0
    db.execute(f"UPDATE tools SET {field} = ? WHERE slug = ?", (value, slug))
    db.commit()
    print(f"Updated {slug}.{field} = {value}")
    db.close()


def cmd_add_e1659a(slug, score, interpretation):
    db = get_db()
    db.execute("""
        INSERT OR REPLACE INTO e1659a_results (tool_slug, score, interpretation)
        VALUES (?, ?, ?)
    """, (slug, score, interpretation))
    db.commit()
    print(f"E1659A result for {slug}: {score} — {interpretation}")
    db.close()


def check_url(url, timeout=10):
    """Check a URL, return (http_code, response_ms, error)."""
    if not url or url.startswith(('repo does not exist', 'GitHub', 'SITE', 'BROKEN', '404', '500', 'REDIRECTS', 'SHINY')):
        return 0, 0, 'invalid_url'
    try:
        start = time.time()
        req = Request(url, headers={'User-Agent': 'Mozilla/5.0 STRC-Audit/1.0'})
        resp = urlopen(req, timeout=timeout)
        ms = int((time.time() - start) * 1000)
        return resp.getcode(), ms, None
    except HTTPError as e:
        ms = int((time.time() - start) * 1000)
        return e.code, ms, None
    except URLError as e:
        return 0, 0, str(e.reason)[:200]
    except Exception as e:
        return 0, 0, str(e)[:200]


def cmd_check_urls():
    db = get_db()
    tools = db.execute("SELECT slug, name, url, status FROM tools ORDER BY slug").fetchall()
    print(f"Checking {len(tools)} URLs...")

    results = {'ok': 0, 'warn': 0, 'fail': 0}
    status_updates = []

    for t in tools:
        code, ms, error = check_url(t['url'])

        # Record check
        db.execute("""
            INSERT INTO url_checks (tool_slug, http_code, response_ms, error)
            VALUES (?, ?, ?, ?)
        """, (t['slug'], code, ms, error))

        # Determine result
        if code == 200:
            icon = '✅'
            results['ok'] += 1
            # If was site-down, suggest status change
            if t['status'] == 'site-down':
                status_updates.append((t['slug'], t['name'], 'site-down', 'available'))
        elif code in (301, 302, 303, 403):
            icon = '⚠️'
            results['warn'] += 1
        else:
            icon = '❌'
            results['fail'] += 1
            # If was available/tested, suggest status change
            if t['status'] in ('available', 'tested', 'verified') and code == 0 and error != 'invalid_url':
                status_updates.append((t['slug'], t['name'], t['status'], 'site-down'))

        if code != 200:
            print(f"  {icon} {t['name']:35s} HTTP {code:3d} {ms:5d}ms {error or ''}")

    db.commit()

    print(f"\nResults: {results['ok']} ok, {results['warn']} warnings, {results['fail']} failures")

    if status_updates:
        print(f"\nSuggested status changes:")
        for slug, name, old, new in status_updates:
            print(f"  {name}: {old} -> {new}")
            print(f"    python3 scripts/tools-cli.py set {slug} status {new}")

    db.close()


def cmd_health():
    db = get_db()
    print("Health Check")
    print("=" * 40)

    # Tools with empty descriptions
    empty_desc = db.execute("SELECT slug, name FROM tools WHERE description = '' OR description IS NULL").fetchall()
    print(f"\nEmpty descriptions: {len(empty_desc)}")
    for r in empty_desc:
        print(f"  - {r['name']} [{r['slug']}]")

    # Tools with no URL
    no_url = db.execute("SELECT slug, name FROM tools WHERE url = '' OR url IS NULL").fetchall()
    print(f"\nMissing URLs: {len(no_url)}")
    for r in no_url:
        print(f"  - {r['name']} [{r['slug']}]")

    # Tools with invalid URLs (contain spaces or no http)
    bad_url = db.execute("SELECT slug, name, url FROM tools WHERE url != '' AND url NOT LIKE 'http%'").fetchall()
    print(f"\nInvalid URLs: {len(bad_url)}")
    for r in bad_url:
        print(f"  - {r['name']}: {r['url']}")

    # Stale verified tools (no URL check in 7+ days)
    stale = db.execute("""
        SELECT t.slug, t.name, MAX(u.checked_at) as last_check
        FROM tools t LEFT JOIN url_checks u ON t.slug = u.tool_slug
        WHERE t.status IN ('verified', 'tested')
        GROUP BY t.slug
        HAVING last_check IS NULL OR last_check < datetime('now', '-7 days')
    """).fetchall()
    print(f"\nVerified/tested tools without recent URL check: {len(stale)}")

    # Tools without Brain notes
    brain = Path.home() / "Library/Mobile Documents/iCloud~md~obsidian/Documents/Brain/notes"
    no_brain = []
    for r in db.execute("SELECT slug, brain_note FROM tools"):
        if not (brain / f"{r['brain_note']}.md").exists():
            no_brain.append(r['brain_note'])
    print(f"\nMissing Brain notes: {len(no_brain)}")
    for name in no_brain[:15]:
        print(f"  - {name}")
    if len(no_brain) > 15:
        print(f"  ... and {len(no_brain) - 15} more")

    db.close()


def cmd_generate():
    """Run db-to-ts.py to regenerate tools.ts."""
    script = ROOT / "scripts/db-to-ts.py"
    subprocess.run([sys.executable, str(script)], check=True)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return

    cmd = sys.argv[1]

    if cmd == 'status':
        cmd_status()
    elif cmd == 'check-urls':
        cmd_check_urls()
    elif cmd == 'list':
        status = None
        for arg in sys.argv[2:]:
            if arg.startswith('--status='):
                status = arg.split('=')[1]
        cmd_list(status)
    elif cmd == 'get' and len(sys.argv) >= 3:
        cmd_get(sys.argv[2])
    elif cmd == 'set' and len(sys.argv) >= 5:
        cmd_set(sys.argv[2], sys.argv[3], sys.argv[4])
    elif cmd == 'add-e1659a' and len(sys.argv) >= 5:
        cmd_add_e1659a(sys.argv[2], sys.argv[3], ' '.join(sys.argv[4:]))
    elif cmd == 'health':
        cmd_health()
    elif cmd == 'generate':
        cmd_generate()
    else:
        print(__doc__)


if __name__ == '__main__':
    main()
