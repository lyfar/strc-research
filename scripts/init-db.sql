-- STRC Tools Database Schema
-- Short descriptions + operational data here. Full content in Brain vault.
-- Run: sqlite3 data/tools.db < scripts/init-db.sql

CREATE TABLE IF NOT EXISTS tools (
  slug TEXT PRIMARY KEY,
  brain_note TEXT NOT NULL,           -- Brain vault note filename (without .md)
  name TEXT NOT NULL,                 -- display name on site
  category TEXT NOT NULL,
  description TEXT NOT NULL DEFAULT '',-- short 1-2 sentence description
  url TEXT NOT NULL DEFAULT '',       -- primary tool URL
  status TEXT NOT NULL DEFAULT 'stub'
    CHECK (status IN ('verified', 'tested', 'available', 'stub', 'site-down', 'link-broken')),
  is_free INTEGER NOT NULL DEFAULT 1,
  has_api INTEGER NOT NULL DEFAULT 0,
  acmg_relevance TEXT,
  sort_order INTEGER DEFAULT 0,
  created_at TEXT NOT NULL DEFAULT (datetime('now')),
  updated_at TEXT NOT NULL DEFAULT (datetime('now'))
);

CREATE TABLE IF NOT EXISTS tags (
  tool_slug TEXT NOT NULL REFERENCES tools(slug) ON DELETE CASCADE,
  tag TEXT NOT NULL,
  PRIMARY KEY (tool_slug, tag)
);

CREATE TABLE IF NOT EXISTS e1659a_results (
  tool_slug TEXT PRIMARY KEY REFERENCES tools(slug) ON DELETE CASCADE,
  score TEXT NOT NULL,
  interpretation TEXT NOT NULL,
  raw_data TEXT,
  tested_at TEXT NOT NULL DEFAULT (datetime('now')),
  notes TEXT
);

CREATE TABLE IF NOT EXISTS url_checks (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  tool_slug TEXT NOT NULL REFERENCES tools(slug) ON DELETE CASCADE,
  http_code INTEGER,
  response_ms INTEGER,
  checked_at TEXT NOT NULL DEFAULT (datetime('now')),
  error TEXT
);

CREATE TABLE IF NOT EXISTS categories (
  name TEXT PRIMARY KEY,
  sort_order INTEGER NOT NULL
);

INSERT OR IGNORE INTO categories (name, sort_order) VALUES
  ('AI & Foundation Models', 1),
  ('Structural Biology', 2),
  ('Variant Effect Prediction', 3),
  ('Splicing Prediction', 4),
  ('Regulatory & Non-Coding', 5),
  ('Population Databases', 6),
  ('Clinical Databases', 7),
  ('Gene-Level Resources', 8),
  ('Hearing Loss & Inner Ear', 9),
  ('Conservation & Evolution', 10),
  ('Structural Variants & CNV', 11),
  ('Nomenclature & Validation', 12),
  ('Literature Mining', 13),
  ('Workflow Platforms', 14);

CREATE TRIGGER IF NOT EXISTS tools_updated_at
AFTER UPDATE ON tools
FOR EACH ROW
BEGIN
  UPDATE tools SET updated_at = datetime('now') WHERE slug = NEW.slug;
END;

CREATE INDEX IF NOT EXISTS idx_tools_status ON tools(status);
CREATE INDEX IF NOT EXISTS idx_tools_category ON tools(category);
CREATE INDEX IF NOT EXISTS idx_url_checks_slug ON url_checks(tool_slug);
