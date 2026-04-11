#!/bin/bash
# Cleanup script: remove old static tool pages after migration to dynamic [slug].astro
# Run AFTER verifying the dynamic route works: npm run build && npm run preview
#
# This removes:
# 1. All 143+ individual .astro tool pages (replaced by [slug].astro)
# 2. Category sub-pages (replaced by #anchors on main /tools page)
# 3. Old ToolsGrid component (replaced by ToolCard)

set -e

SITE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
TOOL_DIR="$SITE_DIR/src/pages/tools/tool"
TOOLS_DIR="$SITE_DIR/src/pages/tools"

echo "=== STRC Site: Cleanup Old Static Tool Pages ==="
echo ""

# Count files to remove
static_count=$(find "$TOOL_DIR" -name "*.astro" ! -name "[slug].astro" 2>/dev/null | wc -l | tr -d ' ')
category_count=$(find "$TOOLS_DIR" -maxdepth 1 -name "*.astro" ! -name "index.astro" 2>/dev/null | wc -l | tr -d ' ')

echo "Static tool pages to remove: $static_count"
echo "Category sub-pages to remove: $category_count"
echo ""

if [ "$1" != "--yes" ]; then
  echo "Dry run. Pass --yes to actually delete."
  echo ""
  echo "Static tool pages:"
  find "$TOOL_DIR" -name "*.astro" ! -name "[slug].astro" | sort | head -20
  echo "... (showing first 20)"
  echo ""
  echo "Category pages:"
  find "$TOOLS_DIR" -maxdepth 1 -name "*.astro" ! -name "index.astro" | sort
  exit 0
fi

echo "Removing static tool pages..."
find "$TOOL_DIR" -name "*.astro" ! -name "[slug].astro" -delete
echo "  Removed $static_count files"

echo "Removing category sub-pages..."
for f in "$TOOLS_DIR"/conservation.astro "$TOOLS_DIR"/hearing-loss.astro "$TOOLS_DIR"/literature.astro "$TOOLS_DIR"/population-genetics.astro "$TOOLS_DIR"/protein-analysis.astro "$TOOLS_DIR"/structural-biology.astro "$TOOLS_DIR"/gene-therapy.astro "$TOOLS_DIR"/variant-scoring.astro; do
  [ -f "$f" ] && rm "$f" && echo "  Removed $(basename $f)"
done

echo "Removing old ToolsGrid component..."
[ -f "$SITE_DIR/src/components/ToolsGrid.astro" ] && rm "$SITE_DIR/src/components/ToolsGrid.astro" && echo "  Removed ToolsGrid.astro"

echo ""
echo "Done! Now run: npm run build && npm run preview"
