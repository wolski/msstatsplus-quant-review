#!/usr/bin/env bash
# Stage rendered HTML and force-push it as an orphan gh-pages branch.
#
# Scope: this script ONLY stages and publishes. It does not render.
# Rendering lives in src/render_vignettes.sh; chain them explicitly:
#   make render-vignettes && make gh-pages
#
# Behaviour:
#   - Aborts if any PUBLISH slug is missing vignettes/<slug>.html.
#   - Copies <slug>.html and (if present) <slug>_files/ into a temp dir.
#   - Adds .nojekyll so directories starting with _ are served.
#   - Creates an orphan worktree on gh-pages, commits, force-pushes to origin.
#   - The original working tree on main is never modified.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
QUANT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
VIG_DIR="$QUANT_DIR/vignettes"

PUBLISH=(
  index
  review_supplement
  CSF_Spectronaut_data_vis
  CSF_Spectronaut_sample_swap_data_vis
  CSF_Spectronaut_protein_swap_data_vis
  swap_pairs_before_after
)
# Intentionally omitted: diagnostics, swap_visualization
# These qmds are part of the parameterized matrix-render pipeline; their
# YAML defaults do not produce a meaningful standalone render.

BRANCH=gh-pages
REMOTE=origin

# --- sanity checks ---------------------------------------------------------

if ! command -v git >/dev/null 2>&1; then
  echo "ERROR: git not found on \$PATH." >&2
  exit 2
fi

cd "$QUANT_DIR"

if ! git rev-parse --git-dir >/dev/null 2>&1; then
  echo "ERROR: $QUANT_DIR is not a git repository." >&2
  exit 2
fi

if ! git remote get-url "$REMOTE" >/dev/null 2>&1; then
  echo "ERROR: git remote '$REMOTE' is not configured." >&2
  exit 2
fi

MAIN_SHA="$(git rev-parse --short HEAD)"
MAIN_REF="$(git rev-parse --abbrev-ref HEAD)"

MISSING=()
for slug in "${PUBLISH[@]}"; do
  if [ ! -f "$VIG_DIR/$slug.html" ]; then
    MISSING+=("$slug.html")
  fi
done
if [ "${#MISSING[@]}" -gt 0 ]; then
  echo "ERROR: required HTML files are missing in $VIG_DIR :" >&2
  for m in "${MISSING[@]}"; do echo "  - $m" >&2; done
  echo "Run 'make render-vignettes' first." >&2
  exit 1
fi

# --- stage site dir --------------------------------------------------------

TMP="$(mktemp -d)"
SITE="$TMP/site"
mkdir -p "$SITE"

echo "Staging ${#PUBLISH[@]} page(s) into $SITE"
for slug in "${PUBLISH[@]}"; do
  cp "$VIG_DIR/$slug.html" "$SITE/$slug.html"
  if [ -d "$VIG_DIR/${slug}_files" ]; then
    cp -R "$VIG_DIR/${slug}_files" "$SITE/"
  fi
done

# .nojekyll so directories starting with _ are served by GH Pages.
touch "$SITE/.nojekyll"

# --- publish via orphan worktree ------------------------------------------

WT="$TMP/wt"

cleanup() {
  # Remove the worktree (force, in case the commit step failed mid-way).
  if [ -d "$WT" ]; then
    git worktree remove --force "$WT" >/dev/null 2>&1 || rm -rf "$WT"
  fi
  rm -rf "$TMP"
}
trap cleanup EXIT

# Make the worktree creation idempotent across re-runs:
#   1. Prune any stale worktree entries from earlier aborted runs.
#   2. Delete the local gh-pages branch if it exists. We always force-push,
#      so the local ref is disposable; this avoids `worktree add --orphan`
#      rejecting an already-existing branch (observed on git 2.50.x).
git worktree prune >/dev/null 2>&1 || true
if git show-ref --verify --quiet "refs/heads/$BRANCH"; then
  git branch -D "$BRANCH" >/dev/null
fi

echo "Creating orphan worktree on $BRANCH at $WT"
git worktree add --orphan -b "$BRANCH" "$WT" >/dev/null

# Clear anything the worktree inherited (orphan starts from current index).
( cd "$WT" && git rm -rf --quiet . >/dev/null 2>&1 || true )
find "$WT" -mindepth 1 -maxdepth 1 ! -name '.git' -exec rm -rf {} +

# Copy the staged site into the worktree.
cp -R "$SITE/." "$WT/"

cd "$WT"

# Identity falls back to the user's global git config when unset.
git add -A
COMMIT_MSG="Publish vignettes $(date -u +%FT%TZ) — $MAIN_REF@$MAIN_SHA"
git commit -m "$COMMIT_MSG"

PUBLISH_SHA="$(git rev-parse --short HEAD)"

echo "Force-pushing $BRANCH to $REMOTE"
git push -f "$REMOTE" "$BRANCH"

cd "$QUANT_DIR"

# --- final report ----------------------------------------------------------

REMOTE_URL="$(git remote get-url "$REMOTE")"
echo
echo "Published:"
echo "  branch:     $BRANCH"
echo "  remote:     $REMOTE_URL"
echo "  commit:     $PUBLISH_SHA"
echo "  from main:  $MAIN_REF@$MAIN_SHA"
echo
echo "If this was the first push, enable GitHub Pages once:"
echo "  Settings -> Pages -> Source: Deploy from a branch -> Branch: $BRANCH / (root)"
