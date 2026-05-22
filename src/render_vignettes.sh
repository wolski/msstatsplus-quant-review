#!/usr/bin/env bash
# Render all vignettes/*.qmd to HTML, best-effort.
#
# Scope: this script ONLY renders. It does not stage, copy, commit, or push.
# Site staging and gh-pages publishing live in src/publish_gh_pages.sh.
#
# Usage:
#   make render-vignettes
#   # or directly:
#   bash src/render_vignettes.sh
#
# Behaviour:
#   - Each entry in RENDER is rendered with `quarto render <f> --to html`.
#   - A failure on one file does not abort the others.
#   - Outputs land next to the source in vignettes/<name>.html (and
#     vignettes/<name>_files/ where Quarto needs it).
#   - Exit code: 0 if all renders succeeded, 1 otherwise.

set -u
set -o pipefail

# Resolve quant/ as the directory above this script (src/ lives in quant/).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
QUANT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
VIG_DIR="$QUANT_DIR/vignettes"

RENDER=(
  index.qmd
  review.qmd
  review2.qmd
  review_main.qmd
  review_supplement.qmd
  CSF_Spectronaut_data_vis.qmd
  CSF_Spectronaut_sample_swap_data_vis.qmd
  CSF_Spectronaut_protein_swap_data_vis.qmd
  swap_pairs_before_after.qmd
)
# Intentionally omitted: diagnostics.qmd, swap_visualization.qmd
# These are parameterized matrix-render targets; their YAML defaults do
# not produce a meaningful standalone render and they fail every run.
# They are rendered elsewhere via per-cell param overrides.

if ! command -v quarto >/dev/null 2>&1; then
  echo "ERROR: quarto not found on \$PATH." >&2
  exit 2
fi

if [ ! -d "$VIG_DIR" ]; then
  echo "ERROR: vignettes dir not found: $VIG_DIR" >&2
  exit 2
fi

OK=()
FAILED=()

echo "Rendering ${#RENDER[@]} vignette(s) from $VIG_DIR"
for f in "${RENDER[@]}"; do
  if [ ! -f "$VIG_DIR/$f" ]; then
    echo "  SKIP (missing source): $f"
    FAILED+=("$f (missing)")
    continue
  fi
  echo "  rendering: $f"
  if (cd "$VIG_DIR" && quarto render "$f" --to html); then
    OK+=("$f")
  else
    FAILED+=("$f")
  fi
done

echo
echo "Render summary:"
echo "  ok:     ${#OK[@]}"
for f in "${OK[@]}";     do echo "    + $f"; done
echo "  failed: ${#FAILED[@]}"
for f in "${FAILED[@]}"; do echo "    - $f"; done

if [ "${#FAILED[@]}" -gt 0 ]; then
  exit 1
fi
