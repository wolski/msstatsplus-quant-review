# GitHub Pages publish for `quant/vignettes/`

Plan for a single bash script that renders all `vignettes/*.qmd` and force-pushes
a selected subset to an orphan `gh-pages` branch on `origin`.

## Context (discovered, not assumed)

- `RMSV000000701.3-rerun/quant/` is its own git repo.
  - Remote: `https://github.com/wolski/msstatsplus-quant-review.git`
  - Current branch: `main`, up-to-date with `origin/main`.
  - No `gh-pages` branch yet (`git branch -a` shows only `main` and `remotes/origin/main`).
- Vignette sources in `vignettes/`:
  `review.qmd`, `review2.qmd`, `review_main.qmd`, `review_supplement.qmd`,
  `diagnostics.qmd`, `CSF_Spectronaut_data_vis.qmd`,
  `CSF_Spectronaut_sample_swap_data_vis.qmd`,
  `CSF_Spectronaut_protein_swap_data_vis.qmd`, `swap_visualization.qmd`.
- Existing Makefile target `review-html` runs `cd vignettes && quarto render review.qmd --to html`. We will reuse Quarto via the script but not depend on the heavy cell stamps.

## Design principle: render and staging are separate

Two scripts, two Makefile targets, no implicit dependency between them.
Running the site is an explicit two-command workflow:

```
make render-vignettes   # qmd -> html, including index.qmd. No git, no push.
make gh-pages           # stage selected HTML into a temp orphan worktree, force-push.
```

The staging step assumes the HTML it needs already exists in `vignettes/`. If
something is missing it logs a clear error and refuses to publish that file;
it does NOT silently re-render. This means you can iterate on rendering
without ever publishing, and you can re-publish (e.g. after editing
`index.qmd` or the `PUBLISH` list) without re-rendering anything heavy.

## Files to add (inside `quant/`)

1. **`vignettes/index.qmd`** — Quarto-themed landing page (Option C).
   Standalone qmd (no `_quarto.yml`), renders with the same per-file flow
   the other vignettes use, so it cannot interfere with the existing
   per-file render commands. Hand-curated link list to the `PUBLISH`
   set. The list of links in `index.qmd` and the `PUBLISH` array in the
   publish script are the two places to update when adding/removing a
   page — kept deliberately in sync by the human, not auto-generated.

2. **`src/render_vignettes.sh`** — rendering only. No git interaction.
   Configurable array at the top:
   ```bash
   RENDER=(index.qmd
           review.qmd review2.qmd review_main.qmd review_supplement.qmd
           diagnostics.qmd CSF_Spectronaut_data_vis.qmd
           CSF_Spectronaut_sample_swap_data_vis.qmd
           CSF_Spectronaut_protein_swap_data_vis.qmd
           swap_visualization.qmd)
   ```
   - Each entry is rendered best-effort: one qmd failing does not abort
     the others; failures are summarised at the end with a non-zero
     exit code.
   - Output lands next to the source in `vignettes/<name>.html` (and
     `vignettes/<name>_files/` where Quarto needs it).
   - The script does nothing else — no copy, no temp dirs, no git.

3. **`src/publish_gh_pages.sh`** — staging + push only. Does NOT render.
   Configurable array at the top:
   ```bash
   PUBLISH=(index review_supplement diagnostics
            CSF_Spectronaut_data_vis
            CSF_Spectronaut_sample_swap_data_vis
            CSF_Spectronaut_protein_swap_data_vis
            swap_visualization)
   ```
   - For each slug, expects `vignettes/<slug>.html` to already exist;
     if any required HTML is missing, the script aborts before touching
     git and prints the missing files.
   - `review.qmd`, `review2.qmd`, `review_main.qmd` are intentionally
     absent from `PUBLISH` — they may be rendered (so failures surface)
     but they are not on the site.
   - Copies each `<slug>.html` and (if present) `<slug>_files/` into a
     temp build dir, adds `.nojekyll`, then publishes the temp dir as
     an orphan `gh-pages` branch (see step detail below).

4. **`Makefile`** — add two independent targets:
   ```make
   .PHONY: render-vignettes gh-pages
   render-vignettes:
   	bash src/render_vignettes.sh

   gh-pages:
   	bash src/publish_gh_pages.sh
   ```
   `gh-pages` deliberately does NOT depend on `render-vignettes`. The
   user chains them explicitly when they want a fresh render before
   publishing. The publish script's own missing-file check is the safety
   net against accidentally publishing a stale build.

5. **`.gitignore`** — already ignores some generated HTML in
   `vignettes/`. Add the remaining ones we don't want in `main` history
   (the HTML lives on `gh-pages` only). Minimal change.

## `src/render_vignettes.sh` — steps

1. Sanity check: `quarto` on `$PATH`.
2. Render loop:
   ```bash
   for f in "${RENDER[@]}"; do
     (cd vignettes && quarto render "$f" --to html) || FAILED+=("$f")
   done
   ```
3. Print summary: rendered OK / failed list. Exit non-zero if any failed.
   No copy, no git, no temp dir. The script's only side effect is files
   in `vignettes/`.

## `src/publish_gh_pages.sh` — steps

1. Sanity checks:
   - `git` on `$PATH`, inside the `quant/` repo, `origin` configured.
   - Working tree on `main` clean enough that the orphan worktree
     can be created (warn but allow if `main` is dirty — the script
     only reads from `vignettes/`, it doesn't `git add` on `main`).
   - For each `slug` in `PUBLISH`: `vignettes/<slug>.html` must
     exist. Abort with a clear error if any are missing.
2. Build the site dir in `$(mktemp -d)/site/`:
   - For each slug: copy `vignettes/<slug>.html` and, if present,
     `vignettes/<slug>_files/` (Quarto's external-resource dir;
     required when HTML isn't fully self-contained).
   - Write `.nojekyll` (tells GH Pages not to run Jekyll, so dirs
     starting with `_` are served).
3. Publish via a throwaway worktree:
   ```bash
   git worktree add --orphan -B gh-pages "$TMP/wt"
   rm -rf "$TMP/wt"/*
   cp -R "$TMP/site/." "$TMP/wt/"
   cd "$TMP/wt"
   git add -A
   git -c user.name="…" -c user.email="…" \
       commit -m "Publish vignettes ($(date -u +%FT%TZ)) — main@<sha>"
   git push -f origin gh-pages
   cd -
   git worktree remove --force "$TMP/wt"
   ```
   - `--orphan` makes `gh-pages` history-independent from `main`. The
     branch only ever holds the latest build.
   - `git push -f` is normal here (no useful history is lost — the
     previous build is just an older snapshot of the same site). We
     never force-push `main`.
   - The original working tree on `main` is untouched.
4. Final report: pushed commit SHA / site URL.

## After the first push (one-time GitHub setup)

GitHub → repo Settings → Pages → "Build and deployment":
- Source: *Deploy from a branch*
- Branch: `gh-pages` / `/ (root)` → Save.

Site URL: `https://wolski.github.io/msstatsplus-quant-review/`.

## Index page — decision: Option C (Quarto-rendered `index.qmd`)

`vignettes/index.qmd` is a standalone qmd with its own `format: html:`
YAML block, no `_quarto.yml`. It renders with the same per-file Quarto
invocation the other vignettes use, so Quarto never enters project mode
and the existing `make review-html` etc. continue to write `review.html`
next to `review.qmd` as today.

Content: title, one short paragraph explaining the site, and a bullet
list of links to the published vignettes (each `<slug>.html`). Titles
inside the bullets are copied from the target qmds' YAML `title:` once
when authoring `index.qmd` and updated by hand when pages are
added/removed.

Risk we considered and rejected: a full Quarto website project (`_quarto.yml`
at `vignettes/`) would put Quarto in project mode and change the output
path of `quarto render review.qmd --to html` in `Makefile`. `make
review-html` would silently break. Option C avoids this entirely because
nothing changes for the other vignettes.

## What the script will NOT do

- Touch `main` history (the only commit is on the orphan `gh-pages` worktree).
- Push to any branch other than `gh-pages`.
- Add a GitHub Actions workflow.
- Change anything in the existing Makefile review targets.

## Resolved decisions

- Approach: orphan `gh-pages` branch on `origin`, force-pushed.
- Index: Option C (`vignettes/index.qmd`, rendered standalone).
- Separation: render and staging are two scripts and two Makefile
  targets; `gh-pages` does not depend on `render-vignettes`.
- Repo: only the `quant/` subfolder is the site repo (already wired up).
- `RENDER` excludes `diagnostics.qmd` and `swap_visualization.qmd`
  because their YAML defaults are placeholders for the parameterized
  matrix-render pipeline (per-cell `--execute-param` overrides) and
  they fail when rendered standalone. The remaining 7 vignettes plus
  `index.qmd` are rendered best-effort.
- `PUBLISH` is the working subset that renders standalone with no
  parameter overrides: `index`, `review_supplement`,
  `CSF_Spectronaut_data_vis`, `CSF_Spectronaut_sample_swap_data_vis`,
  `CSF_Spectronaut_protein_swap_data_vis`. `review`, `review2`,
  `review_main` are rendered (so their failures surface in the render
  summary) but not on the site.

## Open / to confirm during implementation

- Script location: `src/render_vignettes.sh` and `src/publish_gh_pages.sh`
  (alongside the existing Python helpers in `src/`). Confirm before
  writing if a different location is preferred.
