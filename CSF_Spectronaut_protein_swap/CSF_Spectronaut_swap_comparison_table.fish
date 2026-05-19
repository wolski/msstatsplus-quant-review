#!/usr/bin/env fish
# Build the comparison_table.{csv,txt} for each named run directory.
# Pure summary step — reads existing model + timing CSVs, doesn't re-fit.
#
# Usage:
#   ./CSF_Spectronaut_swap_comparison_table.fish [TAG ...]
#
# Defaults to building tables for both all_dilutions and no_high_dilutions.
# Examples:
#   ./CSF_Spectronaut_swap_comparison_table.fish
#   ./CSF_Spectronaut_swap_comparison_table.fish no_high_dilutions
#   ./CSF_Spectronaut_swap_comparison_table.fish all_dilutions custom_run

set TAGS $argv
if test (count $TAGS) -eq 0
    set TAGS all_dilutions no_high_dilutions
end

for tag in $TAGS
    if not test -d $tag
        echo "[skip] $tag/ does not exist" >&2
        continue
    end
    echo "[build] $tag/comparison_table.{csv,txt}"
    env OUT_TAG=$tag Rscript CSF_Spectronaut_swap_comparison_table.R
end
