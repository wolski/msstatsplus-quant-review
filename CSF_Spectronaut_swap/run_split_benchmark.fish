#!/usr/bin/env fish
# Run the precursor-swap benchmark with the split pipeline (run_msstats.R
# and run_nonmsstats.R in parallel within each variant×swap-state cell).
#
# Usage:
#   ./run_split_benchmark.fish [OUT_TAG] [EXCLUDE_DILUTIONS]
#
# OUT_TAG is the parent run directory (e.g. "all_dilutions" or
# "no_high_dilutions"); outputs land in <OUT_TAG>/<variant>/<method>{,_preswap}/.
#
# Examples:
#   ./run_split_benchmark.fish all_dilutions ""
#   ./run_split_benchmark.fish no_high_dilutions "1to32,1to64"
#
# Variants:
#   V1_log2     : NORMALIZATION=none for both scripts.
#   v2_vsn      : MSstats=equalizeMedians; non-MSstats=vsn.
#   v3_quantile : MSstats=quantile; non-MSstats=quantile (log2 + limma's
#                 normalizeBetweenArrays); MSstats's built-in "quantile"
#                 mode matches that on its summarized intensities.

set TAG  "all_dilutions"
set EXCL ""
if test (count $argv) -ge 1
    set TAG $argv[1]
end
if test (count $argv) -ge 2
    set EXCL $argv[2]
end
mkdir -p $TAG

set SWAP "20250130_163144_CSF dilutions Jan 2025 no normalization_Report.tsv"
set ORIG "../CSF_Spectronaut/20250130_163144_CSF dilutions Jan 2025 no normalization_Report.tsv"

echo "TAG=$TAG  EXCLUDE_DILUTIONS=$EXCL"

# Run one variant×swap-state cell: MSstats and non-MSstats in parallel.
function run_cell --argument variant suffix report norm_msstats norm_nonms
    set ms_log  $TAG/run_(string trim "$variant")_msstats(string trim "$suffix").log
    set non_log $TAG/run_(string trim "$variant")_nonmsstats(string trim "$suffix").log
    echo "  -> $ms_log  +  $non_log  (parallel)"

    env OUT_TAG=$TAG EXCLUDE_DILUTIONS=$EXCL REPORT_PATH=$report \
        VARIANT=$variant OUT_SUFFIX=$suffix NORMALIZATION=$norm_msstats \
        Rscript run_msstats.R > $ms_log 2>&1 &
    set ms_pid $last_pid

    env OUT_TAG=$TAG EXCLUDE_DILUTIONS=$EXCL REPORT_PATH=$report \
        VARIANT=$variant OUT_SUFFIX=$suffix NORMALIZATION=$norm_nonms \
        Rscript run_nonmsstats.R > $non_log 2>&1 &
    set non_pid $last_pid

    wait $ms_pid; set ms_status $status
    wait $non_pid; set non_status $status

    if test $ms_status -ne 0
        echo "  ! MSstats failed (status $ms_status) — see $ms_log" >&2
    end
    if test $non_status -ne 0
        echo "  ! non-MSstats failed (status $non_status) — see $non_log" >&2
    end
end

echo "== V1_log2     post-swap"
run_cell V1_log2     ""         $SWAP none            none
echo "== V1_log2     pre-swap"
run_cell V1_log2     "_preswap" $ORIG none            none
echo "== v2_vsn      post-swap"
run_cell v2_vsn      ""         $SWAP equalizeMedians vsn
echo "== v2_vsn      pre-swap"
run_cell v2_vsn      "_preswap" $ORIG equalizeMedians vsn
echo "== v3_quantile post-swap"
run_cell v3_quantile ""         $SWAP quantile        quantile
echo "== v3_quantile pre-swap"
run_cell v3_quantile "_preswap" $ORIG quantile        quantile

echo "Building comparison table..."
env OUT_TAG=$TAG Rscript CSF_Spectronaut_swap_comparison_table.R
