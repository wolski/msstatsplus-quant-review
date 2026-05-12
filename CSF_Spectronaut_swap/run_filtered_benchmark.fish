#!/usr/bin/env fish
# Run the precursor-swap benchmark with a custom dilution filter.
#
# Usage:
#   ./run_filtered_benchmark.fish <OUT_TAG> <EXCLUDE_DILUTIONS>
#
# Examples:
#   ./run_filtered_benchmark.fish _no_high_dilutions "1to32,1to64"
#
# The tag is appended to the variant directories (V1_log2<TAG>/ and
# v2_vsn<TAG>/) so canonical outputs are not overwritten. Pass the same
# OUT_TAG to CSF_Spectronaut_swap_comparison_table.R afterwards.

if test (count $argv) -lt 2
    echo "Usage: $argv[0] <OUT_TAG> <EXCLUDE_DILUTIONS>" >&2
    echo "  e.g.  $argv[0] _no_high_dilutions \"1to32,1to64\"" >&2
    exit 2
end

set TAG    $argv[1]
set EXCL   $argv[2]
set SWAP   "20250130_163144_CSF dilutions Jan 2025 no normalization_Report.tsv"
set ORIG   "../CSF_Spectronaut/20250130_163144_CSF dilutions Jan 2025 no normalization_Report.tsv"

echo "TAG=$TAG  EXCLUDE_DILUTIONS=$EXCL"

set BASE_ENV "OUT_TAG=$TAG" "EXCLUDE_DILUTIONS=$EXCL"

# Step 1: V1_log2 post-swap (full main script) ----------------------
env $BASE_ENV REPORT_PATH=$SWAP VARIANT=V1_log2 OUT_SUFFIX="" \
    Rscript CSF_Spectronaut_swap_processing.R > run_V1_log2$TAG.log 2>&1
or exit $status

# Step 2: V1_log2 pre-swap ------------------------------------------
env $BASE_ENV REPORT_PATH=$ORIG VARIANT=V1_log2 OUT_SUFFIX="_preswap" \
    Rscript CSF_Spectronaut_swap_processing.R > run_V1_log2_preswap$TAG.log 2>&1
or exit $status

# Step 3: v2_vsn post-swap (skips MSstats; main script branches on VARIANT) -
env $BASE_ENV REPORT_PATH=$SWAP VARIANT=v2_vsn OUT_SUFFIX="" \
    Rscript CSF_Spectronaut_swap_processing.R > run_v2_vsn$TAG.log 2>&1
or exit $status

# Step 4: v2_vsn pre-swap -------------------------------------------
env $BASE_ENV REPORT_PATH=$ORIG VARIANT=v2_vsn OUT_SUFFIX="_preswap" \
    Rscript CSF_Spectronaut_swap_processing.R > run_v2_vsn_preswap$TAG.log 2>&1
or exit $status

# Step 5: comparison table -------------------------------------------
env $BASE_ENV Rscript CSF_Spectronaut_swap_comparison_table.R
