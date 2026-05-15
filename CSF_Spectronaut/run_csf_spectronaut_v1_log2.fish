#!/usr/bin/env fish
# Canonical V1/log2 runner for the original CSF Spectronaut benchmark.
#
# Usage:
#   ./run_csf_spectronaut_v1_log2.fish [OUT_TAG] [EXCLUDE_DILUTIONS]
#
# Outputs land in <OUT_TAG>/V1_log2/<method>/.

set script_dir (dirname (status --current-filename))
cd $script_dir

set tag "all_dilutions"
set exclude_dilutions ""

if test (count $argv) -ge 1
    set tag $argv[1]
end
if test (count $argv) -ge 2
    set exclude_dilutions $argv[2]
end

mkdir -p $tag

function run_step --argument-names label script log_file
    echo "== $label"
    echo "   log: $log_file"
    env OUT_TAG=$tag EXCLUDE_DILUTIONS=$exclude_dilutions \
        VARIANT=V1_log2 NORMALIZATION=none \
        Rscript $script > $log_file 2>&1

    set step_status $status
    if test $step_status -ne 0
        echo "!! $label failed with status $step_status; see $log_file" >&2
        exit $step_status
    end
end

echo "OUT_TAG=$tag"
if test -n "$exclude_dilutions"
    echo "EXCLUDE_DILUTIONS=$exclude_dilutions"
else
    echo "EXCLUDE_DILUTIONS=(none)"
end

run_step "MSstats/MSstats+" \
    run_msstats.R \
    $tag/run_V1_log2_msstats.log

run_step "MSstats/MSstats+ pre-swap" \
    run_msstats_preswap.R \
    $tag/run_V1_log2_msstats_preswap.log

run_step "non-MSstats" \
    run_nonmsstats.R \
    $tag/run_V1_log2_nonmsstats.log

echo "V1/log2 run finished."
