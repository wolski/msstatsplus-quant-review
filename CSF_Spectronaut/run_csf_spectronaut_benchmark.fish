#!/usr/bin/env fish

# Execution order for the CSF Spectronaut benchmark:
# 1. CSF_Spectronaut_processing.R writes method inputs, model tables, and RDA files.
# 2. CSF_Spectronaut_analysis.R reads those outputs to compute benchmark summaries/plots.
#
# The upstream R scripts are intentionally left unchanged. Temporary runnable copies
# are created because the restored scripts expect both data_folder="" and a helper
# source file in the current working directory.

set script_dir (realpath (dirname (status --current-filename)))
set helper "$script_dir/../benchmark_experiments_functions.R"
set tmp_dir "$script_dir/.runner_tmp"

if not test -f "$helper"
    echo "Missing helper file: $helper" >&2
    exit 1
end

mkdir -p "$tmp_dir"; or exit $status
cd "$script_dir"; or exit $status

function prepare_r_script
    set src $argv[1]
    set dst $argv[2]

    sed \
        -e 's/source("benchmark_experiments_functions.R")/source("..\/benchmark_experiments_functions.R")/' \
        -e 's/data_folder = ""/data_folder = "."/' \
        "$src" > "$dst"
end

set processing_tmp "$tmp_dir/01_CSF_Spectronaut_processing.R"
set analysis_tmp "$tmp_dir/02_CSF_Spectronaut_analysis.R"

prepare_r_script "$script_dir/CSF_Spectronaut_processing.R" "$processing_tmp"; or exit $status
prepare_r_script "$script_dir/CSF_Spectronaut_analysis.R" "$analysis_tmp"; or exit $status

echo "Step 1/2: running CSF_Spectronaut_processing.R"
Rscript "$processing_tmp"; or exit $status

echo "Step 2/2: running CSF_Spectronaut_analysis.R"
Rscript "$analysis_tmp"; or exit $status
