#!/usr/bin/env fish
# Run the precursor-swap benchmark.
#   V1_log2: log2 normalization, 7 methods (MSstats+/MSstats/msqrob2/limma/
#            limpa/DEqMS/prolfqua), post+pre swap.
#   v2_vsn:  vsn::justvsn normalization, 5 methods (skip MSstats+/MSstats),
#            post+pre swap.
# V1 and V2 run in parallel; each post-swap → pre-swap sequentially.

set SWAP_REPORT "20250130_163144_CSF dilutions Jan 2025 no normalization_Report.tsv"
set ORIG_REPORT "../CSF_Spectronaut/20250130_163144_CSF dilutions Jan 2025 no normalization_Report.tsv"

function run_variant
    set variant $argv[1]
    set logfile run_$variant.log
    echo "[$variant] post-swap" > $logfile
    env REPORT_PATH=$SWAP_REPORT VARIANT=$variant OUT_SUFFIX="" \
        Rscript CSF_Spectronaut_swap_processing.R >> $logfile 2>&1
    or return $status
    echo "[$variant] pre-swap" >> $logfile
    env REPORT_PATH=$ORIG_REPORT VARIANT=$variant OUT_SUFFIX="_preswap" \
        Rscript CSF_Spectronaut_swap_processing.R >> $logfile 2>&1
end

echo "Launching V1_log2 and v2_vsn in parallel..."
run_variant V1_log2 &
set v1_pid $last_pid
run_variant v2_vsn &
set v2_pid $last_pid

wait $v1_pid; set v1_status $status
wait $v2_pid; set v2_status $status

if test $v1_status -ne 0
    echo "V1_log2 failed (status $v1_status)" >&2
    exit $v1_status
end
if test $v2_status -ne 0
    echo "v2_vsn failed (status $v2_status)" >&2
    exit $v2_status
end

echo "Building comparison table..."
Rscript CSF_Spectronaut_swap_comparison_table.R
