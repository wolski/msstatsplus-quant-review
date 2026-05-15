#!/usr/bin/env fish

echo "V1 log2: running Mixture_of_proteomes_processing_v1_log2.R"
Rscript Mixture_of_proteomes_processing_v1_log2.R; or exit $status
