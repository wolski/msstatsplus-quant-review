#!/usr/bin/env fish

echo "V2 VSN: running Mixture_of_proteomes_processing_v2_vsn.R"
Rscript Mixture_of_proteomes_processing_v2_vsn.R; or exit $status
