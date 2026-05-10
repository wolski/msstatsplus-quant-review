#!/usr/bin/env fish

echo "V1 log2: running CSF_Spectronaut_model_nonmsstats_variant.R log2"
Rscript CSF_Spectronaut_model_nonmsstats_variant.R log2; or exit $status
