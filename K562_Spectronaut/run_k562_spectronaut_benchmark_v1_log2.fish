#!/usr/bin/env fish

echo "K562 Spectronaut V1 log2: running K562_Spectronaut_model_nonmsstats_variant.R log2"
Rscript K562_Spectronaut_model_nonmsstats_variant.R log2; or exit $status
