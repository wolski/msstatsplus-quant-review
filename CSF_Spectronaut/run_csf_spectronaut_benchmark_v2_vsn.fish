#!/usr/bin/env fish

echo "V2 VSN: running CSF_Spectronaut_model_nonmsstats_variant.R vsn"
Rscript CSF_Spectronaut_model_nonmsstats_variant.R vsn; or exit $status
