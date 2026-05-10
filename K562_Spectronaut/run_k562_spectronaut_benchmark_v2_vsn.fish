#!/usr/bin/env fish

echo "K562 Spectronaut V2 VSN: running K562_Spectronaut_model_nonmsstats_variant.R vsn"
Rscript K562_Spectronaut_model_nonmsstats_variant.R vsn; or exit $status
