#!/usr/bin/env fish

echo "Preprocessing non-MSstats K562 Spectronaut inputs"
Rscript K562_Spectronaut_preprocess_nonmsstats.R; or exit $status
