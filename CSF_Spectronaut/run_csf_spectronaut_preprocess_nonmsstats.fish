#!/usr/bin/env fish

echo "Preprocessing non-MSstats CSF Spectronaut inputs"
Rscript CSF_Spectronaut_preprocess_nonmsstats.R; or exit $status
