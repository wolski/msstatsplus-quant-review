#!/usr/bin/env fish

echo "Step 1/2: running CSF_Spectronaut_processing.R"
Rscript CSF_Spectronaut_processing.R; or exit $status

echo "Step 2/2: running CSF_Spectronaut_analysis.R"
Rscript CSF_Spectronaut_analysis.R; or exit $status
