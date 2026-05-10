#!/usr/bin/env fish

echo "Step 1/3: running CSF_Spectronaut_processing.R"
Rscript CSF_Spectronaut_processing.R; or exit $status

echo "Step 2/3: running CSF_Spectronaut_analysis.R"
Rscript CSF_Spectronaut_analysis.R; or exit $status

echo "Step 3/3: running CSF_Spectronaut_comparison_table.R"
Rscript CSF_Spectronaut_comparison_table.R; or exit $status
