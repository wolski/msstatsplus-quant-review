# all_dilutions — TPR / PPV per method × variant

Cells are TPR or PPV at p < 0.05. Column suffix is the data transformation/normalization variant.
`fail` = method errored out;  `see v2_median` = same MSstats internal normalization as v2_median;  `—` = method not run for that variant.

### post-swap, raw p < 0.05

| Method | TPR_V1_log2 | PPV_V1_log2 | TPR_v2_median | PPV_v2_median | TPR_v2_vsn | PPV_v2_vsn | TPR_v3_quantile | PPV_v3_quantile |
|---|---|---|---|---|---|---|---|---|
| MSstats+ | 0.867 | 0.880 | 0.891 | 0.817 | see v2_median | see v2_median | fail | fail |
| MSstats | 0.755 | 0.894 | 0.821 | 0.916 | see v2_median | see v2_median | fail | fail |
| limpa | 0.670 | 0.973 | 0.825 | 0.845 | — | — | 0.717 | 0.916 |
| MaxLFQ + limma | 0.849 | 0.923 | 0.915 | 0.805 | 0.840 | 0.824 | 0.910 | 0.785 |
| msqrob2 | 0.929 | 0.785 | 0.929 | 0.767 | 0.882 | 0.862 | 0.873 | 0.889 |
| DEqMS | 0.783 | 0.907 | 0.854 | 0.823 | 0.854 | 0.733 | 0.882 | 0.716 |
| prolfqua | 0.892 | 0.815 | 0.892 | 0.794 | 0.849 | 0.807 | 0.868 | 0.893 |

### post-swap, FDR (adj.pvalue) < 0.05

| Method | TPR_V1_log2 | PPV_V1_log2 | TPR_v2_median | PPV_v2_median | TPR_v2_vsn | PPV_v2_vsn | TPR_v3_quantile | PPV_v3_quantile |
|---|---|---|---|---|---|---|---|---|
| MSstats+ | 0.678 | 1.000 | 0.806 | 0.977 | see v2_median | see v2_median | fail | fail |
| MSstats | 0.274 | 1.000 | 0.627 | 1.000 | see v2_median | see v2_median | fail | fail |
| limpa | 0.165 | 1.000 | 0.458 | 1.000 | — | — | 0.278 | 0.983 |
| MaxLFQ + limma | 0.472 | 1.000 | 0.797 | 0.988 | 0.675 | 0.993 | 0.731 | 0.987 |
| msqrob2 | 0.840 | 0.978 | 0.854 | 0.984 | 0.750 | 1.000 | 0.608 | 0.985 |
| DEqMS | 0.368 | 1.000 | 0.689 | 0.993 | 0.708 | 0.974 | 0.759 | 0.994 |
| prolfqua | 0.774 | 0.988 | 0.788 | 0.982 | 0.722 | 0.994 | 0.604 | 0.992 |

### pre-swap, raw p < 0.05

| Method | TPR_V1_log2 | PPV_V1_log2 | TPR_v2_median | PPV_v2_median | TPR_v2_vsn | PPV_v2_vsn | TPR_v3_quantile | PPV_v3_quantile |
|---|---|---|---|---|---|---|---|---|
| MSstats+ | 0.019 | 0.133 | 0.028 | 0.136 | see v2_median | see v2_median | fail | fail |
| MSstats | 0.009 | 0.095 | 0.014 | 0.158 | see v2_median | see v2_median | fail | fail |
| limpa | 0.005 | 0.200 | 0.024 | 0.135 | — | — | 0.009 | 0.125 |
| MaxLFQ + limma | 0.005 | 0.062 | 0.024 | 0.100 | 0.033 | 0.121 | 0.024 | 0.091 |
| msqrob2 | 0.024 | 0.083 | 0.019 | 0.063 | 0.028 | 0.162 | 0.014 | 0.052 |
| DEqMS | 0.000 | 0.000 | 0.028 | 0.130 | 0.038 | 0.108 | 0.033 | 0.090 |
| prolfqua | 0.009 | 0.043 | 0.009 | 0.038 | 0.005 | 0.023 | 0.014 | 0.061 |

### pre-swap, FDR (adj.pvalue) < 0.05

| Method | TPR_V1_log2 | PPV_V1_log2 | TPR_v2_median | PPV_v2_median | TPR_v2_vsn | PPV_v2_vsn | TPR_v3_quantile | PPV_v3_quantile |
|---|---|---|---|---|---|---|---|---|
| MSstats+ | 0.000 | NA | 0.000 | NA | see v2_median | see v2_median | fail | fail |
| MSstats | 0.000 | NA | 0.000 | NA | see v2_median | see v2_median | fail | fail |
| limpa | 0.000 | NA | 0.000 | NA | — | — | 0.000 | NA |
| MaxLFQ + limma | 0.000 | NA | 0.000 | NA | 0.000 | NA | 0.000 | NA |
| msqrob2 | 0.000 | NA | 0.000 | NA | 0.000 | NA | 0.000 | NA |
| DEqMS | 0.000 | NA | 0.000 | NA | 0.000 | NA | 0.000 | NA |
| prolfqua | 0.000 | NA | 0.000 | NA | 0.000 | NA | 0.000 | NA |
