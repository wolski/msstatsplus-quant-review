# small — TPR / PPV per method × variant

Cells are TPR or PPV at p < 0.05. Column suffix is the data transformation/normalization variant.
`fail` = method errored out;  `see v2_median` = same MSstats internal normalization as v2_median;  `—` = method not run for that variant.

### post-swap

| Method | TPR_V1_log2 | PPV_V1_log2 | TPR_v2_median | PPV_v2_median | TPR_v2_vsn | PPV_v2_vsn | TPR_v3_quantile | PPV_v3_quantile |
|---|---|---|---|---|---|---|---|---|
| MSstats+ | 0.387 | 0.796 | 0.703 | 0.797 | see v2_median | see v2_median | fail | fail |
| MSstats | 0.231 | 0.721 | 0.580 | 0.788 | see v2_median | see v2_median | fail | fail |
| limpa | 0.203 | 1.000 | 0.561 | 0.815 | — | — | 0.495 | 0.868 |
| MaxLFQ + limma | 0.415 | 0.926 | 0.703 | 0.706 | 0.627 | 0.787 | 0.675 | 0.681 |
| msqrob2 | 0.741 | 0.724 | 0.745 | 0.702 | 0.637 | 0.718 | 0.708 | 0.658 |
| DEqMS | 0.321 | 0.883 | 0.571 | 0.665 | 0.613 | 0.663 | 0.642 | 0.654 |
| prolfqua | 0.684 | 0.829 | 0.693 | 0.826 | 0.608 | 0.787 | 0.703 | 0.801 |

### pre-swap

| Method | TPR_V1_log2 | PPV_V1_log2 | TPR_v2_median | PPV_v2_median | TPR_v2_vsn | PPV_v2_vsn | TPR_v3_quantile | PPV_v3_quantile |
|---|---|---|---|---|---|---|---|---|
| MSstats+ | 0.009 | 0.083 | 0.019 | 0.098 | see v2_median | see v2_median | fail | fail |
| MSstats | 0.005 | 0.050 | 0.019 | 0.108 | see v2_median | see v2_median | fail | fail |
| limpa | 0.000 | NA | 0.024 | 0.156 | — | — | 0.009 | 0.111 |
| MaxLFQ + limma | 0.000 | 0.000 | 0.028 | 0.087 | 0.024 | 0.132 | 0.028 | 0.083 |
| msqrob2 | 0.038 | 0.118 | 0.047 | 0.132 | 0.028 | 0.098 | 0.052 | 0.121 |
| DEqMS | 0.000 | 0.000 | 0.024 | 0.075 | 0.038 | 0.107 | 0.042 | 0.118 |
| prolfqua | 0.009 | 0.056 | 0.009 | 0.054 | 0.005 | 0.027 | 0.005 | 0.026 |
