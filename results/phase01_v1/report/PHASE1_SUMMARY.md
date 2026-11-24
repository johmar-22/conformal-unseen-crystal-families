# Phase 1 summary (alpha = 0.1)

Repeats found: 25  |  units: 54

## bandgap: mean coverage over repeats

| group   |   cvplus_lattice |   cvplus_random |   split_global |   split_mondrian |   split_normalized_knn |   cqr |
|:--------|-----------------:|----------------:|---------------:|-----------------:|-----------------------:|------:|
| ALL     |            0.985 |           0.907 |          0.887 |            0.897 |                  0.892 | 0.897 |
| C2/m    |            1     |           0.978 |          0.96  |            0.89  |                  0.967 | 0.881 |
| Pna21   |            1     |           0.998 |          0.995 |            0.913 |                  0.99  | 0.932 |
| R-3c    |            0.997 |           0.964 |          0.949 |            0.89  |                  0.952 | 0.868 |
| P63/mmc |            0.988 |           0.905 |          0.867 |            0.902 |                  0.86  | 0.901 |
| Ia-3    |            0.996 |           0.983 |          0.974 |            0.885 |                  0.971 | 0.935 |
| Fd-3m   |            0.932 |           0.6   |          0.543 |            0.898 |                  0.577 | 0.859 |

## bandgap: median interval width

| group   |   cvplus_lattice |   cvplus_random |   split_global |   split_mondrian |   split_normalized_knn |   cqr |
|:--------|-----------------:|----------------:|---------------:|-----------------:|-----------------------:|------:|
| ALL     |            1.43  |           0.657 |          0.626 |            0.406 |                  0.599 | 0.994 |
| C2/m    |            1.535 |           0.653 |          0.626 |            0.398 |                  0.554 | 0.574 |
| Pna21   |            1.376 |           0.656 |          0.626 |            0.326 |                  0.54  | 0.84  |
| R-3c    |            1.585 |           0.666 |          0.626 |            0.421 |                  0.755 | 1.19  |
| P63/mmc |            1.403 |           0.668 |          0.626 |            0.718 |                  0.586 | 1.055 |
| Ia-3    |            1.386 |           0.653 |          0.626 |            0.275 |                  0.502 | 0.817 |
| Fd-3m   |            1.654 |           0.688 |          0.626 |            1.361 |                  0.649 | 1.369 |

## formation: mean coverage over repeats

| group   |   cvplus_lattice |   cvplus_random |   split_global |   split_mondrian |   split_normalized_knn |   cqr |
|:--------|-----------------:|----------------:|---------------:|-----------------:|-----------------------:|------:|
| ALL     |            0.997 |           0.919 |          0.896 |            0.894 |                  0.893 | 0.903 |
| C2/m    |            1     |           0.968 |          0.939 |            0.881 |                  0.929 | 0.878 |
| Pna21   |            1     |           0.989 |          0.979 |            0.911 |                  0.966 | 0.926 |
| R-3c    |            1     |           0.975 |          0.952 |            0.902 |                  0.987 | 0.904 |
| P63/mmc |            0.993 |           0.942 |          0.924 |            0.902 |                  0.916 | 0.898 |
| Ia-3    |            0.997 |           0.992 |          0.982 |            0.854 |                  0.965 | 0.93  |
| Fd-3m   |            0.989 |           0.646 |          0.582 |            0.906 |                  0.586 | 0.87  |

## formation: median interval width

| group   |   cvplus_lattice |   cvplus_random |   split_global |   split_mondrian |   split_normalized_knn |   cqr |
|:--------|-----------------:|----------------:|---------------:|-----------------:|-----------------------:|------:|
| ALL     |            0.333 |           0.096 |          0.091 |            0.068 |                  0.082 | 0.133 |
| C2/m    |            0.357 |           0.096 |          0.091 |            0.066 |                  0.08  | 0.144 |
| Pna21   |            0.326 |           0.096 |          0.091 |            0.056 |                  0.078 | 0.106 |
| R-3c    |            0.336 |           0.096 |          0.091 |            0.077 |                  0.103 | 0.11  |
| P63/mmc |            0.338 |           0.097 |          0.091 |            0.084 |                  0.081 | 0.176 |
| Ia-3    |            0.327 |           0.095 |          0.091 |            0.036 |                  0.069 | 0.083 |
| Fd-3m   |            0.365 |           0.1   |          0.091 |            0.176 |                  0.091 | 0.233 |

## Pre-registered endpoint

```
{
  "bandgap": {
    "primary_endpoint_fd3m_cov_mean": 0.6002819628609325,
    "decision_H1a": "H1a supported",
    "H1b_supported": true
  },
  "formation": {
    "primary_endpoint_fd3m_cov_mean": 0.6464319435697737,
    "decision_H1a": "H1a supported",
    "H1b_supported": false
  }
}
```

## Hyperparameter sensitivity (own vs split-0 tuning)

|                                       |   mean |    min |   max |
|:--------------------------------------|-------:|-------:|------:|
| ('bandgap', 'cqr')                    |  0.012 | -0.023 | 0.062 |
| ('bandgap', 'cvplus_lattice')         | -0.005 | -0.033 | 0     |
| ('bandgap', 'cvplus_random')          |  0     | -0.014 | 0.022 |
| ('bandgap', 'split_global')           |  0.022 |  0     | 0.109 |
| ('bandgap', 'split_mondrian')         | -0.004 | -0.076 | 0.062 |
| ('bandgap', 'split_normalized_knn')   | -0.003 | -0.044 | 0.028 |
| ('formation', 'cqr')                  | -0.008 | -0.043 | 0.027 |
| ('formation', 'cvplus_lattice')       |  0     |  0     | 0     |
| ('formation', 'cvplus_random')        |  0.002 |  0     | 0.031 |
| ('formation', 'split_global')         |  0.002 | -0.014 | 0.031 |
| ('formation', 'split_mondrian')       |  0.005 | -0.094 | 0.058 |
| ('formation', 'split_normalized_knn') |  0.002 | -0.031 | 0.031 |
