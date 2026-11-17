# CD4/CD8 ratio and mathematical coupling: Phase 2 re-analysis

Data and R code for Phase 2 of the HDAT9901 dissertation (UNSW Sydney, 2026), *Assessing the Impact of Mathematical Coupling on CD4/CD8 Ratio Studies in HIV: A Systematic Review and Re-analysis*. The code reproduces Tables 4 and 5 of the dissertation.

## Files

| File | What it is |
| --- | --- |
| `data/phase2_data.csv` | The data. One row per study group: the published values used, and where each came from (paper DOI, table or page). |
| `phase2_analysis.R` | The analysis, with comments. Builds Table 4 from the data file and calculates every value in Table 5. |
| `cd4cd8-coupling-reanalysis.Rproj` | RStudio project file (optional). |

Running the script creates a `results/` folder: `table4.csv`, `phase2_results.csv` (Table 5), the sensitivity and rounding checks, and `session_info.txt`.

## Data

**These are real published values, not simulated data.** All values are summary statistics (means and SDs, or medians and IQRs) copied from tables and text in five published papers. No values were read from figures, and there are no individual-level data. Every value was checked a second time against its source.

| Study | DOI | Source in the paper |
| --- | --- | --- |
| Monsalvo et al. 2019 | 10.1177/0956462419834129 | p. 659 |
| Taramasso et al. 2022 | 10.3390/v14112315 | Tables 1 and 2 |
| Trizzino et al. 2025 | 10.3390/v17111449 | Table 3 |
| Gidari et al. 2023 | 10.3390/idr15060069 | Table 4 |
| Ergün et al. 2025 | 10.1177/23259582251404545 | Table 6 |

Section 2 of the script also makes a small simulated dataset (fixed random seed). It is used only to check that the formulas give the same answers as calculating from individual values. It is not used for any result.

## How to run

1. Download the repository (green **Code** button, then **Download ZIP**) and unzip it.
2. Open `cd4cd8-coupling-reanalysis.Rproj` in RStudio and run `source("phase2_analysis.R")`. Or, from the repository folder, run `Rscript phase2_analysis.R`.

The script prints Table 4, the Table 5 results and a check of 31 values against the dissertation (OK or CHECK for each). It stops with an error if any internal check fails.

## Software

- R 4.3.3 (tested).
- No add-on packages. Only base R is used (the `stats` and `utils` packages that come with R).
- Each run writes the R version and session details to `results/session_info.txt`.

## Methods in brief

s1 = baseline SD, s2 = follow-up SD, r = correlation between baseline and follow-up, n = number of people.

- **Variance ratio:** s2² / s1².
- **Oldham's method:** correlation between the change and the average of baseline and follow-up.
- **Pitman–Morgan test:** tests whether the baseline and follow-up variances are equal; t on n − 2 degrees of freedom.
- **Unreported r:** recovered from the SD of the change where reported (Monsalvo 2019); otherwise varied from 0.3 to 0.9.
- **Median and IQR:** SDs estimated with Wan et al. (2014).
- **Rounding checks:** Gidari 2023 (SDs to one decimal) and Monsalvo 2019 (125 combinations of the three IQR widths).

## References

- Morgan WA (1939). A test for the significance of the difference between the two variances in a sample from a normal bivariate population. *Biometrika* 31:13–19.
- Oldham PD (1962). A note on the analysis of repeated measurements of the same subjects. *Journal of Chronic Diseases* 15:969–977.
- Pitman EJG (1939). A note on normal correlation. *Biometrika* 31:9–12.
- Tu YK, Gilthorpe MS (2007). Revisiting the relation between change and initial value: a review and evaluation. *Statistics in Medicine* 26:443–457.
- Wan X, Wang W, Liu J, Tong T (2014). Estimating the sample mean and standard deviation from the sample size, median, range and/or interquartile range. *BMC Medical Research Methodology* 14:135.
