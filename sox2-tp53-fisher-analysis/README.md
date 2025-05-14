# SOX2 Copy Number Variation and TP53 Mutation Analysis

This project investigates the potential association between **SOX2 copy number variation** and **TP53 mutation status** across multiple cancer types using **Fisher’s Exact Test**. The analysis is performed in R using data obtained from [cBioPortal](https://www.cbioportal.org), with results visualized using `ggstatsplot`, base R plots, and `RColorBrewer`.


## 📊 Objective

To statistically assess whether **SOX2 gene amplification (Gain)** is associated with **TP53 missense mutations** in cancer patients.

## 📦 Tools & Packages

- `ggstatsplot`
- `RColorBrewer`
- Base R plotting
- Fisher’s Exact Test (`fisher.test()`)

## 📂 Data Source

- Dataset: `pancancer_combined.txt`
- Origin: [cBioPortal](https://www.cbioportal.org)
- Columns:  
  - `ID`: Unique patient identifier  
  - `TP53_status`: Mutation status (e.g., "Missense", "No mutation")  
  - `SOX2_status`: Copy number variation status (e.g., "Gain", "Diploid")

## 🔬 Analysis Summary

1. **Data cleaning** and renaming of columns for readability.
2. **Subsetting** the dataset based on mutation and amplification status.
3. **Construction of a 2x2 contingency table**.
4. **Statistical testing** using Fisher’s Exact Test.
5. **Visualization** using mosaic plots, bar charts, and `ggbarstats()`.

## 📈 Key Result

- **P-value**: `< 0.001`
- **Odds ratio**: `2.47`
- **Interpretation**: There is a statistically significant association between SOX2 gain and TP53 mutation across the cancer types analyzed.

## 📌 How to Reproduce

1. Clone this repository or download the files.
2. Open `SOX2_TP53_analysis.Rmd` in RStudio.
3. Run all chunks or knit to HTML to reproduce the analysis and visualizations.

## 📚 References

1. Oliveira, N.L., et al. (2018). *A discussion on significance indices for contingency tables under small sample sizes*. PLoS One, 13(8), e0199102.
2. Kim, H.-Y. (2017). *Statistical notes for clinical researchers: Chi-squared test and Fisher’s exact test*. Restorative Dentistry & Endodontics, 42(2), 152–155.

---