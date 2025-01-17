# BSGS: Bayesian Subset-based Gene Selection  
**A Stochastic Multi-Criteria Approach for High-Dimensional Biomarker Identification**

## Introduction  
The **Bayesian Subset-based Gene Selection (BSGS)** framework is an advanced, customizable tool designed for identifying biomarkers in high-dimensional datasets. By integrating Bayesian subset sampling, Lasso-based preselection, and multi-criteria scoring, BSGS delivers a robust methodology for biomarker discovery.

This R package is specifically crafted to enhance user experience and reproducibility in biomarker research, making it easy to apply and adapt for different datasets.

---

## Key Features  
- **Lasso-based Preselection**: Efficiently reduces the feature space by selecting genes with high importance.  
- **Bayesian Subset Sampling**: Iteratively updates posterior probabilities to refine gene scores.  
- **Customizable Scoring**: Combines statistical metrics (e.g., SNR, variance, correlations) and SHAP values for multi-criteria ranking.  
- **Visualization**: Automatically generates distribution plots for gene scores.  
- **Reproducibility**: Offers clear workflows for reproducible research.

---

## Installation  

To get started, install the BSGS package directly from GitHub:

```r
# Install the devtools package if not already installed
install.packages("devtools")

# Install BSGS from GitHub
devtools::install_github("PaplomatasP/BSGS")

# Load the BSGS package
library(BSGS)

# Load the included example dataset
data(LUAD_MYE_MRS_GSE97168)

# Explore the dataset
str(LUAD_MYE_MRS_GSE97168)

# Convert the dataset to a data frame
full_data <- as.data.frame(LUAD_MYE_MRS_GSE97168)

# Add target labels
full_data$target <- as.factor(LUAD_MYE_MRS_target)

# Run the BSGS framework
bsgs_results <- BSGS(
  data = full_data,
  label_col = "target",
  iterations = 1000,
  alpha = 0.3,
  min_genes = round(floor(sqrt(length(colnames(full_data[,-1])))) / 2),
  max_genes = 5 * round(floor(sqrt(length(colnames(full_data[,-1])))), 2),
  n_genes = 150,
  plot_distribution = TRUE,
  save_plot = TRUE,
  file_path = "./bsgs_results.png"
)

# View selected genes
print(bsgs_results$final_genes)

# View gene scores
print(bsgs_results$final_scores)

```



## Gene Score Distribution Plot

When `plot_distribution = TRUE` is set in the **BSGS** function, a gene score distribution plot is displayed to the user. This plot visualizes the scores of all genes analyzed, allowing the user to identify meaningful thresholds for selecting the final genes based on their scores. The user is prompted to provide a threshold value, and genes with scores above this threshold are retained.

Below is an example of the gene score distribution plot:

```{r, echo=FALSE, out.width="100%", fig.align='center'}
knitr::include_graphics("Experiment/Results_from_the_Bayesian_Gene_Scoring_System_(BSGS)_20000.png")
