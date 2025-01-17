########################
# Support Functions
########################

#' Normalize Gene Expression Data
#'
#' Applies min-max normalization and power scaling to gene expression data.
#'
#' @param x Numeric vector of gene expression data.
#' @param gamma Power scaling factor (default = 3).
#' @return A numeric vector with normalized values.
normalize <- function(x, gamma = 3) {
  x_no_na <- x
  x_no_na[is.na(x_no_na)] <- 0

  if (length(x_no_na) == 0 || max(x_no_na) == min(x_no_na)) {
    return(rep(0, length(x_no_na)))
  }

  x_norm <- (x_no_na - min(x_no_na)) / (max(x_no_na) - min(x_no_na))
  x_power <- x_norm ^ gamma
  return(x_power)
}

#' Sample Genes Based on Posterior Probabilities
#'
#' Selects a subset of genes based on posterior probabilities.
#'
#' @param posteriors Numeric vector of posterior probabilities.
#' @param num_samples Number of samples to select.
#' @return Indices of selected genes.
sample_genes <- function(posteriors, num_samples) {
  eps <- 1e-15
  positive_indices <- which(posteriors > eps)
  num_samples <- min(num_samples, length(positive_indices))
  sample(seq_along(posteriors), size = num_samples, prob = posteriors, replace = FALSE)
}

#' Update Posterior Probabilities
#'
#' Updates posterior probabilities using priors and likelihoods.
#'
#' @param priors Numeric vector of prior probabilities.
#' @param likelihoods Numeric vector of likelihood values.
#' @param floor_val Minimum allowable posterior value (default = 1e-6).
#' @return Updated posterior probabilities.
update_posteriors <- function(priors, likelihoods, floor_val = 1e-6) {
  likelihoods <- likelihoods + 1e-15
  posterior <- priors * likelihoods

  if (sum(posterior) == 0) {
    warning("All posterior values are zero. Using uniform priors.")
    posterior <- rep(1 / length(priors), length(priors))
  } else {
    posterior[posterior < floor_val] <- floor_val
    posterior <- posterior / sum(posterior)
  }
  return(posterior)
}

#' Compute Group Score and Per-Gene Scores
#'
#' Calculates group scores and per-gene scores based on statistical measures.
#'
#' @param selected_genes Character vector of selected gene names.
#' @param data Data frame containing gene expression data.
#' @param label_col Column name of labels in the data.
#' @param shap_values Numeric vector of SHAP values for the genes.
#' @return List with group score and per-gene scores.
compute_group_score <- function(selected_genes, data, label_col, shap_values) {
  if (length(selected_genes) == 0) {
    warning("No genes selected, returning default likelihood 1e-6")
    return(list(group_score = 1e-6, per_gene_scores = numeric(0)))
  }

  labels <- data[[label_col]]
  selected_data <- data[, colnames(data) %in% selected_genes, drop = FALSE]

  if (is.null(selected_data) || ncol(selected_data) == 0) {
    warning("No matching genes found in data, returning default likelihood 1e-6")
    return(list(group_score = 1e-6, per_gene_scores = numeric(0)))
  }

  selected_data[[label_col]] <- labels

  unique_labels <- unique(labels)
  if (length(unique_labels) < 2) {
    warning("Not enough unique groups in label_col, returning default likelihood 1e-6")
    return(list(group_score = 1e-6, per_gene_scores = numeric(0)))
  }

  group_means <- sapply(unique_labels, function(g) {
    colMeans(selected_data[selected_data[[label_col]] == g, -ncol(selected_data), drop = FALSE],
             na.rm = TRUE)
  })

  group_sds <- sapply(unique_labels, function(g) {
    apply(selected_data[selected_data[[label_col]] == g, -ncol(selected_data), drop = FALSE],
          2, sd, na.rm = TRUE)
  })

  snr_vector <- abs(group_means[, 1] - group_means[, 2]) / (group_sds[, 1] + group_sds[, 2] + 1e-9)
  gene_correlations <- apply(selected_data[, -ncol(selected_data), drop = FALSE], 2, function(g) {
    cor(g, as.numeric(labels), use = "complete.obs")
  })

  gene_variances <- apply(selected_data[, -ncol(selected_data), drop = FALSE], 2, var, na.rm = TRUE)

  norm_snr <- normalize(snr_vector)
  norm_cor <- normalize(abs(gene_correlations))
  norm_var <- normalize(gene_variances)

  gene_weights <- shap_values[selected_genes]
  gene_weights[is.na(gene_weights)] <- 0
  norm_shap <- normalize(gene_weights)

  subset_scores <- 0.3 * norm_snr + 0.2 * norm_cor + 0.4 * norm_var + 0.2 * norm_shap

  group_score <- mean(subset_scores, na.rm = TRUE)
  if (is.na(group_score) || group_score <= 0) {
    group_score <- 1e-6
  }

  return(list(
    group_score = group_score,
    per_gene_scores = subset_scores
  ))
}

#' Lasso Preselection of Genes
#'
#' Performs Lasso regression for gene selection and assigns SHAP-like weights.
#'
#' @param data A data frame containing gene expression data.
#' @param label_col The name of the column with labels for classification.
#' @param alpha Regularization parameter (default = 1 for Lasso).
#' @return A list with selected genes and their corresponding SHAP-like weights.
#' @examples
#' \dontrun{
#'   lasso_results <- lasso_preselection(data, label_col = "labels", alpha = 1)
#'   print(lasso_results$genes)
#'   print(lasso_results$shap_values)
#' }
lasso_preselection <- function(data, label_col, alpha = 1) {
  x <- as.matrix(data[, -which(colnames(data) == label_col), drop = FALSE])
  y <- as.numeric(data[[label_col]])

  if (ncol(x) == 0 || length(unique(y)) < 2) {
    warning("Data not suitable for Lasso. Returning all genes.")
    return(list(
      genes = colnames(x),
      shap_values = setNames(rep(0, ncol(x)), colnames(x))
    ))
  }

  lasso_model <- tryCatch({
    glmnet::cv.glmnet(x, y, alpha = alpha, family = "binomial", standardize = TRUE)
  }, error = function(e) {
    warning("Lasso model failed: ", e$message)
    return(NULL)
  })

  if (is.null(lasso_model)) {
    return(list(
      genes = colnames(x),
      shap_values = setNames(rep(0, ncol(x)), colnames(x))
    ))
  }

  selected_raw <- coef(lasso_model, s = "lambda.min")
  selected_raw <- as.matrix(selected_raw)

  # Αφαιρούμε το intercept
  beta <- selected_raw[-1, , drop = FALSE]
  rownames_beta <- rownames(beta)

  # Ποια έχουν μη-μηδενικό βάρος
  nonzero_indices <- which(beta != 0)
  if (length(nonzero_indices) == 0) {
    warning("No nonzero coefficients found by Lasso.")
    return(list(genes = character(0), shap_values = numeric(0)))
  }

  selected_genes <- rownames_beta[nonzero_indices]
  selected_scores <- abs(beta[nonzero_indices])
  total_score <- sum(selected_scores, na.rm = TRUE)

  shap_vec <- selected_scores / total_score
  names(shap_vec) <- selected_genes

  # Επιστρέφουμε μόνο τα επιλεγμένα γονίδια
  return(list(
    genes = selected_genes,
    shap_values = shap_vec
  ))
}

#' Plot Final Scores Distribution
#'
#' Generates a histogram of gene scores and overlays key statistics.
#'
#' @param gene_scores A numeric vector of gene scores.
#' @param save_plot Logical. Whether to save the plot (default: FALSE).
#' @param file_path Character. File path for saving the plot (default: NULL).
#' @return A ggplot object of the distribution of gene scores.
plot_final_scores_distribution <- function(gene_scores, save_plot = FALSE,
                                           file_path ) {
  library(ggplot2)
  library(scales)
  # Φιλτράρισμα μη έγκυρων τιμών
  valid_scores <- !is.na(gene_scores) & gene_scores > 0
  filtered_scores <- gene_scores[valid_scores]
  # Check if gene_scores is valid
  if (length(gene_scores) == 0 || all(gene_scores == gene_scores[1])) {
    stop("Invalid input: gene_scores is empty or contains identical values.")
  }

  # Create a dataframe
  scores_data <- data.frame(
    Gene = names(gene_scores),
    Score = gene_scores
  )

  # Calculate statistics
  stats <- boxplot.stats(scores_data$Score)
  quartiles <- quantile(scores_data$Score, probs = c(0.25, 0.5, 0.75))
  mean_score <- mean(scores_data$Score) * 100
  scores_data <-na.omit(scores_data)
  # Create the plot
  p <- ggplot(scores_data, aes(x = Score)) +
    geom_histogram(binwidth = 0.01,
                   fill = "#99B4D1",
                   color = "white",
                   alpha = 1) +
    scale_y_log10()+
    geom_point(data = data.frame(x = stats$out),
               aes(x = x, y = -1),
               color = "red",
               size = 1,
               alpha = 0.6) +
    geom_segment(aes(x = stats$stats[1], xend = stats$stats[5],
                     y = -1, yend = -1),
                 color = "black",
                 size = 0.5) +
    geom_vline(xintercept = mean_score / 100,
               color = "red",
               linetype = "dashed",
               size = 0.5) +
    geom_vline(xintercept = quartiles,
               color = "darkblue",
               linetype = "dotted",
               size = 0.5,
               alpha = 0.5) +
    annotate("text",
             x = quartiles,
             y = c(60, 50, 40),
             label = sprintf("Q%d: %.1f%%", 1:3, quartiles * 100),
             color = "darkblue",
             size = 3,
             hjust = -0.1) +
    annotate("text",
             x = mean_score / 100,
             y = 65,
             label = sprintf("Mean: %.1f%%", mean_score),
             color = "red",
             size = 3.5,
             hjust = -0.1) +
    labs(title = "Distribution of Gene Scores",
         subtitle = "Analysis of Potential Biomarkers",
         x = "Gene Scores (Normalized)",
         y = "Frequency",
         caption = "Results from the Bayesian Gene Scoring System (BSGS)") +
    scale_x_continuous(labels = percent,
                       breaks = seq(0, 1, by = 0.25),
                       limits = c(0, 1)) +
    #scale_y_continuous(breaks = seq(0, 500, by = 20)) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9, color = "grey30"),
      panel.grid.major.y = element_line(color = "grey90", size = 0.3),
      panel.grid.major.x = element_line(color = "grey90", size = 0.3),
      panel.grid.minor = element_blank(),
      plot.caption = element_text(size = 8, hjust = 1, color = "grey50")
    )

  # Print distribution statistics in the console
  cat("\nDistribution Statistics:\n")
  cat(sprintf("1st Quartile (Q1): %.1f%%\n", quartiles[1] * 100))
  cat(sprintf("Median (Q2): %.1f%%\n", quartiles[2] * 100))
  cat(sprintf("3rd Quartile (Q3): %.1f%%\n", quartiles[3] * 100))
  cat(sprintf("Mean: %.2f%%\n", mean_score))

  # Show the plot
  print(p)

  # Save the plot if requested
  if (save_plot) {
    ggsave(filename = file_path, plot = p, dpi = 300, width = 10, height = 8, bg = "white")
  }
}
