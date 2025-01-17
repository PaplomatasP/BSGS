#' Bayesian Gene Scoring System (BSGS)
#'
#' Implements a Bayesian-based feature selection algorithm for identifying significant genes.
#' This function performs iterative Bayesian updates to calculate gene scores and select
#' the most informative genes for classification tasks.
#'
#' @param data A data frame containing gene expression data.
#' @param label_col A string specifying the column name for the class labels in `data`.
#' @param iterations An integer specifying the number of iterations for the Bayesian updates (default = 10).
#' @param min_genes The minimum number of genes to sample in each iteration
#' (default = `round(floor(sqrt(number of genes)) / 2)`).
#' @param max_genes The maximum number of genes to sample in each iteration
#' (default = `5 * round(floor(sqrt(number of genes)))`).
#' @param alpha Weighting factor for the Bayesian update, balancing per-gene scores and subset scores (default = 0.3).
#' @param n_genes Number of top-ranked genes to return if no threshold is provided (default = 150).
#' @param plot_distribution Logical; if `TRUE`, generates a plot of the gene score distribution (default = FALSE).
#' @param save_plot Logical; if `TRUE`, saves the gene score distribution plot to a file (default = FALSE).
#' @param file_path A string specifying the file path for saving the gene score distribution plot
#' (default = "Results from the Bayesian Gene Scoring System (BSGS).png").
#'
#' @return A list with the following components:
#' \describe{
#'   \item{`final_genes`}{A character vector of the selected genes.}
#'   \item{`final_scores`}{A numeric vector of the scores corresponding to the selected genes.}
#'   \item{`posteriors`}{A numeric vector of the final posterior probabilities for all genes.}
#'   \item{`gene_scores`}{A numeric vector of the final gene scores for all genes.}
#' }
#'
#' @examples
#' \dontrun{
#' # Example usage
#' full_data <- as.data.frame(LUAD_MYE_MRS_GSE97168)
#' full_data$target <- as.factor(LUAD_MYE_MRS_target)
#' bsgs_results <- BSGS(
#'   data = full_data,
#'   label_col = "target",
#'   iterations = 1000,
#'   alpha = 0.3,
#'   min_genes = round(floor(sqrt(length(colnames(full_data[,-1])))) / 2),
#'   max_genes = 5 * round(floor(sqrt(length(colnames(full_data[,-1])))), 2),
#'   n_genes = 150,
#'   plot_distribution = TRUE,
#'   save_plot = TRUE,
#'   file_path = "./bsgs_results.png"
#' )
#' print(bsgs_results$final_genes)
#' }
#' @export
BSGS <- function(
    data,
    label_col,
    iterations = 1000,
    min_genes = round(floor(sqrt(length(colnames(data)[-1]))) / 2),
    max_genes = 5 * round(floor(sqrt(length(colnames(data)[-1])))),
    alpha = 0.3,
    n_genes = 150,
    plot_distribution = FALSE,
    save_plot = FALSE,
    file_path = "Results from the Bayesian Gene Scoring System (BSGS).png") {
  set.seed(1987)
  use_lasso_preselection <- TRUE

  # Optional LASSO Preselection
  shap_values <- numeric(0)
  if (use_lasso_preselection) {
    lasso_results <- lasso_preselection(data, label_col)
    initial_genes <- lasso_results$genes
    shap_values <- lasso_results$shap_values

    # Update dataset to include only selected genes
    labels <- data[[label_col]]
    data <- data[, colnames(data) %in% initial_genes, drop = FALSE]
    data[[label_col]] <- labels
  } else {
    # Use all genes if no LASSO preselection
    initial_genes <- setdiff(colnames(data), label_col)
    shap_values <- setNames(rep(0, length(initial_genes)), initial_genes)
  }

  gene_names <- setdiff(colnames(data), label_col)
  num_genes <- length(gene_names)

  # Initialize uniform priors
  priors <- rep(1 / num_genes, num_genes)

  # Initialize gene scores
  gene_scores <- setNames(rep(0, num_genes), gene_names)

  # Main loop
  for (iter in seq_len(iterations)) {
    num_genes_selected <- sample(min_genes:max_genes, 1)
    selected_indices <- sample_genes(priors, num_genes_selected)
    selected_genes <- gene_names[selected_indices]

    result <- compute_group_score(selected_genes, data, label_col, shap_values)
    set_likelihood <- result$group_score
    per_gene_scores <- result$per_gene_scores

    likelihoods <- rep(1e-6, num_genes)
    likelihoods[selected_indices] <- alpha * per_gene_scores + (1 - alpha) * set_likelihood

    posteriors <- update_posteriors(priors, likelihoods)

    for (i in seq_along(selected_indices)) {
      idx <- selected_indices[i]
      gene_scores[idx] <- max(gene_scores[idx], per_gene_scores[i])
    }

    cat(sprintf("Iteration %d: Chosen subset = %d genes, Group Score = %.4f\n",
                iter, num_genes_selected, set_likelihood))
  }

  sorted_indices <- order(gene_scores, decreasing = TRUE)
  final_genes <- gene_names[sorted_indices]
  final_scores <- gene_scores[sorted_indices]

  if (plot_distribution) {
    plot_final_scores_distribution(gene_scores, save_plot, file_path)

    threshold_input <- readline(prompt = "Provide the score threshold for gene selection: ")
    threshold <- as.numeric(threshold_input)

    if (is.na(threshold)) stop("Invalid input. Please enter a number.")

    final_scores <- final_scores * 100
    selected_genes <- final_genes[final_scores >= threshold]
    selected_scores <- final_scores[final_scores >= threshold]
  } else {
    selected_genes <- final_genes[1:n_genes]
    selected_scores <- final_scores[1:n_genes]
  }

  return(list(
    final_genes = selected_genes,
    final_scores = selected_scores,
    posteriors = posteriors,
    gene_scores = gene_scores
  ))
}
