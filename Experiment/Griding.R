library(caret)
library(e1071)


# Ορισμός του grid search
param_grid <- expand.grid(
  gamma = c(1, 2, 3),
  alpha = c(0.3, 0.4, 0.5),
  min_genes = c(50, 100, 150),
  max_genes = c(300, 400, 500),
  n_genes = c(50, 100, 150)
)

# Ορισμός των παραμέτρων για 10-fold CV
set.seed(1987)
train_control <- trainControl(method = "cv", number = 10, savePredictions = TRUE)

best_params <- NULL
best_score <- -Inf

for (i in 1:nrow(param_grid)) {
  gamma <- param_grid$gamma[i]
  alpha <- param_grid$alpha[i]
  min_genes <- param_grid$min_genes[i]
  max_genes <- param_grid$max_genes[i]
  n_genes <- param_grid$n_genes[i]
  
  # Εκτέλεση του BSGS αλγορίθμου
  bsgs_results <- BSGS(
    data = full_data,
    alpha = alpha,
    label_col = "labels",
    iterations = 1000,
    min_genes = min_genes,
    max_genes = max_genes,
    n_genes = n_genes
  )
  
  # Απομόνωση των επιλεγμένων γονιδίων και του στόχου
  selected_genes <- bsgs_results$final_genes
  data_selected <- full_data[, selected_genes]
  data_selected$labels <- full_data$labels
  
  # Δημιουργία μοντέλου με 10-fold CV
  rf_cv_model <- train(
    labels ~ .,
    data = data_selected,
    method = "rf",
    trControl = train_control,
    importance = TRUE,
    ntree = 100
  )
  
  # Υπολογισμός των μέτρων απόδοσης (Accuracy, F1, Precision, Recall)
  conf_matrix <- confusionMatrix(rf_cv_model$pred$pred, rf_cv_model$pred$obs)
  accuracy <- conf_matrix$overall['Accuracy']
  f1_score <- mean(conf_matrix$byClass['F1'])
  precision <- mean(conf_matrix$byClass['Precision'])
  recall <- mean(conf_matrix$byClass['Recall'])
  
  # Εμφάνιση των αποτελεσμάτων
  cat("Gamma:", gamma, "Alpha:", alpha, "Min Genes:", min_genes, "Max Genes:", max_genes, "N Genes:", n_genes, "\n")
  cat("Accuracy:", accuracy, "\n")
  cat("F1 Score (mean):", f1_score, "\n")
  cat("Precision (mean):", precision, "\n")
  cat("Recall (mean):", recall, "\n")
  
  # Αποθήκευση των καλύτερων παραμέτρων και σκορ
  if (f1_score > best_score) {
    best_score <- f1_score
    best_params <- c(gamma, alpha, min_genes, max_genes, n_genes)
  }
}
saveRDS(best_params,"best_params.rds")
# Εμφάνιση των καλύτερων παραμέτρων
cat("Best Parameters:\n")
cat("Gamma:", best_params[1], "\n")
cat("Alpha:", best_params[2], "\n")
cat("Min Genes:", best_params[3], "\n")
cat("Max Genes:", best_params[4], "\n")
cat("N Genes:", best_params[5], "\n")
cat("Best F1 Score:", best_score, "\n")



##################### GRID metrics parameter ###############################


library(caret)
library(e1071)

# Δημιουργία grid για τα βάρη
# Θα χρησιμοποιήσουμε βήματα του 0.1 και θα διασφαλίσουμε ότι αθροίζουν στο 1
weights_grid <- expand.grid(
  snr_weight = seq(0.1, 0.5, by = 0.1),
  cor_weight = seq(0.1, 0.5, by = 0.1),
  var_weight = seq(0.1, 0.5, by = 0.1),
  shap_weight = seq(0.1, 0.5, by = 0.1)
)

# Κρατάμε μόνο τους συνδυασμούς που αθροίζουν στο 1
weights_grid$sum <- weights_grid$snr_weight + weights_grid$cor_weight + 
  weights_grid$var_weight + weights_grid$shap_weight
weights_grid <- subset(weights_grid, sum == 1)
weights_grid$sum <- NULL

# Ρύθμιση για 10-fold CV
set.seed(1987)
train_control <- trainControl(method = "cv", number = 10, savePredictions = TRUE)

# Αποθήκευση καλύτερων παραμέτρων
best_weights <- NULL
best_score <- -Inf

# Grid Search
for (i in 1:nrow(weights_grid)) {
  # Τροποποίηση της compute_group_score function με τα τρέχοντα βάρη
  compute_group_score_modified <- function(selected_genes, data, label_col, shap_values) {
    result <- compute_group_score(selected_genes, data, label_col, shap_values)
    
    # Επαναϋπολογισμός των subset_scores με τα τρέχοντα βάρη
    subset_scores <- weights_grid$snr_weight[i] * result$norm_snr + 
      weights_grid$cor_weight[i] * result$norm_cor + 
      weights_grid$var_weight[i] * result$norm_var + 
      weights_grid$shap_weight[i] * result$norm_shap
    
    return(list(
      group_score = mean(subset_scores, na.rm = TRUE),
      per_gene_scores = subset_scores
    ))
  }
  
  # Εκτέλεση BSGS με τα τροποποιημένα βάρη
 bsgs_results <- BSGS(
  data = full_data,
  alpha = 0.3,
  label_col = "labels",
  iterations = 1000,
  min_genes =  round(floor(sqrt(length(colnames(full_data[,-1])))) / 2) ,
  max_genes = 5 * round(floor(sqrt(length(colnames(full_data[,-1])))) , 2) , 
  n_genes = 150
)
  
  # Επιλογή γονιδίων και δημιουργία dataset
  selected_genes <- bsgs_results$final_genes
  data_selected <- full_data[, selected_genes]
  data_selected$labels <- full_data$labels
  
  # Εκπαίδευση και αξιολόγηση με Random Forest
  rf_cv_model <- train(
    labels ~ .,
    data = data_selected,
    method = "rf",
    trControl = train_control,
    importance = TRUE,
    ntree = 100
  )
  
  # Υπολογισμός μετρικών
  conf_matrix <- confusionMatrix(rf_cv_model$pred$pred, rf_cv_model$pred$obs)
  f1_score <- mean(conf_matrix$byClass['F1'])
  
  # Εκτύπωση τρέχοντων αποτελεσμάτων
  cat("Weights combination", i, "of", nrow(weights_grid), ":\n")
  cat("SNR:", weights_grid$snr_weight[i], 
      "Cor:", weights_grid$cor_weight[i],
      "Var:", weights_grid$var_weight[i],
      "SHAP:", weights_grid$shap_weight[i], "\n")
  cat("F1 Score:", f1_score, "\n\n")
  
  # Ενημέρωση καλύτερων παραμέτρων
  if (f1_score > best_score) {
    best_score <- f1_score
    best_weights <- c(
      snr = weights_grid$snr_weight[i],
      cor = weights_grid$cor_weight[i],
      var = weights_grid$var_weight[i],
      shap = weights_grid$shap_weight[i]
    )
  }
}

# Εκτύπωση τελικών αποτελεσμάτων
cat("Best Weights Found:\n")
cat("SNR weight:", best_weights["snr"], "\n")
cat("Correlation weight:", best_weights["cor"], "\n")
cat("Variance weight:", best_weights["var"], "\n")
cat("SHAP weight:", best_weights["shap"], "\n")
cat("Best F1 Score:", best_score, "\n")

# Αποθήκευση των καλύτερων βαρών
saveRDS(best_weights, "best_weights.rds")