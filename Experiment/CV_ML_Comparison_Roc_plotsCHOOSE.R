classification_comparison_CV <- function(full_data, 
                                         selected_genes, 
                                         label_col,
                                         models = c("GBM", "RandomForest", "kNN", "NaiveBayes"), 
                                         folds = 10,  
                                         NZV = TRUE,
                                         cores = NULL,
                                         run_full = TRUE,  
                                         run_subset = TRUE 
) {
  # -- Βιβλιοθήκες --------------------------------------------------------- #
  library(doParallel)
  library(caret)
  library(randomForest)
  library(pROC)
  library(e1071)
  library(gbm)
  library(class)  # για kNN (knn function)
  
  # -- Έλεγχος ότι υπάρχει η στήλη-ετικέτα (labels) ------------------------- #
  if (!label_col %in% colnames(full_data)) {
    stop("Label column not found in data")
  }
 
  # -- Oρισμός cores για παράλληλη εκτέλεση -------------------------------- #
  if (is.null(cores)) {
    cores <- max(1, parallel::detectCores() - 1)
  }
  cl <- makeCluster(cores)
  on.exit(stopCluster(cl))  # τερματισμός cluster σε περίπτωση σφάλματος
  registerDoParallel(cl)
  
  set.seed(2111987)
  colnames(full_data) <- make.names(colnames(full_data))
  # -- Αφαίρεση zero/near-zero variance features αν ζητηθεί ---------------- #
  if (NZV) {
    nzv <- nearZeroVar(full_data[, !colnames(full_data) %in% label_col])
    if (length(nzv) > 0) {
      warning(
        "Removed zero/near-zero variance predictors: ", 
        paste(colnames(full_data)[nzv], collapse = ", ")
      )
      full_data <- full_data[, -nzv]
    }
  }
  
  # -- Δημιουργία folds για cross-validation -------------------------------- #
  fold_indices <- createFolds(full_data[[label_col]], k = folds, list = TRUE)
  
  # -- Προετοιμασία μεταβλητής "subsets_list" ------------------------------ #
  if (is.list(selected_genes)) {
    # Ο χρήστης έδωσε ήδη λίστα από subsets
    subsets_list <- selected_genes
  } else {
    # Έδωσε μόνο ένα subset (vector)
    subsets_list <- list(single_subset = selected_genes)
  }
  
  # -- Συνάρτηση υπολογισμού metrics --------------------------------------- #
  calculate_metrics <- function(true_labels, pred_probs, pred_labels) {
    # true_labels: factor (π.χ. c("diseases","healthy"))
    # pred_probs : numeric, οι πιθανότητες για το "θετικό" class
    # pred_labels: factor, οι τελικές προβλέψεις ως κλάση
    
    positive_class <- levels(true_labels)[2]  # υποθέτουμε ότι το [2] είναι το "θετικό"
    true_numeric <- as.numeric(true_labels == positive_class)  # 0/1
    pred_numeric <- as.numeric(pred_labels == positive_class)  # 0/1
    
    cm <- confusionMatrix(
      factor(pred_numeric, levels = c(0, 1)),
      factor(true_numeric, levels = c(0, 1))
    )
    
    list(
      AUC         = roc(true_numeric, pred_probs)$auc,
      Accuracy    = cm$overall["Accuracy"],
      Sensitivity = cm$byClass["Sensitivity"],
      Specificity = cm$byClass["Specificity"],
      F1          = cm$byClass["F1"]
    )
  }
  
  # -- Συνάρτηση εκπαίδευσης/πρόβλεψης για κάθε μοντέλο -------------------- #
  train_and_predict <- function(model_type, train_data, test_data, label_col) {
    # Επιστρέφει λίστα με:
    #  - probs: numeric vector (πιθανότητες ή σκορ)
    #  - labels: factor (προβλεπόμενη κλάση)
    tryCatch({
      switch(model_type,
             "GBM" = {
               model <- gbm(
                 as.formula(paste(label_col, "~ .")), 
                 data = train_data, 
                 distribution = "bernoulli", 
                 n.trees = 100,
                 interaction.depth = 3,
                 verbose = FALSE
               )
               probs <- predict(model, newdata = test_data, type = "response", n.trees = 100)
               labels <- ifelse(probs > 0.5,
                                levels(train_data[[label_col]])[2],
                                levels(train_data[[label_col]])[1])
               list(
                 probs = probs,
                 labels = factor(labels, levels = levels(train_data[[label_col]]))
               )
             },
             "RandomForest" = {
               model <- randomForest(as.formula(paste(label_col, "~ .")), data = train_data)
               probs <- predict(model, newdata = test_data, type = "prob")[, 2]
               labels <- predict(model, newdata = test_data)
               list(probs = probs, labels = labels)
             },
             "kNN" = {
               train_x <- train_data[, !colnames(train_data) %in% label_col, drop = FALSE]
               train_y <- train_data[[label_col]]
               test_x  <- test_data[,  !colnames(test_data) %in% label_col,  drop = FALSE]
               
               pred_labels <- knn(
                 train = train_x,
                 test  = test_x,
                 cl    = train_y,
                 k = 5
               )
               positive_class <- levels(train_y)[2]
               probs <- as.numeric(pred_labels == positive_class)
               
               list(
                 probs = probs,
                 labels = factor(pred_labels, levels = levels(train_y))
               )
             },
             "NaiveBayes" = {
               model <- naiveBayes(as.formula(paste(label_col, "~ .")), data = train_data)
               probs <- predict(model, newdata = test_data, type = "raw")[, 2]
               labels <- predict(model, newdata = test_data)
               list(probs = probs, labels = labels)
             }
      )
    }, error = function(e) {
      warning(paste("Error in model", model_type, ":", e$message))
      NULL
    })
  }
  
  # Αρχικοποιούμε δομή αποτελεσμάτων
  # Θα σώζουμε και metrics και predictions
  results <- list(
    Full = list(),
    Subsets = list()
  )
  
  # ------------------------------------------------------------------------ #
  # 1) FULL DATASET (αν run_full=TRUE)
  # ------------------------------------------------------------------------ #
  if (run_full) {
    cat("\nRunning cross-validation on the FULL dataset...\n")
    for (model in models) {
      cat(sprintf("  - Evaluating %s model...\n", model))
      
      model_folds_metrics <- list()  # για τα metrics
      model_folds_preds   <- list()  # για τις προβλέψεις
      
      for (i in seq_along(fold_indices)) {
        test_indices <- fold_indices[[i]]
        train_data <- full_data[-test_indices, ]
        test_data  <- full_data[test_indices, ]
        
        # Προαιρετικό scaling
        preProc <- preProcess(
          train_data[, !colnames(train_data) %in% label_col, drop=FALSE],
          method = c("center", "scale")
        )
        train_scaled <- predict(preProc, train_data)
        test_scaled  <- predict(preProc, test_data)
        
        # Εκπαίδευση & πρόβλεψη
        preds <- train_and_predict(model, train_scaled, test_scaled, label_col)
        
        if (!is.null(preds)) {
          # Yπολογίζουμε metrics
          fold_metrics <- calculate_metrics(
            true_labels = test_data[[label_col]],
            pred_probs  = preds$probs,
            pred_labels = preds$labels
          )
          model_folds_metrics[[i]] <- fold_metrics
          
          # Αποθηκεύουμε τις προβλέψεις (για MCC, ROC κ.λπ. αργότερα)
          # Θέλουμε να θυμόμαστε τα "αληθινά" labels & τα numeric scores
          model_folds_preds[[i]] <- list(
            labels       = test_data[[label_col]], # factor
            pred_scores  = preds$probs,            # numeric (probabilities)
            pred_classes = preds$labels            # factor
          )
        }
      }
      
      # Υπολογίζουμε average metrics
      if (length(model_folds_metrics) > 0) {
        metric_names <- names(model_folds_metrics[[1]])
        avg_metrics <- sapply(metric_names, function(mn) {
          mean(sapply(model_folds_metrics, `[[`, mn), na.rm = TRUE)
        }, simplify = FALSE)
        avg_metrics <- setNames(unlist(avg_metrics), metric_names)
        
        # Αποθηκεύουμε *και* τις μέσες μετρικές, *και* τις fold-level προβλέψεις
        results$Full[[model]] <- list(
          metrics     = avg_metrics,
          predictions = model_folds_preds
        )
      } else {
        # Αν δεν βρήκαμε τίποτα για αυτό το μοντέλο
        results$Full[[model]] <- NULL
      }
    }
  }
  
  # ------------------------------------------------------------------------ #
  # 2) SUBSETS (αν run_subset=TRUE)
  # ------------------------------------------------------------------------ #
  if (run_subset) {
    cat("\nRunning cross-validation on all SUBSETS...\n")
    
    for (subset_name in names(subsets_list)) {
      cat(sprintf("\n  Subset: %s\n", subset_name))
      
      gene_subset <- subsets_list[[subset_name]]
      
      # Αρχικοποιούμε μία λίστα για τα αποτελέσματα αυτού του subset
      subset_results <- list()
      
      # Για κάθε μοντέλο
      for (model in models) {
        cat(sprintf("    - Evaluating %s model...\n", model))
        
        model_folds_metrics <- list()
        model_folds_preds   <- list()
        
        for (i in seq_along(fold_indices)) {
          test_indices <- fold_indices[[i]]
          train_data <- full_data[-test_indices, ]
          test_data  <- full_data[test_indices, ]
          
          # Κρατάμε μόνο τα columns που ανήκουν στο subset + την ετικέτα
          valid_genes <- intersect(gene_subset, colnames(train_data))
          if (length(valid_genes) == 0) {
            warning(sprintf(
              "No valid genes found in subset '%s' for fold %d. Skipping this subset.",
              subset_name, i
            ))
            next
          }
          
          train_subset <- train_data[, c(valid_genes, label_col), drop=FALSE]
          test_subset  <- test_data[,  c(valid_genes, label_col), drop=FALSE]
          
          # Optional scaling
          preProc <- preProcess(
            train_subset[, !colnames(train_subset) %in% label_col, drop=FALSE],
            method = c("center", "scale")
          )
          train_scaled <- predict(preProc, train_subset)
          test_scaled  <- predict(preProc, test_subset)
          
          # Εκπαίδευση & πρόβλεψη
          preds <- train_and_predict(model, train_scaled, test_scaled, label_col)
          if (!is.null(preds)) {
            fold_metrics <- calculate_metrics(
              true_labels = test_subset[[label_col]],
              pred_probs  = preds$probs,
              pred_labels = preds$labels
            )
            model_folds_metrics[[i]] <- fold_metrics
            
            model_folds_preds[[i]] <- list(
              labels       = test_subset[[label_col]],
              pred_scores  = preds$probs,
              pred_classes = preds$labels
            )
          }
        }
        
        # Μέσος όρος
        if (length(model_folds_metrics) > 0) {
          metric_names <- names(model_folds_metrics[[1]])
          avg_metrics <- sapply(metric_names, function(mn) {
            mean(sapply(model_folds_metrics, `[[`, mn), na.rm = TRUE)
          }, simplify = FALSE)
          avg_metrics <- setNames(unlist(avg_metrics), metric_names)
          
          subset_results[[model]] <- list(
            metrics     = avg_metrics,
            predictions = model_folds_preds
          )
        } else {
          subset_results[[model]] <- NULL
        }
      }
      
      # Αποθήκευση
      results$Subsets[[subset_name]] <- subset_results
    }
  }
  
  return(results)
}



full_data <- as.data.frame(h5ad_df)
colnames(full_data) <- make.names(colnames(full_data))

full_data$labels <- as.factor(labels)

# Subsets
list_of_subsets <- list(
  subset100    = bsgs_results100$final_genes,
  subset1000   = bsgs_results1000$final_genes,
  subset5000   = bsgs_results5000$final_genes,
  subset10000  = bsgs_results10000$final_genes,
  subset20000  = bsgs_results20000$final_genes,
  subset50000  = bsgs_results50000$final_genes
)

results_all <- classification_comparison_CV(
  full_data      = full_data,
  selected_genes = list_of_subsets,
  label_col      = "labels",
  models         = c( "kNN", "NaiveBayes","RandomForest"),
  folds          = 5,
  NZV            = FALSE,
  cores          = 28,
  run_full       = TRUE,   # ή FALSE ανάλογα τι θέλεις
  run_subset     = TRUE
)


str(results_all)

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(yardstick)

# Assuming results_all_ML is already defined in your environment

###############################
## 1) Βιβλιοθήκες και ορισμοί
###############################

# Υποθέτουμε ότι στο περιβάλλον υπάρχει το αντικείμενο results_all_ML
# (ή π.χ. results_all), που έχει τη δομή:
# results_all_ML$Full$kNN$metrics, results_all_ML$Full$kNN$predictions, κ.λπ.

###############################
## 2) Υπολογισμός MCC από τα fold-level predictions
###############################
# Παραδείγματα ονομάτων μοντέλων και datasets
models <- c("kNN", "NaiveBayes", "RandomForest")

# Εδώ, απλώς θα ενώσουμε όλα τα folds στο "Full" σε ένα data frame,
# και το ίδιο για κάθε subset, ώστε να υπολογίζουμε έναν MCC συνολικά.
mcc_results <- list()

for (model in models) {
  
  # 2.1) FULL dataset
  if (!is.null(results_all_ML$Full[[model]]$predictions)) {
    all_labels <- c()
    all_scores <- c()
    
    # Συνένωση όλων των folds
    for (i in seq_along(results_all_ML$Full[[model]]$predictions)) {
      fold_pred <- results_all_ML$Full[[model]]$predictions[[i]]
      all_labels <- c(all_labels, as.character(fold_pred$labels))
      all_scores <- c(all_scores, fold_pred$pred_scores)
    }
    
    # Δημιουργούμε predicted_classes με threshold=0.5
    predicted_classes <- ifelse(all_scores > 0.5, "tumor", "normal")
    
    df_temp <- data.frame(
      truth    = factor(all_labels, levels=c("normal","tumor")),
      estimate = factor(predicted_classes, levels=c("normal","tumor"))
    )
    
    # Υπολογισμός MCC (μέσω yardstick)
    mcc_val <- mcc(df_temp, truth, estimate)
    mcc_results[[paste(model, "Full", sep="_")]] <- mcc_val$.estimate
  } else {
    mcc_results[[paste(model, "Full", sep="_")]] <- NA
  }
  
  # 2.2) SUBSETS
  if (!is.null(results_all_ML$Subsets)) {
    # Για κάθε subset_name
    for (subset_name in names(results_all_ML$Subsets)) {
      if (!is.null(results_all_ML$Subsets[[subset_name]][[model]]$predictions)) {
        all_labels <- c()
        all_scores <- c()
        for (i in seq_along(results_all_ML$Subsets[[subset_name]][[model]]$predictions)) {
          fold_pred <- results_all_ML$Subsets[[subset_name]][[model]]$predictions[[i]]
          all_labels <- c(all_labels, as.character(fold_pred$labels))
          all_scores <- c(all_scores, fold_pred$pred_scores)
        }
        predicted_classes <- ifelse(all_scores > 0.5, "tumor", "normal")
        
        df_temp <- data.frame(
          truth    = factor(all_labels, levels=c("normal","tumor")),
          estimate = factor(predicted_classes, levels=c("normal","tumor"))
        )
        mcc_val <- mcc(df_temp, truth, estimate)
        
        mcc_results[[paste(model, subset_name, sep="_")]] <- mcc_val$.estimate
      } else {
        # Αν δεν έχει predictions για αυτό το subset
        mcc_results[[paste(model, subset_name, sep="_")]] <- NA
      }
    }
  }
}

# Μετατρέπουμε τη λίστα mcc_results σε data frame
mcc_df <- do.call(rbind, lapply(names(mcc_results), function(name) {
  parts <- strsplit(name, "_")[[1]]   # π.χ. c("kNN","Full") ή c("RF","subset100")
  this_model   <- parts[1]
  this_dataset <- parts[2:length(parts)] %>% paste(collapse="_")  # χειρισμός αν έχει underscores
  
  data.frame(
    Model   = this_model,
    Dataset = this_dataset,
    Metric  = "MCC",
    Value   = mcc_results[[name]],
    stringsAsFactors = FALSE
  )
}))

###############################
## 3) Δημιουργία data frame με “κλασικές” μετρικές (AUC, Accuracy, κ.λπ.)
###############################
make_metrics_dataframe <- function(results_list) {
  rows <- list()
  
  # 3.1) Full
  if (!is.null(results_list$Full)) {
    for (model_name in names(results_list$Full)) {
      model_entry <- results_list$Full[[model_name]]
      if (!is.null(model_entry) && !is.null(model_entry$metrics)) {
        metrics_vec <- model_entry$metrics
        df_temp <- data.frame(
          Model   = model_name,
          Dataset = "Full",
          Metric  = names(metrics_vec),
          Value   = as.numeric(metrics_vec),
          stringsAsFactors = FALSE
        )
        rows[[paste0("Full_", model_name)]] <- df_temp
      }
    }
  }
  
  # 3.2) Subsets
  if (!is.null(results_list$Subsets)) {
    for (subset_name in names(results_list$Subsets)) {
      subset_obj <- results_list$Subsets[[subset_name]]
      for (model_name in names(subset_obj)) {
        model_entry <- subset_obj[[model_name]]
        if (!is.null(model_entry) && !is.null(model_entry$metrics)) {
          metrics_vec <- model_entry$metrics
          df_temp <- data.frame(
            Model   = model_name,
            Dataset = subset_name,
            Metric  = names(metrics_vec),
            Value   = as.numeric(metrics_vec),
            stringsAsFactors = FALSE
          )
          rows[[paste0(subset_name, "_", model_name)]] <- df_temp
        }
      }
    }
  }
  
  final_df <- do.call(rbind, rows)
  rownames(final_df) <- NULL
  return(final_df)
}

metrics_df <- make_metrics_dataframe(results_all_ML)


###############################
## 4) Συνένωση MCC με υπάρχουσες μετρικές
###############################
# Εφόσον στο metrics_df έχουμε π.χ. (Model, Dataset, Metric, Value)
# και στο mcc_df έχουμε κάτι αντίστοιχο, μπορούμε να κάνουμε rbind.

all_metrics_df <- rbind(metrics_df, mcc_df)

###############################
## 5) Δημιουργία γραφήματος σύγκρισης
###############################
library(ggplot2)
library(dplyr)

ggplot(all_metrics_df, aes(x = Model, y = Value, fill = Dataset)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  # Προσθήκη ετικετών μετρικής πάνω από τις μπάρες
  geom_text(
    aes(label = round(Value, 2), y = Value),
    position = position_dodge(width = 0.8),
    vjust = -0.3,
    size = 3
  ) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 2) +
  labs(
    title = "Performance Metrics by Model and Dataset",
    x = "Model",
    y = "Metric Value",
    fill = "Dataset"
  ) +
  # Αν θέλεις χρωματική παλέτα π.χ. Set2 (θα υπερισχύσει των χειροκίνητων χρωμάτων)
  scale_fill_brewer(palette = "Set2") +
  # Εναλλακτικά, αν επιμένεις σε custom χειροκίνητα χρώματα:
  # scale_fill_manual(values = c("Full" = "blue", "subset100" = "red", 
  #                              "subset1000" = "green", "subset10000" = "purple",
  #                              "subset20000" = "orange","subset50000" = "pink")) +
  theme_bw(base_size = 14) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x     = element_text(angle = 45, hjust = 1),
    strip.text      = element_text(face = "bold", size = 14),
    legend.position = "top"
  )


############# HEAT MAP ################

library(dplyr)
library(ggplot2)

# Suppose you have all_metrics_df with columns: Model, Dataset, Metric, Value
# We'll rename "Dataset" -> "Subset" and "Metric" -> "variable"
library(dplyr)

# Modify the 'Subset' column to rename 'Full' and replace subset names
heatmap_data <- heatmap_data %>%
  mutate(
    Subset = case_when(
      Subset == "Full" ~ "Complete Dataset",  # Rename 'Full' to 'Complete Dataset'
      grepl("subset", Subset) ~ gsub("subset", "BSGS", Subset),  # Replace 'subset' with 'BSGS'
      TRUE ~ Subset  # For any other cases, keep as is
    )
  )

# Verify the changes
head(heatmap_data)


# Now we create the heatmap:
ggplot(heatmap_data, aes(x = Subset, y = variable, fill = Value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", Value)), color = "black", size = 4) +
  scale_fill_viridis_c(option = "plasma", name = "Value") +
  facet_wrap(~ Model, ncol = 1) +
  labs(
    title = "Performance Evaluation of BSGS and Traditional Feature Selection Methods",
    x = "Subset",
    y = "Metric"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "none",
    panel.grid.major = element_line(size = 0.5, color = "grey85"),
    panel.grid.minor = element_blank()
  )
ggsave("Heatmap_subsetsANDfull_performance_comparison.png", dpi = 300, width = 10, height = 8, bg = "white")
