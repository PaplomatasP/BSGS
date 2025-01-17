# Load necessary libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("STRINGdb")

library(STRINGdb)

# Function to create a PPI network analysis for the isolated genes
PPInetwork <- function(Genes, organismus = "Human", Score_Threshold_PPI = 400, GenesN = NULL) {

  
  # Determine the species based on the organismus input
  if (organismus == "Human") {
    PPI_species <- 9606
  } else if (organismus == "Mouse") {
    PPI_species <- 10090
  } else {
    stop("Unsupported organismus. Please choose 'Human' or 'Mouse'.")
  }
  
  # Create a STRINGdb object
  string_db <- STRINGdb$new(
    version = "11.5",
    species = PPI_species,
    score_threshold = Score_Threshold_PPI,
    input_directory = ""
  )
  
  # Create a data frame with the genes and their scores
  dfGenes <- data.frame(Genes = Genes)
  
  # Map the genes to STRING IDs
  Mapped <- string_db$map(dfGenes, "Genes", removeUnmappedRows = TRUE)
  hits <- Mapped$STRING_id
  
  # Plot the PPI network
  string_db$plot_network(hits)
}
final_genes <- bsgs_results20000$final_genes

PPInetwork(final_genes,Score_Threshold_PPI = 200)
