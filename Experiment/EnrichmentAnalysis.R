
########### Enrichemnt Analysis 

enrichment_analysis <- function(genes,dbs = "KEGG_2021_Human",multiPlot=TRUE ) {
  websiteLive=TRUE
  library(enrichR)
  library(ggplot2)
  library("gridExtra")
  # dbs <- c(
  #   "KEGG_2021_Human",
  #   "WikiPathway_2021_Human",
  #   "BioPlanet_2019",
  #   "BioCarta_2016",
  #   "Reactome_2016",
  #   "MSigDB_Hallmark_2020",
  #   "GO_Biological_Process_2021",
  #   "GO_Molecular_Function_2021",
  #   "GO_Cellular_Component_2021",
  #   "MGI_Mammalian_Phenotype_Level_4_2021",
  #   "Human_Phenotype_Ontology",
  #   "Jensen_DISEASES",
  #   "DisGeNET",
  #   "DSigDB",
  #   "DrugMatrix",
  #   "OMIM_Disease",
  #   "HDSigDB_Human_2021",
  #   "COVID-19_Related_Gene_Sets_2021"
  # ),
  
  enriched <- enrichr(genes, dbs)
  
  # Check if any databases returned no results
  empty_dbs <- which(sapply(enriched, nrow) == 0)
  if (length(empty_dbs) > 0) {
    for (i in empty_dbs) {
      enriched[[i]] <- data.frame(matrix(rbinom(10 * 10, 1, 0), ncol = 9))
      colnames(enriched[[i]]) <- c(
        "Term",
        "Overlap",
        "P.value",
        "Adjusted.P.value",
        "Old.P.value",
        "Old.Adjusted.P.value",
        "Odds.Ratio",
        "Combined.Score",
        "Genes"
      )
    }
  }
  
  # Create barplots for each database
  BBP <- list()
  Table_from_Enrichment <- list()
  for (i in seq_along(enriched)) {
    allPlotData <- as.data.frame(enriched[[i]][1:10,])
    allPlotData <- allPlotData[order(allPlotData$Combined.Score, decreasing = FALSE), drop = FALSE, ]
    Table_from_Enrichment[[i]] <- allPlotData
    BBP[[i]] <- ggplot(data = allPlotData, aes(
      x = reorder(Term, +Combined.Score),
      y = Combined.Score
    )) +
      ggtitle(dbs[i]) +
      geom_bar(stat = "identity",
               position = "dodge",
               aes(fill = Combined.Score)) + coord_flip() + theme_gray() +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none"
      ) +
      geom_text(
        aes(label = Term),
        hjust = -0.002,  # Ορίζει την αριστερή ευθυγράμμιση
        size = 3.5,
        position = position_stack(vjust = 0),
        inherit.aes = TRUE
      )+
      scale_fill_distiller(name = "Value",
                           palette = "Reds",
                           direction = 1)
  }
  if (multiPlot){
    BBP <-do.call("grid.arrange", c(BBP, ncol = 1))    # Apply do.call & grid.arrange
    
  }
  
  return(BBP)
}



enrichment_analysis(bsgs_results20000$final_genes,dbs = c("GO_Biological_Process_2021",
                                                     "GO_Molecular_Function_2021") )
