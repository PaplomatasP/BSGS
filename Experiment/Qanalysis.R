# Load necessary libraries
library(ggplot2)
library(tidyr)

# Example data
iterations <- c(100, 1000, 5000, 10000, 20000, 50000)
q1 <- c(2.0, 25.2, 69.6, 89.5, 90.0, 90.0) # Q1 values
q2 <- c(6.7, 40.2, 70.3, 89.6, 90.0, 90.1) # Median (Q2) values
q3 <- c(15.2, 58.2, 89.4, 90.1, 90.3, 90.5) # Q3 values
mean <- c(12.4, 42.8, 76.2, 88.2, 90.5, 91.7) # Mean values

data <- data.frame(Iterations = iterations, Q1 = q1, Q2 = q2, Q3 = q3, Mean = mean)

# Convert data to long format for ggplot
data_long <- pivot_longer(data, cols = c("Q1", "Q2", "Q3", "Mean"), names_to = "Metric", values_to = "Value")

# Create the fancy plot
ggplot(data_long, aes(x = Iterations, y = Value, color = Metric, group = Metric)) +
  geom_line(size = 1.2) + # Bold lines for clarity
  geom_point(size = 3) + # Distinct points
  scale_x_log10(breaks = iterations, labels = scales::comma_format()) + # Log scale with clean labels
  scale_color_manual(values = c("Q1" = "#1b9e77", "Q2" = "#d95f02", "Q3" = "#7570b3", "Mean" = "#e7298a")) + # Custom colors
  labs(
    title = "Stabilization of Metrics Across Iterations",
    subtitle = "Q1, Q2, Q3, and Mean Values Over Iterations",
    x = "Iterations ",
    y = "Metric Values (%)",
    color = "Metric"
  ) +
  theme_minimal(base_size = 14) + # Clean minimal theme with larger font size
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid.major = element_line(size = 0.5, linetype = "dotted", color = "gray"),
    panel.grid.minor = element_blank() # Remove minor gridlines for clarity
  )


ggsave("Q1_Q2_Q3_mean_AcrossIter.png", dpi = 300, width = 10, height = 8, bg = "white")

library(UpSetR)


# Μετονομάζουμε τα στοιχεία της λίστας, προσθέτοντας "(n=...)" στο τέλος
names(list_of_subsets) <- paste0(
  names(list_of_subsets),
  " (n=", sapply(list_of_subsets, length), ")"
)



list_of_subsets <- list(
  subset100    = bsgs_results100$final_genes,
  subset1000   = bsgs_results1000$final_genes,
  subset5000   = bsgs_results5000$final_genes,
  subset10000  = bsgs_results10000$final_genes,
  subset20000  = bsgs_results20000$final_genes,
  subset50000  = bsgs_results50000$final_genes
)

upset_data <- fromList(list_of_subsets)

# Ορίζουμε τη σειρά των sets από το μικρότερο στο μεγαλύτερο.
my_order <- c("subset100", "subset1000", "subset5000", 
              "subset10000", "subset20000", "subset50000")

upset(
  upset_data,
  sets       = my_order,    # Τα σετ εμφανίζονται με αυτή την ακολουθία
  keep.order = TRUE,        # Διατηρεί την ακολουθία που δηλώσαμε
  order.by   = "freq",      # Πώς ταξινομούνται οι τομές (διατηρείται για την main bar)
  main.bar.color = "cornflowerblue",
  sets.bar.color = "darkorange"
)

ggsave("UpSetR.png", dpi = 300, width = 10, height = 8, bg = "white")
