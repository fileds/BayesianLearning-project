library(tidyverse)
library(ggrepel)
library(mybnlearn)
library(knitr)
library(kableExtra)

rm(list=ls())
source("R/utils.R")

group = "rank"
#group = "shdem"

# Asia true network
true = model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]")

# Load the object as an .rds file with the timestamped filename
# Specify the directory path
directory_path <- "simulations/iss"

# List all .rds files in the directory
rds_files <- list.files(directory_path, pattern = "\\.rds$", full.names = TRUE)
rds_files
filenames <- grep("iss_30_5000", rds_files, value = TRUE)
filenames <- grep("iss_30_5000_1000_", rds_files, value = TRUE)
for (filename in filenames)
{
  result <- readRDS(file = filename)
  n_experiments <- max(unique(result$experiment))
  N <- result$count[[1]] / result$freq[[1]]
  sample_size <- as.numeric(sub(".*_(\\d+)_.*", "\\1", filename))
  
  # Correlation
  correlation <- result %>%
    group_by(iss, experiment) %>%
    summarise(
      score_corr = cor(shdem, score),
      rank_corr = cor(shdem, rank)) %>%
    summarise(
      score_corr = round(mean(score_corr), 2),
      rank_corr = round(mean(rank_corr), 2)
    ) %>%
    mutate(sample_size = sample_size) %>%
    select(sample_size, everything())
  correlation
  
  formatted_table <- kable(correlation, format = "latex", booktabs = TRUE) %>%
    kable_styling()
  
  # Tibble for visualization
  if (group == "rank")
  {
    # Find the minmax rank across all ISS
    minmax_rank <- 
      result %>%
      group_by(experiment) %>%
      summarise(rank = max(rank)) %>%
      pull(rank) %>% 
      min()
    
    plot_tbl <- result %>%
      filter(rank <= minmax_rank) %>%
      group_by(iss, rank) %>%
      summarise(
        prior = first(prior),
        score = mean(score),
        shdem = mean(shdem)) %>%
      mutate(
        normalized_score = (score - min(score)) / (max(score) - min(score))) %>%
      ungroup()
  }
  else if (group == "shdem")
  {
    plot_tbl <- result %>%
      group_by(iss, shdem) %>%
      summarise(
        prior = first(prior),
        score = mean(score),
        rank = mean(rank),
        N = n()) %>%
      mutate(
        normalized_score = (score - min(score)) / (max(score) - min(score))) %>%
      ungroup()
  }
  else
  {
    print("No group selected")  
  }
  
  # Plot score vs SHD Edge Marks
  score_plot <- plot_tbl %>%
    ggplot(aes(x = shdem, y = score, col = as.factor(iss), shape = prior)) +
    geom_point(alpha = 0.5, size = 5) +
    scale_x_continuous(breaks = floor(unique(plot_tbl$shdem))) +
    scale_color_manual(values = tab10,
                       labels = c(prior = "Prior", freq = "Frequency"),
                       name = "ISS") +
    scale_shape_manual(values = c(16,17),
                       name = "Prior") +
    labs(
      y = "Score",
      x = "SHD EM",
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.position = "top"
    )
  score_plot
  
  # Plot normalized score vs SHD Edge Marks
  norm_score_plot <- plot_tbl %>%
    ggplot(aes(x = shdem, y = normalized_score, col = as.factor(iss), shape = prior)) +
    geom_point(alpha = 0.5, size = 5) +
    scale_x_continuous(breaks = floor(unique(plot_tbl$shdem))) +
    scale_color_manual(values = tab10,
                       name = "ISS") +
    scale_shape_manual(values = c(16,17),
                       name = "Prior") +
    labs(
      y = "Normalized score",
      x = "SHD EM",
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.position = "top"
    )
  norm_score_plot
  
  # Plot SHD vs rank
  rank_plot <- plot_tbl %>%
    ggplot(aes(x = shdem, y = rank, col = as.factor(iss), shape = prior)) +
    geom_point(alpha = 0.5, size = 3) +
    scale_x_continuous(breaks = floor(unique(plot_tbl$shdem))) +
    scale_color_manual(values = tab10,
                       name = "ISS") +
    scale_shape_manual(values = c(16,17),
                       name = "Prior") +
    labs(
      y = "Score rank",
      x = "SHD EM"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.position = "top"
    )
  rank_plot
  
  # Create directory if it doesn't exist
  #output_directory <- paste0("figures/iss/", n_experiments, "-", N, "/", group)
  output_directory <- paste0("figures/iss/", n_experiments, "-", N)
  if (!file.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
  }
  
  # Construct the filename
  filename <- paste0(output_directory, "/score-", sample_size, ".pdf")
  ggsave(filename, plot = score_plot, width = 8, height = 6)
  filename <- paste0(output_directory, "/norm_score-", sample_size, ".pdf")
  ggsave(filename, plot = norm_score_plot, width = 8, height = 6)
  filename <- paste0(output_directory, "/rank-", sample_size, ".pdf")
  ggsave(filename, plot = rank_plot, width = 8, height = 6)  
  
  output_directory <- paste0(output_directory, "/tables")
  if (!file.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
  }
  filename <- paste0(output_directory, "/", sample_size, ".tex")
  writeLines(formatted_table, filename)
}
