library(tidyverse)
library(ggrepel)
library(mybnlearn)

rm(list=ls())
source("R/utils.R")

# Asia true network
true = model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]")

# Load the object as an .rds file with the timestamped filename
# Specify the directory path
directory_path <- "simulations/mcmc"

# List all .rds files in the directory
rds_files <- list.files(directory_path, pattern = "\\.rds$", full.names = TRUE)
rds_files
filenames <- grep("08-20", rds_files, value = TRUE)
filename <- rds_files[1]
for (filename in filenames)
{
  result <- readRDS(file = filename)
  n_chains <- max(unique(result$chain))
  N <- 5000
  sample_size <- as.numeric(sub(".*_(\\d+)_.*", "\\1", filename))
  
  # Remove starting DAGs
  result <- result %>% filter(phase != "Start") %>% filter(phase != "Burn")
  
  # Visual inspection of convergence of MCMC chains.
  result %>%
    ggplot(aes(x = step, y = score, linetype = as.factor(chain), col = as.factor(iss))) +
    geom_line(alpha = 0.5, linewidth = 1) + 
    scale_color_manual(values = tab10)
  
  # Visualizing the most frequents DAGs.
  plot_tbl <- result %>%
    group_by(iss, hash) %>%
    summarise(
      dag = list(first(dag)),
      prior = first(prior),
      iss = first(iss),
      rel_iss = first(iss) / sample_size,
      count = n(),
      freq = count / (n_chains * N),
      score = first(score)) %>%
    mutate(
      rank = dense_rank(desc(score)),
      normalized_score = case_when(
        min(score) == max(score) ~ 1,
        TRUE ~ (score - min(score)) / (max(score) - min(score))),
      shdem = purrr::map_dbl(dag, ~ shd_edge_marks(true, .x)) / 2) %>%
    filter(freq > quantile(freq, 0.9)) %>%
    ungroup() %>%
    mutate(
      iss_level = match(iss, unique(result$iss)) - 1,
      offset = (iss_level / max(iss_level)) * 0.8,
      shdem = shdem + offset) %>%
    select(!c(iss_level, offset))
  
  
  # Find the DAG with the highest score within each iss.
  highest_freq <- plot_tbl %>%
    group_by(iss) %>%
    filter(freq == max(freq)) %>%
    filter(shdem == min(shdem)) %>%
    mutate(
      lbl1 = paste0("ISS ", iss),
      lbl2 = paste0("ISS ", iss, ", # ", ceiling(rank))) %>%
    ungroup()
  
  # Plot score vs SHD Edge Marks
  score_plot <- plot_tbl %>%
    ggplot(aes(y = shdem, x = score, col = as.factor(iss), shape = prior, size = freq)) +
    geom_point(alpha = 0.5) +
    geom_label_repel(data = highest_freq, aes(y = shdem, x = score, label = lbl1), 
                     color = "black", size = 5, label.size = NA, nudge_x = -4, nudge_y = 0.5) +
    scale_y_continuous(limits = c(0, 8), 
                       breaks = unique(floor(plot_tbl$shdem))) +
    scale_color_manual(values = tab10,
                       labels = c(prior = "Prior", freq = "Frequency"),
                       name = "ISS") +
    scale_shape_manual(values = c(16,17),
                       name = "Prior") +
    scale_size(range = c(2, 8),
               name = "Frequency") +
    labs(
      x = "Score",
      y = "SHD EM"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.position = "top"
    ) +
    guides(
      shape = guide_legend(nrow = 2, override.aes = list(size = 8)),
      size = guide_legend(nrow = 2),
      col = guide_legend(nrow = 2, override.aes = list(size = 8))
    )
  score_plot
  
  
  # Plot score vs SHD Edge Marks
  norm_score_plot <- plot_tbl %>%
    ggplot(aes(y = shdem, x = normalized_score, col = as.factor(iss), shape = prior, size = freq)) +
    geom_point(alpha = 0.5) +
    geom_label_repel(data = highest_freq, aes(y = shdem, x = normalized_score, label = lbl1), 
                     color = "black", size = 5, label.size = NA, nudge_x = -0.2, nudge_y = 0.5) +
    scale_y_continuous(limits = c(0, 8), 
                       breaks = unique(floor(plot_tbl$shdem))) +
    scale_color_manual(values = tab10,
                       labels = c(prior = "Prior", freq = "Frequency"),
                       name = "ISS") +
    scale_shape_manual(values = c(16,17),
                       name = "Prior") +
    scale_size(range = c(2, 8),
               name = "Frequency") +
    labs(
      x = "Normalized Score",
      y = "SHD EM"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.position = "top"
    ) +
    guides(
      shape = guide_legend(nrow = 2, override.aes = list(size = 8)),
      size = guide_legend(nrow = 2),
      col = guide_legend(nrow = 2, override.aes = list(size = 8))
    )
  norm_score_plot
  
  # Plot SHD vs rank
  nudge_x <- 50
  rank_plot <- plot_tbl %>%
    ggplot(aes(x = rank, y = shdem, col = as.factor(iss), shape = prior, size = freq)) +
    geom_point(alpha = 0.5) +
    geom_label_repel(data = highest_freq, aes(x = rank, y = shdem, label = lbl2), 
                     color = "black", size = 5, label.size = NA, 
                     nudge_x = nudge_x,
                     nudge_y = 1) +
    xlim(0, max(highest_freq$rank) + nudge_x + 50) +
    scale_y_continuous(limits = c(0, 8),
                       breaks = seq(0, 8)) +
    scale_color_manual(values = tab10,
                       labels = c(prior = "Prior", freq = "Frequency"),
                       name = "ISS") +
    scale_shape_manual(values = c(16,17),
                       name = "Prior") +
    scale_size(range = c(2, 8),
               name = "Frequency") +
    labs(
      x = "Score rank",
      y = "SHD EM"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.position = "top"
    ) +
    guides(
      shape = guide_legend(nrow = 2, override.aes = list(size = 8)),
      size = guide_legend(nrow = 2),
      col = guide_legend(nrow = 2, override.aes = list(size = 8))
    )
  rank_plot
  
  # Create directory if it doesn't exist
  output_directory <- paste0("figures/mcmc/", n_chains, "-", N)
  if (!file.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
  }
  
  # Construct the filename
  filename <- paste0(output_directory, "/score-", n_chains, "-", sample_size, ".pdf")
  ggsave(filename, plot = score_plot, width = 8, height = 6)
  filename <- paste0(output_directory, "/norm_score-", n_chains, "-", sample_size, ".pdf")
  ggsave(filename, plot = norm_score_plot, width = 8, height = 6)
  filename <- paste0(output_directory, "/rank-", n_chains, "-", sample_size, ".pdf")
  ggsave(filename, plot = rank_plot, width = 8, height = 6)  
  
  rm(result)
}
