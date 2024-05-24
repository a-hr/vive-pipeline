#!/usr/bin/env Rscript

library(tidyverse)
library(ggplot2)
library(patchwork)

# load data
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    stop("No input files provided")
}

ss5_scores <- read_delim(args[1], delim = "\t")
ss3_scores <- read_delim(args[2], delim = "\t")

# create a matrix with the scores for each k-mer
ss3.mat <- matrix(0, nrow = nrow(ss3_scores), ncol = max(ss3_scores$location2))
for (w in 1:nrow(ss3_scores)) {
    a <- ss3_scores$location[w]
    b <- ss3_scores$location2[w]
    ss3.mat[w, a:b] <- ss3_scores$ss3_scores[w]
}

# get the number of non-zero elements per column
ss3.nonzero <- colSums(ss3.mat != 0)

# average the scores for each position across all k-mers it appears in
ss3.vec <- colSums(ss3.mat) / ss3.nonzero

# add position index
ss3.df <- data.frame(position = 1:length(ss3.vec), score = ss3.vec)

# normalize the scores
# calculate the gap between the minimum value and 0 (if any)
gap <- max(0, -min(ss3.df$score))
# add the gap and subtract the minimum value
ss3.df$score <- ss3.df$score + gap
# divide by the maximum value
ss3.df$score <- ss3.df$score / max(ss3.df$score)


# plot 1D tile
p3 <- ggplot(ss3.df, aes(x = position, y = score)) +
    geom_bar(stat = "identity", color = "tomato3") +
    ggtitle("3' splice site scores") +
    xlab("Position in reference fasta") +
    ylab("Average score per position") +
    theme_minimal()


# calculate the gap between the minimum value and 0 (if any)
gap <- max(0, -min(ss5_scores$ss5_scores))
# add the gap and subtract the minimum value
ss5_scores$ss5_scores <- ss5_scores$ss5_scores + gap
# divide by the maximum value
ss5_scores$ss5_scores <- ss5_scores$ss5_scores / max(ss5_scores$ss5_scores)

p5 <- ggplot(ss5_scores, aes(x = location2, y = ss5_scores)) +
    geom_bar(stat = "identity", color = "navy") +
    ggtitle("5' splice site scores") +
    xlab("Position in reference fasta") +
    ylab("Score") +
    theme_minimal()

wrap_plots(p3, p5, nrow = 2, ncol = 1) %>%
    ggsave(filename = "splice_site_scores.png", width = 10, height = 6, dpi = 300)
