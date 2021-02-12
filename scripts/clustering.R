# Copyright 2020 Dzuba, Krylov
#
# Licensed under the MIT License
# Permission is hereby granted, free of charge,
# to any person obtaining a copy of this software and associated documentation
# files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to permit
# persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.
#
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#
# This file contains the entire pipeline of cluster analysis of 2030 corporations financial data
# It is structured as follows:
#   1) Attaching libraries, raw data input;
#   2) Data cleaning, outlier removal;
#   3) Correlation analysis;
#   4) PCA analysis;
#   5) Data normalization techniques tests;
#   6) Finding the "most" suitable algorithm for data clustering;
#   7) K-means clustering and clusters visualization with parallel coordinates plots
#   8) Output tables


# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#
#                                                                              #
#                           Libraries, data input-----------------------------
#                                                                              #
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#


library("dplyr")
library("RColorBrewer")
library("Hmisc")
library("corrplot")
library("factoextra")
library("cluster")
library("ggplot2")
library("openxlsx")
library("parallelDist")
library("fastcluster")
library("summarytools")
library("egg")
library("tibble")
library("grid")
library("scales")
library("tidyr")
library("glue")

source("src/functions.R", encoding = "UTF-8")

PLOT_ALL_VIOLINPLOTS <- FALSE
FIND_BEST_CLUSTERING_ALGORITHM <- FALSE
USE_CALCULATED_BEST_K_MEANS_CENTERS <- TRUE
PATH_TO_CALCULATED_BEST_K_MEANS_CENTERS <- "data/centers/best_k_means_centers.RDS"
POSSIBLE_NUMBER_OF_CLUSTERS <- 2:15

VISUALIZE_CLUSTERINGS <- 6:8
BEST_CLUSTERING <- 7

PLOT_COUNTER <- 0

cluster_names <- list(
  "k6" = c(
    "P+ L-",
    "F+",
    "R+ F-",
    "L+ G+ T+",
    "L+ G+",
    "RPL"
  ),
  "k7" = c(
    "P+ L-",
    "RPL",
    "L+ G+ T+",
    "R- P- L+",
    "F+",
    "R+ F-",
    "L+ G+"
  ),
  "k8" = c(
    "F+",
    "L+ G+",
    "R+ F-",
    "P+ L-",
    "L++ G+ T++",
    "T+",
    "RPL",
    "R- P- L+"
  )
)


df <- df_raw <- read.csv("data/raw/corp_data.csv")
var_index <- c("P", "R", "L", "G", "F", "T")


# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#
#                                                                              #
#                       Outlier analysis-------------------------------
#                                                                              #
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#


out1 <- which(is.na(df$L) | is.infinite(df$L) | is.na(df$R) | is.infinite(df$R) |
  is.na(df$G) | is.infinite(df$G) | is.na(df$F) | is.infinite(df$F) |
  is.na(df$P) | is.infinite(df$P) | is.na(df$T) | is.infinite(df$T))


df[out1, var_index] <- NA


# Outliers were removed based on visual inspection of box_violin plots.
# This process was performed in several iterations,
# as smaller outliers are not visible if there are strong ones.
# Boundries iterations:
# (0)
plot.box_violin(na.omit(df), bnd = NULL, arrange_p = T)

boundries_names <- c(
  "P_low",
  "P_high",
  "R_low",
  "R_high",
  "L_low",
  "L_high",
  "G_low",
  "G_high",
  "F_low",
  "F_high",
  "T_low",
  "T_high"
)

if (PLOT_ALL_VIOLINPLOTS) {
  # (1)
  bnd1 <- c(-50, 50, -200, 100, -50, 10, -50, 100, -1, 50, -20, 10) %>%
    `names<-`(boundries_names)
  plot.box_violin(na.omit(df), bnd = bnd1, arrange_p = T)

  # (2)
  bnd2 <- c(-4, 3, -10, 50, -2, 4, -5, 5, -0.5, 5, -5, 4) %>%
    `names<-`(boundries_names)
  plot.box_violin(na.omit(df), bnd = bnd2, arrange_p = T)

  # (3)
  bnd3 <- c(-1.5, 1.5, -3, 2, -0.25, 2.75, -3, 3.5, -0.25, 1.2, -1.25, 2) %>%
    `names<-`(boundries_names)
  plot.box_violin(na.omit(df), bnd = bnd3, arrange_p = T)

  # (4)
  bnd4 <- c(-0.7, 0.9, -1.8, 1.2, -0.15, 1.75, -1.05, 1.95, -0.25, 1.05, -0.55, 1.35) %>%
    `names<-`(boundries_names)
  plot.box_violin(na.omit(df), bnd = bnd4, arrange_p = T)

  # (5)
  bnd5 <- c(-0.6, 0.725, -1.35, 1, -0.14, 1.75, -0.75, 1.45, -0.25, 1.05, -0.25, 1.35) %>%
    `names<-`(boundries_names)
  plot.box_violin(na.omit(df), bnd = bnd5, arrange_p = T)
}

# (6)
bnd6 <- c(
  "P_low" = -0.45,
  "P_high" = 0.6,
  "R_low" = -0.65,
  "R_high" = 0.94,
  "L_low" = 0,
  "L_high" = 1.17,
  "G_low" = 0,
  "G_high" = 1.45,
  "F_low" = -0.25,
  "F_high" = 1.05,
  "T_low" = -0.25,
  "T_high" = 1.013
) # 6 sigma-like distribution

PLOT_COUNTER <- PLOT_COUNTER + 1

svg(glue("figs/{PLOT_COUNTER} boxplot.svg"), height = 3, family = "serif")
plot.box_violin(na.omit(df),
  bnd = bnd6,
  arrange_p = T,
  base_p_size = 8
)
dev.off()

out2 <- which(df$P < bnd6["P_low"] | df$P > bnd6["P_high"] | df$R < bnd6["R_low"] | df$R > bnd6["R_high"] |
  df$L < bnd6["L_low"] | df$L > bnd6["L_high"] | df$G < bnd6["G_low"] | df$G > bnd6["G_high"] |
  df$F < bnd6["F_low"] | df$F > bnd6["F_high"] | df$T < bnd6["T_low"] | df$T > bnd6["T_high"])

df[out2, var_index] <- NA

which_not_NA <- which(!is.na(df$P))


# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#
#                                                                              #
#                         Correlation analysis------------------------
#                                                                              #
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#


cor_matrix <- Hmisc::rcorr(as.matrix(df[which_not_NA, var_index]), type = "pearson")

PLOT_COUNTER <- PLOT_COUNTER + 1

svg(glue("figs/{PLOT_COUNTER} corrplot.svg"), family = "serif")
corrplot::corrplot(
  cor_matrix$r,
  method = "ellipse", type = "lower", addCoef.col = "black",
  tl.col = "black", tl.srt = 0, p.m = cor_matrix$P
)
dev.off()


# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#
#                                                                              #
#                         Principal components analysis------------------------
#                                                                              #
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#


pca_res <- prcomp(df[which_not_NA, var_index], scale. = TRUE)
fviz_screeplot(pca_res)
plot(pca_res$x[, 1], pca_res$x[, 2])


# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#
#                                                                              #
#                         Data normalization analysis------------------------
#                                                                              #
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#


test_z_score <- df[which_not_NA, c("R", "T")] %>% scale_data("z-score")
test_min_max <- df[which_not_NA, c("R", "T")] %>% scale_data("min-max")

plot_z_score <- ggplot(test_z_score, aes(x = R, y = T)) +
  geom_point(alpha = 0.05, size = 1) +
  theme_bw(base_size = 11) +
  labs(x = "R", y = "T")

plot_min_max <- ggplot(test_min_max, aes(x = R, y = T)) +
  geom_point(alpha = 0.05, size = 1) +
  theme_bw(base_size = 11) +
  labs(x = "R", y = "T")

PLOT_COUNTER <- PLOT_COUNTER + 1

svg(glue("figs/{PLOT_COUNTER} z-score vs min-max.svg"), 7, 3.5, family = "serif")
egg::ggarrange(plot_z_score,
  plot_min_max,
  ncol = 3,
  labels = letters[1:2],
  label.args = list(
    gp = gpar(cex = 1)
  )
)
dev.off()


# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#
#                                                                              #
#                           Finding the best clustering algorithm---------------
#                                                                              #
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#


df[, var_index] <- scale_data(df[, var_index], method = "z-score")

if (!USE_CALCULATED_BEST_K_MEANS_CENTERS) {
  kmeans_centers_silh_scores <- calculate_kmeans_centers_silh_scores(
    df = df[which_not_NA, ],
    var_index = var_index,
    k = POSSIBLE_NUMBER_OF_CLUSTERS,
    iterations = 100,
    save_best_centers_to_RDS = TRUE,
    path_to_RDS = "data/centers/best_k_means_centers.RDS"
  )
} else {
  if (!file.exists(PATH_TO_CALCULATED_BEST_K_MEANS_CENTERS)) {
    stop("No data with calculated best k-means centers has been found with
         provided path variable 'PATH_TO_CALCULATED_BEST_K_MEANS_CENTERS'")
  } else {
    kmeans_centers_silh_scores <- readRDS(PATH_TO_CALCULATED_BEST_K_MEANS_CENTERS)
  }
}


if (FIND_BEST_CLUSTERING_ALGORITHM) {
  find_good_algos <- list()

  best_kmeans_centers <- lapply(
    kmeans_centers_silh_scores[["best_centers"]],
    function(x) x[, var_index]
  )

  dist_mat <- parallelDist::parallelDist(as.matrix(df[which_not_NA, var_index]))

  find_good_algos[["kmeans"]] <- cluster_data(
    df = df[which_not_NA, var_index],
    k = POSSIBLE_NUMBER_OF_CLUSTERS,
    init_clus_centers = best_kmeans_centers[glue("k {POSSIBLE_NUMBER_OF_CLUSTERS}")],
    cluster_method = "kmeans",
    hclust_method = NULL,
    dist_mat = dist_mat,
    hclust_tree = NULL,
    internal_metrics = TRUE,
    verbose = TRUE
  )

  # README
  # pam is slow - run if needed

  find_good_algos[["pam"]] <- cluster_data(
    df[which_not_NA, var_index],
    k = POSSIBLE_NUMBER_OF_CLUSTERS,
    cluster_method = "pam",
    hclust_method = "ward.D2",
    dist_mat = dist_mat,
    verbose = TRUE,
    internal_metrics = TRUE,
    pam_pamonce = 3
  )

  find_good_algos[["ward.D2"]] <- cluster_data(
    df[which_not_NA, var_index],
    k = POSSIBLE_NUMBER_OF_CLUSTERS,
    cluster_method = "hclust",
    hclust_method = "ward.D2",
    dist_mat = dist_mat,
    verbose = TRUE,
    internal_metrics = TRUE
  )

  compare_algos <- data.frame(
    "k" = POSSIBLE_NUMBER_OF_CLUSTERS,
    "kmeans" = unlist(find_good_algos[["kmeans"]][["avg.silhouette.width"]]),
    "k-medoids_pam" = unlist(find_good_algos[["pam"]][["avg.silhouette.width"]]),
    "agnes ward.D2" = unlist(find_good_algos[["ward.D2"]][["avg.silhouette.width"]])
  ) %>%
    `names<-`(c("k", "k-means", "k-medoids (PAM)", "agnes ward.D2")) %>%
    pivot_longer(cols = 2:ncol(.))

  PLOT_COUNTER <- PLOT_COUNTER + 1

  svg(glue("figs/{PLOT_COUNTER} algorithms average silhouette.svg"), height = 2.5, family = "serif")
  ggplot(compare_algos, aes(x = as.factor(k), y = value, shape = name, group = name)) +
    geom_point(size = 1) +
    geom_line(size = 0.1) +
    guides(shape = guide_legend(
      override.aes = list(alpha = 1, size = 1),
      title = "Clustering algorithms\n(Euclidean Distance)",
      labels = unique(compare_algos[["name"]])
    )) +
    theme_bw(base_size = 8) +
    theme(
      plot.margin = margin(0, 0, 0, 0, "cm"),
      legend.position = c(0.87, 0.75),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.background = element_rect(fill = "transparent")
    ) +
    ylab("Average silhouette width") +
    xlab("Number of clusters k")
  dev.off()
} else {
  PLOT_COUNTER <- PLOT_COUNTER + 1
}


# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#
#                                                                              #
#                           k-means clustering-----------------
#                                                                              #
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#


df[, var_index] <- scale_data(df[, var_index], method = "z-score")

best_kmeans_centers <- lapply(
  kmeans_centers_silh_scores[["best_centers"]],
  function(x) x[, var_index]
)

dist_mat <- parallelDist::parallelDist(as.matrix(df[which_not_NA, var_index]), "euclidean")

# silhoutte calculations costs a lot of RAM, set internal_metrics to FALSE if needed
kmeans_res <- cluster_data(
  df = df[which_not_NA, var_index],
  k = POSSIBLE_NUMBER_OF_CLUSTERS,
  init_clus_centers = best_kmeans_centers,
  cluster_method = "kmeans",
  hclust_method = NULL,
  dist_mat = dist_mat,
  hclust_tree = NULL,
  internal_metrics = TRUE,
  verbose = TRUE
)

print("Average silhouette width for different clusterings:")
print(unlist(kmeans_res[["avg.silhouette.width"]]))

clusterings <- list()
df_average_cluster <- list()
silhouette_res <- list()

for (i in VISUALIZE_CLUSTERINGS) {
  ci <- glue("k{i}") # cluster index

  temp_cluster <- kmeans_res[["clustering"]][[ci]]

  if (ci %in% names(cluster_names)) {
    for (j in seq_len(length(unique(temp_cluster)))) {
      temp_cluster[temp_cluster == j] <- cluster_names[[ci]][j]
    }
  }

  clusterings[[ci]] <- temp_cluster
}

rm(temp_cluster)

for (i in VISUALIZE_CLUSTERINGS) {
  ci <- glue("k{i}") # cluster index

  PLOT_COUNTER <- PLOT_COUNTER + 1

  svg(glue("figs/{PLOT_COUNTER} k{i} clusters in parallel coord.svg"), width = 10, height = 6, family = "serif")
  print(
    plot.par_coord(df[which_not_NA, var_index],
      clusterings[[ci]],
      alpha_trasp = 0.05,
      scale_y_step = 1,
      base_p_size = 8,
      line_size = 0.08
    ) +
      geom_hline(
        yintercept = 0,
        color = "grey30",
        linetype = "twodash"
      )
  )
  dev.off()
}

for (i in VISUALIZE_CLUSTERINGS) {
  ci <- glue("k{i}") # cluster index

  df_average_cluster[[ci]] <- aggregate.data.frame(
    df[which_not_NA, var_index],
    list("cluster" = clusterings[[ci]]),
    mean
  )

  PLOT_COUNTER <- PLOT_COUNTER + 1

  svg(glue("figs/{PLOT_COUNTER} k{i} mean cluster values.svg"), height = 2.5, family = "serif")
  print(
    plot.par_coord(df_average_cluster[[ci]][, var_index],
      df_average_cluster[[ci]][, "cluster"],
      clust_plot_scheme = "shape",
      alpha_trasp = 1,
      scale_y_step = 1,
      base_p_size = 8,
      line_size = 0.08, legend_pos = c(0.9, 0.2)
    ) +
      scale_y_continuous(
        name = "Z-score",
        limits = c(min(df_average_cluster[[ci]][, var_index]), max(df_average_cluster[[ci]][, var_index])),
        breaks = seq(round(min(df_average_cluster[[ci]][, var_index])),
          round(max(df_average_cluster[[ci]][, var_index])),
          by = 0.5
        ),
        minor_breaks = NULL,
        expand = c(0, 0.1)
      ) +
      theme(plot.title = element_blank())
  )
  dev.off()
}

for (i in VISUALIZE_CLUSTERINGS) {
  ci <- glue("k{i}") # cluster index

  silhouette_res[[ci]] <- as.data.frame(kmeans_res[["silhouette"]][[ci]][, 1:3])
  silhouette_res[[ci]][, "cluster"] <- clusterings[[ci]]

  PLOT_COUNTER <- PLOT_COUNTER + 1

  svg(glue("figs/{PLOT_COUNTER} k{i} silhoutte plot k{i}.svg"), height = 3, family = "serif")
  print(
    plot.silhouette(silhouette_res[[ci]], base_p_size = 8)
  )
  dev.off()
}


# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#
#                                                                              #
#                                 Output---------------
#                                                                              #
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#


bci <- glue("k{BEST_CLUSTERING}") # best cluster index

best_clustering <- clusterings[[bci]]

# return unnormalized data
df[, var_index] <- df_raw[, var_index]

df[, "cluster"] <- 0
df[, "silhouette"] <- NA

df[, "cluster"][which_not_NA] <- best_clustering
df[, "silhouette"][which_not_NA] <- silhouette_res[[bci]][, "sil_width"]

df[, "cluster"][out1] <- "InfNaNNA"
df[, "cluster"][out2] <- "Boxplot outlier"

for (i in unique(silhouette_res[[bci]][, "neighbor"])) {
  silhouette_res[[bci]][, "neighbor"][silhouette_res[[bci]][, "neighbor"] == i] <- cluster_names[[bci]][i]
}

df[, "silhouette_neighbor"] <- NA
df[which_not_NA, "silhouette_neighbor"] <- silhouette_res[[bci]][, "neighbor"]

write.table(df,
  "output/data_clustering_results.tsv",
  quote = TRUE,
  sep = "\t",
  na = "NA",
  row.names = FALSE,
  col.names = TRUE
)

write.xlsx(df,
  "output/data_clustering_results.xlsx",
  keepNA = TRUE
)

summary_clusters <- by(
  data = df[which_not_NA, var_index],
  INDICES = list(cluster = df[which_not_NA, "cluster"]),
  FUN = descr,
  stats = c("min", "q1", "med", "mean", "sd", "q3", "max")
) %>%
  lapply(function(x) cbind.data.frame(x, rownames(x))) %>%
  do.call(rbind.data.frame, args = .) %>%
  cbind.data.frame(rep(levels(factor(df[which_not_NA, "cluster"])), each = 7))

cluster_box_violins_plots <- by(
  data = df[which_not_NA, var_index],
  INDICES = list(cluster = df[which_not_NA, "cluster"]),
  FUN = function(x) {
    plist <- list()
    for (i in 1:ncol(x)) {
      plist[[i]] <- ggplot(x, aes_string(y = names(x)[i], x = rep(1, nrow(x)))) +
        geom_violin(trim = FALSE, fill = "gray", colour = "grey") +
        geom_boxplot(width = 0.1, outlier.size = 0.5, outlier.alpha = 0.15) +
        theme_bw(base_size = 11) +
        xlab(names(x)[i]) +
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(1, 5, 5, 5, "pt")
        ) +
        scale_y_continuous(
          breaks = scales::pretty_breaks(),
          limits = c(
            min(df[which_not_NA, var_index[i]]) - 0.05,
            max(df[which_not_NA, var_index[i]])
          )
        ) +
        theme(
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.title.y = element_blank()
        )
    }
    egg::ggarrange(
      plots = plist,
      ncol = 3,
      nrow = 2,
      # labels = letters[1:ncol(x)],
      # label.args = list(
      #   gp = grid::gpar(cex = 1)
      # ),
      draw = FALSE
    )
  }
)

for (i in names(cluster_box_violins_plots)) {
  PLOT_COUNTER <- PLOT_COUNTER + 1
    svg(glue("figs/{PLOT_COUNTER} {i}, boxplot.svg"), family = "serif")
  print(cluster_box_violins_plots[[i]])
  dev.off()
} 
