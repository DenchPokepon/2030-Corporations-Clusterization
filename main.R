# This file contains the entire pipeline of cluster analysis of 2030 corporations financial data
# It is structured as follows:
#   1) Attaching libraries, raw data input;
#   2) Data cleaning, outlier removal;
#   3) Correlation analysis;
#   4) PCA analysis;
#   5) Data normalization techniques tests;
#   6) Finding the most suitable algorithm for data clustering;
#   7) K-means clustering and clusters visualization with parallel coordinates plots


# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#
#                                                                              #
#                           Libraries, data input-----------------------------
#                                                                              #
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#

library(dplyr)
library(RColorBrewer)
library(Hmisc)
library(corrplot)
library(factoextra)
library(cluster)
library(ggplot2)
library(openxlsx)
library(parallelDist)
library(fastcluster)
library(summarytools)
library(egg)
library(tibble)
library(grid)
library(scales)
library(tidyr)
library(crayon)

source("src/functions.R", encoding = "UTF-8")

colours <- c("cyan", RColorBrewer::brewer.pal(8, name = "Set1"))

df <- read.csv("data/corp data.csv")
var_index <- which(names(df) %in% c("R", "P", "L", "G", "F", "T"))


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
plot.box_violin(df, bnd = NULL, arrange_p = T)

# (1)
bnd1 <- c(-200, 100, -50, 50, -50, 10, -50, 100, -1, 50, -20, 10)
plot.box_violin(df, bnd = bnd1, arrange_p = T)

# (2)
bnd2 <- c(-10, 50, -4, 3, -2, 4, -5, 5, -0.5, 5, -5, 4)
plot.box_violin(df, bnd = bnd2, arrange_p = T)

# (3)
bnd3 <- c(-3, 2, -1.5, 1.5, -0.25, 2.75, -3, 3.5, -0.25, 1.2, -1.25, 2)
plot.box_violin(df, bnd = bnd3, arrange_p = T)

# (4)
bnd4 <- c(-1.8, 1.2, -0.7, 0.9, -0.15, 1.75, -1.05, 1.95, -0.25, 1.05, -0.55, 1.35)
plot.box_violin(df, bnd = bnd4, arrange_p = T)

# (5)
bnd5 <- c(-1.35, 1, -0.6, 0.725, -0.14, 1.75, -0.75, 1.45, -0.25, 1.05, -0.25, 1.35)
plot.box_violin(df, bnd = bnd5, arrange_p = T)

# (6)
bnd6 <- c(
  "R_low" = -0.65,
  "R_high" = 0.94,
  "P_low" = -0.45,
  "P_high" = 0.6,
  "L_low" = 0,
  "L_high" = 1.17,
  "G_low" = 0,
  "G_high" = 1.45,
  "F_low" = -0.25,
  "F_high" = 1.05,
  "T_low" = -0.25,
  "T_high" = 1.013
) # 6 sigma-like distribution

svg("figs/Boxplot.svg", height = 3.5, family = "serif")
plot.box_violin(df,
                bnd = bnd6,
                arrange_p = T,
                base_p_size = 11
)
dev.off()

out2 <- which(df$R < bnd6[1] | df$R > bnd6[2] | df$P < bnd6[3] | df$P > bnd6[4] |
                     df$L < bnd6[5] | df$L > bnd6[6] | df$G < bnd6[7] | df$G > bnd6[8] |
                     df$F < bnd6[9] | df$F > bnd6[10] | df$T < bnd6[11] | df$T > bnd6[12])

df[out2, var_index] <- NA

which_not_NA <- which(!is.na(df$R))


# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#
#                                                                              #
#                         Correlation analysis------------------------
#                                                                              #
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#


cor_matrix <- Hmisc::rcorr(as.matrix(df[which_not_NA, var_index]), type = "pearson")

svg("figs/Corrplot.svg", family = "serif")
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


pca_res <- prcomp(df[which_not_NA, var_index])
fviz_screeplot(pca_res)


# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#
#                                                                              #
#                         Data normalization analysis------------------------
#                                                                              #
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#


test_z_score <- df[, c("R", "T")] %>% scale_data("z-score")
test_min_max <- df[, c("R", "T")] %>% scale_data("min-max")

p1 <- ggplot(test_z_score, aes(x = R, y = T)) +
  geom_point(alpha = 0.05, size = 1) +
  theme_bw(base_size = 11) +
  labs(x = "R", y = "T")

p2 <- ggplot(test_min_max, aes(x = R, y = T)) +
  geom_point(alpha = 0.05, size = 1) +
  theme_bw(base_size = 11) +
  labs(x = "R", y = "T")

svg("figs/z-score vs min-max.svg", 7, 3.5, family = "serif")
egg::ggarrange(p1,
               p2,
               ncol = 3,
               labels = letters[1:2],
               label.args = list(
                 gp = gpar(cex = 1)
               )
)
dev.off()


# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#
#                                                                              #
#                           Finding best clustering algorithm----------------
#                                                                              #
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#
# 
# README
# high memory usage, run with 16 gb RAM at least or use SWAP
# README
#
#
#
#
# df[, var_index] <- scale_data(df[, var_index],
#                                method = "z-score",
#                                optimum_values = c(0.1, 0.1, 0.5, 0.1, 0.2, 0.1)
# )
# 
# # k-means best clusters centers detection
# clusters_res <- list()
# dist_mat <- parallelDist::parallelDist(as.matrix(df[which_not_NA, var_index]))
# 
# for (k in 2:15) {
#   guess_seeds <- round(rnorm(100, mean = 20, sd = 10) * 10000, 0)
#   
#   kmeans_res <- list()
#   
#   for (i in guess_seeds) {
#     set.seed(i)
#     centers <- df[sample(which_not_NA, k), var_index]
#     kmeans_res[[paste(i)]] <- kmeans(df[which_not_NA, var_index], centers)
#   }
# 
#   
#   silhs <- lapply(kmeans_res, function(cluster_res) {
#     silhouette(cluster_res$cluster, dist_mat)
#   })
#   
#   clusters_res[[paste("k =", k)]] <- lapply(silhs, function(x) mean(x[, 3])) %>%
#     do.call(c, args = .)
#   gc()
# }
# lapply(clusters_res, function(x) max(x))
# lapply(clusters_res, function(x) which(x == max(x)))
#
# seeds for best (max silhouette and parcoord) kmeans centers sampling (z-score, euclidean distance)
# seeds for x64 R version 4.0.2 (2020-06-22) os mingw32

best_kmeans_seeds <- c(
  k2 = 93038,
  k3 = 113663,
  k4 = 114088,
  k5 = 101955,
  k6 = 101766,
  k7 = 103747,
  k8 = 89527,
  k9 = 114078,
  k10 = 190583,
  k11 = 116827,
  k12 = 104093,
  k13 = 102839,
  k14 = 88195,
  k15 = 94217
)

find_good_algor <- list()

centers_list <- list()

for (i in 1:length(best_kmeans_seeds)) {
  set.seed(best_kmeans_seeds[i])
  centers_list[[paste(best_kmeans_seeds[i])]] <-
    df[sample(which_not_NA, i + 1), var_index]
}

dist_mat <- parallelDist::parallelDist(as.matrix(df[which_not_NA, var_index]))

find_good_algor[["kmeans"]] <- cluster.data(
  df = df[which_not_NA, var_index],
  k = 2:15,
  init_clus_centers = centers_list,
  cluster_method = "kmeans",
  hclust_method = NULL,
  dist_mat = dist_mat,
  hclust_tree = NULL,
  internal_metrics = TRUE,
  verbose = TRUE
)
# pam is slow algorithm - run if needed
find_good_algor[["pam"]] <- cluster.data(
  df[which_not_NA, var_index],
                                         k = 2:15,
                                         cluster_method = "pam",
                                         hclust_method = "ward.D2",
                                         dist_mat = dist_mat,
                                         verbose = TRUE,
                                         internal_metrics = TRUE
)
find_good_algor[["ward.D2"]] <- cluster.data(
  df[which_not_NA, var_index],
                                             k = 2:15,
                                             cluster_method = "hclust",
                                             hclust_method = "ward.D2",
                                             dist_mat = dist_mat,
                                             verbose = TRUE,
                                             internal_metrics = TRUE
)
find_good_algor[["complete"]] <- cluster.data(
  df[which_not_NA, var_index],
                                              k = 2:15,
                                              cluster_method = "hclust",
                                              hclust_method = "complete",
                                              dist_mat = dist_mat,
                                              verbose = TRUE,
                                              internal_metrics = TRUE
)

compare_algors <- data.frame(
  "k" = 2:15,
  "kmeans" = find_good_algor$kmeans$Avg.silhouette.width[2:15],
  "k-medoids_pam" = find_good_algor$pam$Avg.silhouette.width[2:15],
  "agnes ward.D2" = find_good_algor$ward.D2$Avg.silhouette.width[2:15],
  "agnes complete" = find_good_algor$complete$Avg.silhouette.width[2:15]
) %>%
  pivot_longer(cols = 2:ncol(.))

algors <- c(
  "k-means",
  "k-medoids (pam)",
  "agnes ward.D2",
  "agnes complete"
)

for (i in 1:length(unique(compare_algors$name))) {
  uniq_names <- unique(compare_algors$name)
  compare_algors$name[compare_algors$name == uniq_names[i]] <- algors[i]
}

svg("figs/Algorithms average silhouette.svg", height = 3.5, family = "serif")
ggplot(compare_algors, aes(x = as.factor(k), y = value, shape = name, group = name)) +
  geom_point(size = 1) +
  geom_line(size = 0.1) +
  guides(shape = guide_legend(
    override.aes = list(alpha = 1, size = 1),
    title = "Clustering Algorithm\n(Euclidean Distance)",
    labels = unique(compare_algors$name)
  )) +
  theme_bw(base_size = 11) +
  theme(
    plot.margin = margin(0, 0, 0, 0, "cm"),
    legend.position = c(0.87, 0.75),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 11),
    legend.background = element_rect(fill = "transparent")
  ) +
  ylab("Average silhouette width") +
  xlab("Number of clusters k")
dev.off()


# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#
#                                                                              #
#                           k-means clustering-----------------
#                                                                              #
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#


df[, var_index] <- scale_data(df[, var_index], method = "z-score")

best_kmeans_seeds <- c(
  k2 = 93038,
  k3 = 113663,
  k4 = 114088,
  k5 = 101955,
  k6 = 101766,
  k7 = 103747,
  k8 = 89527,
  k9 = 114078,
  k10 = 190583,
  k11 = 116827,
  k12 = 104093,
  k13 = 102839,
  k14 = 88195,
  k15 = 94217
)

centers_list <- list()

for (i in 1:length(best_kmeans_seeds)) {
  set.seed(best_kmeans_seeds[i])
  centers_list[[paste(best_kmeans_seeds[i])]] <-
    df[sample(which_not_NA, i + 1), var_index]
}

dist_mat <- parallelDist(as.matrix(df[which_not_NA, var_index]), "euclidean")

kmeans_res <- cluster.data(
  df = df[which_not_NA, var_index],
  k = 2:15,
  init_clus_centers = centers_list,
  cluster_method = "kmeans",
  hclust_method = NULL,
  dist_mat = dist_mat,
  hclust_tree = NULL,
  internal_metrics = TRUE,
  verbose = TRUE
)

kmeans_res$Avg.silhouette.width

# cluster sizes
sapply(1:7, function(k) length(which(kmeans_res$clustering$`k = 7` == k)))

clustering <- kmeans_res$clustering$`k = 7`

# names for k = 7
pattern_names <- c(
  "P+ L-",
  "R+ F-",
  "R- P- L+",
  "F+",
  "L+ G+ T+",
  "L+ G+",
  "Unnamed"
)

for (i in unique(clustering)) {
  clustering[clustering == i] <- pattern_names[i]
}

svg("figs/Clusters in parallel coord.svg", width = 10, height = 6, family = "serif")
plot.par_coord(df[which_not_NA, var_index],
               clustering,
               alpha_trasp = 0.05,
               scale_y_step = 1,
               base_p_size = 11,
               line_size = 0.08
) +
  theme(plot.title = element_blank()) +
  geom_hline(
    yintercept = 0,
    color = "grey30",
    linetype = "twodash"
  )
dev.off()

df_average_cl <- aggregate.data.frame(
  df[which_not_NA, var_index],
  list("cluster" = clustering),
  mean
)

svg("figs/Mean cluster values.svg", height = 3.5, family = "serif")
plot.par_coord(df_average_cl[, 2:7],
               df_average_cl$cluster,
               clust_plot_scheme = "shape",
               alpha_trasp = 1,
               scale_y_step = 1,
               base_p_size = 11,
               line_size = 0.08, legend_pos = c(0.3, 0.85)
) +
  scale_y_continuous(
    name = "z-score",
    limits = c(min(df_average_cl[2:7]), max(df_average_cl[2:7])),
    breaks = seq(round(min(df_average_cl[2:7])),
                 round(max(df_average_cl[2:7])),
                 by = 0.5
    ),
    minor_breaks = NULL,
    expand = c(0, 0.1)
  ) +
  theme(plot.title = element_blank())
dev.off()

silhouette_res <- as.data.frame(kmeans_res$Silhouette$`clustering k = 7`[, 1:3])
silhouette_res$cluster <- clustering

svg("figs/Silhoutte plot.svg", height = 4.5, family = "serif")
plot.silhouette(silhouette_res)
dev.off()

# return unnormalized data
df[which_not_NA, var_index] <- read.csv("data/corp data.csv")[which_not_NA, var_index]

df$Cluster <- 0
df$Silhouette <- NA

df$Cluster[which_not_NA] <- clustering
df$Silhouette[which_not_NA] <- silhouette_res$sil_width

df$Cluster[out1] <- "InfNaNNA"
df$Cluster[out2] <- "Boxplot outlier"

for (i in unique(silhouette_res$neighbor)) {
  silhouette_res$neighbor[silhouette_res$neighbor == i] <- pattern_names[i]
}

df$Silhouette_neighbour <- NA
df$Silhouette_neighbour[which_not_NA] <- silhouette_res$neighbor

write.xlsx(df, "output/Data clustering results.xlsx", keepNA = TRUE)

summary_clusters <- by(
  data = df[which_not_NA, var_index],
  INDICES = list(cluster = df$Cluster[which_not_NA]),
  FUN = descr,
  stats = c("min", "q1", "med", "mean", "sd", "q3", "max")
) %>%
  lapply(function(x) cbind.data.frame(x, rownames(x))) %>%
  do.call(rbind.data.frame, args = .) %>%
  cbind.data.frame(rep(levels(factor(df$Cluster[which_not_NA])), each = 7))

cluster_box_violins_plots <- by(
  data = df[which_not_NA, var_index],
  INDICES = list(cluster = df$Cluster[which_not_NA]),
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
            min(df[which_not_NA, i + 1]) - 0.05,
            max(df[which_not_NA, i + 1])
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
      labels = letters[1:ncol(x)],
      label.args = list(
        gp = grid::gpar(cex = 1)
      ),
      draw = FALSE
    )
  }
)

for (i in names(cluster_box_violins_plots)) {
  svg(paste0("figs/", i, " boxplot.svg"), family = "serif")
  print(cluster_box_violins_plots[[i]])
  dev.off()
}
