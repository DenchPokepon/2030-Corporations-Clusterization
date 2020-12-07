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
# This function takes the data and creates the violinplot with the boxplot on top of it
# If bnd is not NULL, then the data is filtered by the specified boundaries and only then is plotted
# If bnd is NULL returns a plot.
# Args:
#   df - data.frame or tibble;
#   bnd - boundaries of variables in df. Specified as follows:
#     lower boundary of variable 1, upper boundary of variable 1,
#     lower boundary of variable 2, upper boundary of variable 2,
#     ...,
#     lower boundary of variable n, upper boundary of variable n.
#   names_vars - vector of string names of variables in df;
#   base_p_size - base font size in ggplot theme;
#   outlier.alpha - alpha transparency for outliers points (0 = fully transparent, 1 = no transparency).


plot.box_violin <- function(df,
                            bnd = NULL,
                            arrange_p = TRUE,
                            names_vars = c("P", "R", "L", "G", "F", "T"),
                            base_p_size = 11,
                            outlier.alpha = 0.05) {
  library(ggplot2)
  library(egg)
  library(tibble)
  library(grid)
  library(scales)

  if (!(is.vector(bnd) & length(bnd) == length(names_vars) * 2 | is.null(bnd))) {
    stop(paste("bnd should be a vector of length", length(names_vars) * 2, "or NULL"))
  }

  plist <- list()
  if (!is.null(bnd)) {
    condition <- TRUE
    for (i in 1:length(names_vars)) {
      condition <- df[names_vars[i]] >= bnd[i + i - 1] &
        df[names_vars[i]] <= bnd[i + i] & condition
    }
    out <- df[which(!condition), ]
    df <- df[which(condition), ]
  } else {
    out <- NULL
  }

  for (i in 1:length(names_vars)) {
    plist[[i]] <- ggplot(df, aes_string(y = names_vars[i], x = rep(1, nrow(df)))) +
      geom_violin(trim = FALSE, fill = "gray", colour = "grey") +
      geom_boxplot(width = 0.1, outlier.size = 0.5, outlier.alpha = outlier.alpha) +
      theme_bw(base_size = base_p_size) +
      xlab(names_vars[i]) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(1, 5, 5, 5, "pt")
      ) +
      scale_y_continuous(breaks = pretty_breaks()) +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_blank()
      )
  }
  p <- egg::ggarrange(
    plots = plist,
    ncol = 3,
    nrow = 2,
    # labels = letters[1:length(names_vars)],
    # label.args = list(
    #   gp = gpar(cex = 1)
    # ),
    draw = FALSE
  )

  if (arrange_p) {
    p <- p
  } else if (!arrange_p) {
    p <- plist
  }

  return(p)
}


# This function is a collection of data normalization / standardization functions
# Input - data frame or tibble with variables to normalize
# Returns an exact copy of the table with normalized / standardized variables


scale_data <- function(df,
                       method = c("z-score", "min-max", "min-optimum-max"),
                       optimum_values = NULL) {
  method <- match.arg(method, method)

  if (method == "z-score") {
    df[, 1:ncol(df)] <- apply(df, 2, scale)
  }
  else if (method == "min-max") {
    df[, 1:ncol(df)] <- apply(
      df,
      2,
      function(x) (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))
    )
  }
  else if (method == "min-optimum-max") {
    if (inherits(optimum_values, "numeric")) opt_x <- t(optimum_values)
    df[, 1:ncol(df)] <- mapply(
      df,
      opt_x,
      FUN = function(x, x_opt) {
        1 - ((x - x_opt) /
          max(
            x_opt - min(x, na.rm = T),
            max(x, na.rm = T) - x_opt
          )
        )
      }
    )
  }
  return(df)
}


# Creates a ggplot parallel coordinates plot
# Args:
#   df - data.frame or tibble;
#   clustering - clustering vector (can be numeric or character/factor);
#   clust_plot_scheme - how to plot different clusters
#   names_vars - vector of string names of variables in df;
#   alpha_transp - alpha transparency for lines (0 = fully transparent, 1 = no transparency);
#   line_size - width of hyperlines;
#   base_p_size - base font size in ggplot theme;
#   y_units - units of measurement of y axis (used for y axis name);
#   scale_y_step - measurement scale step;
#   colours - color palette for cluster highlighting;
#   legend_pos - legend location specification for ggplot theme legend.position


plot.par_coord <- function(df,
                           clustering = NULL,
                           clust_plot_scheme = c(
                             "facet_wrap",
                             "colour",
                             "shape",
                             "highlight"
                           ),
                           names_vars = c("P", "R", "L", "G", "F", "T"),
                           alpha_trasp = 1,
                           line_size = 1,
                           base_p_size = 11,
                           y_units = "Z-score",
                           scale_y_step = 0.5,
                           colours = c("cyan", RColorBrewer::brewer.pal(8, name = "Set1")),
                           legend_pos = c(0.87, 0.2)) {
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(RColorBrewer)

  clust_plot_scheme <- match.arg(clust_plot_scheme, clust_plot_scheme)

  df$dummy <- NA
  if (is.null(clustering)) {
    vars_cols <- which(names(df) %in% c(names_vars, "dummy"))
  } else {
    df$cluster <- clustering
    vars_cols <- which(names(df) %in% c(names_vars, "dummy", "cluster"))
  }
  df <- df[vars_cols]
  dummy_vector <- rep(c(1:length(names_vars), NA), nrow(df))
  pivot_data <- pivot_longer(data = df, cols = which(names(df) %in% c(
    names_vars,
    "dummy"
  ))) %>%
    cbind(dummy_vector)
  v_len <- length(names_vars)

  p <- ggplot(pivot_data, aes(x = dummy_vector, y = value)) +
    scale_x_continuous(
      name = "Variable",
      limits = c(1, v_len),
      breaks = 1:v_len,
      labels = names_vars,
      minor_breaks = NULL,
      expand = c(0.01, 0.01)
    ) +
    scale_y_continuous(
      name = y_units,
      limits = c(min(df[1:v_len]), max(df[1:v_len])),
      breaks = seq(
        round(min(df[1:v_len])),
        round(max(df[1:v_len])),
        by = scale_y_step
      ),
      minor_breaks = NULL,
      expand = c(0.01, 0.01)
    ) +
    theme_bw(base_size = base_p_size) +
    theme(plot.margin = margin(0, 0, 0, 0, "cm"))

  if (!is.null(clustering)) {
    if (clust_plot_scheme == "facet_wrap") {
      p <- p +
        geom_path(
          aes(y = value),
          alpha = alpha_trasp,
          na.rm = TRUE,
          size = line_size
        ) +
        facet_wrap(~ cluster) +
        theme(
          panel.spacing = unit(0.2, "lines"),
          strip.text = element_text(
            size = base_p_size,
            margin = margin(base_p_size / 420, 0, base_p_size / 420, 0, "cm")
          )
        )
    } else if (clust_plot_scheme == "colour") {
      p <- p +
        geom_path(aes(
          y = value,
          colour = cluster
        ),
        alpha = alpha_trasp,
        na.rm = TRUE,
        size = line_size
        ) +
        scale_color_manual(values = colours, labels = cluster, name = "") +
        guides(
          colour = guide_legend(
            override.aes = list(alpha = 1, size = 1),
            title = "",
            labels = cluster
          )
        ) +
        theme(
          legend.key.size = unit(0.1, "cm"),
          plot.margin = margin(0, 0, 0, 0, "cm"),
          legend.position = legend_pos,
          legend.text = element_text(size = base_p_size),
          legend.background = element_rect(fill = "transparent")
        )
    } else if (clust_plot_scheme == "shape") {
      p <- p +
        geom_point(aes(
          y = value,
          shape = cluster
        ),
        size = 1.5,
        alpha = alpha_trasp,
        na.rm = TRUE
        ) +
        geom_path(aes(
          y = value
        ),
        alpha = alpha_trasp,
        na.rm = TRUE,
        size = line_size
        ) +
        scale_shape_manual(values = c(16, 17, 15, 3, 7, 8, 10, 14, 11)) +
        guides(
          shape = guide_legend(
            override.aes = list(alpha = 1, size = 1.5),
            title = "",
            labels = cluster
          )
        ) +
        theme(
          legend.key.size = unit(0.1, "cm"),
          plot.margin = margin(0, 0, 0, 0, "cm"),
          legend.position = legend_pos,
          legend.text = element_text(size = base_p_size),
          legend.background = element_rect(fill = "transparent")
        )
    } else if (clust_plot_scheme == "highlight") {
      highlight <- pivot_data[grep("Cluster", pivot_data$cluster), ]

      p <- p +
        geom_path(aes(
          y = value,
          color = "gray",
          linetype = "solid"
        ),
        alpha = alpha_trasp,
        na.rm = TRUE,
        size = line_size
        ) +
        geom_path(aes(
          y = value,
          color = "black",
          linetype = "longdash",
        ),
        data = highlight,
        alpha = 0.2,
        na.rm = TRUE,
        size = line_size
        ) +
        scale_color_identity(
          name = "Highlighted cluster",
          breaks = c("gray", "black"),
          labels = c("Main", highlight$cluster[1]),
          guide = "legend"
        ) +
        scale_linetype_identity(
          name = "Highlighted cluster",
          breaks = c("solid", "longdash"),
          labels = c("Main", highlight$cluster[1]),
          guide = "legend"
        ) +
        theme(
          legend.key.size = unit(0.8, "cm"),
          plot.margin = margin(0, 0, 0, 0, "cm"),
          legend.position = legend_pos,
          legend.text = element_text(size = base_p_size),
          legend.background = element_rect(fill = "transparent")
        )
    }
  } else {
    p <- p +
      geom_path(aes(y = value),
        alpha = alpha_trasp,
        na.rm = TRUE,
        size = line_size
      )
  }

  return(p)
}


# performs data partitioning into a certain number of clusters and returns basic validation metrics
# Args:
#   df - data.frame or tibble;
#   k - number of clusters, specify vector of k values to split data into different number of clusters in one function call;
#   init_clus_centers - init cluster centers for k-means or medoids for pam;
#   cluster_method - clustering algorithm;
#   hclust_method - hierarchical clustering method of linkage;
#   dist_mat - precomputed distance matrix;
#   hclust_tree - precomputed hierarchical clustering tree;
#   max_iter - maximum number of iterations for k-means to converge;
#   internal_metrics - logical, calculate silhouette and WSS or not;
#   verbose - logical, print clustering progress or not.


cluster_data <- function(df,
                         k,
                         init_clus_centers = NULL,
                         cluster_method = c(
                           "kmeans",
                           "hclust",
                           "pam"
                         ),
                         hclust_method = c(
                           "ward.D2",
                           "ward.D",
                           "single",
                           "complete",
                           "average",
                           "mcquitty",
                           "centroid",
                           "median"
                         ),
                         dist_mat = NULL,
                         hclust_tree = NULL,
                         max_iter = 20,
                         internal_metrics = FALSE,
                         verbose = TRUE) {
  library(factoextra)
  library(cluster)

  # Input check-----------

  if (!all(apply(df, 2, is.numeric))) {
    stop("All columns in df should be numeric type")
  }

  cluster_method <- match.arg(cluster_method, cluster_method)
  hclust_method <- match.arg(hclust_method, hclust_method)

  if (is.null(dist_mat)) {
    dist_mat <- parallelDist::parallelDist(df, method = "euclidean")
  }
  if (is.null(hclust_tree) & cluster_method == "hclust") {
    hclust_tree <- fastcluster::hclust(dist_mat, method = hclust_method)
  }

  cluster_res <- list(
    "cluster_method" = cluster_method,
    "hclust_method" = hclust_method,
    "main_metric" = "euclidean",
    "clustering" = list()
  )

  wss_vec <- sil_vec <- numeric(k[length(k)])
  sil_list <- list()
  names(wss_vec) <- names(sil_vec) <- paste("k =", 1:k[length(k)])
  int_metr <- internal_metrics

  for (kk in k) {
    if (verbose) print(paste("clustering", kk, "clusters"))

    # Classic kmeans, Hartigan-Wong algorithm if init centers not supplied------

    if (cluster_method == "kmeans") {
      if (is.null(init_clus_centers)) {
        clustering <- kmeans(df,
          centers = kk,
          iter.max = max_iter
        )
      } else if (length(init_clus_centers) == length(k)) {
        clustering <- kmeans(df,
          centers = init_clus_centers[[kk - 1]],
          iter.max = max_iter
        )
      } else {
        stop("kmeans_centers length", length(init_clus_centers), "!= k length", length(k))
      }

      if (int_metr) {
        cluster_wss <- clustering$tot.withinss
        cluster_sil <- silhouette(clustering$cluster, dist_mat)
      } else {
        cluster_wss <- NA
        cluster_sil <- NA
      }

      clustering <- clustering$cluster
    }

    # Hierarchical clustering-------------

    else if (cluster_method == "hclust") {
      clustering <- cutree(tree = hclust_tree, kk)

      euclidean_wss <- function(x, clustering) {
        spl <- split(x, clustering)

        wss <- sum(sapply(spl, function(d) sum(scale(d, scale = FALSE)^2)))
        return(wss)
      }
      if (int_metr) {
        cluster_wss <- euclidean_wss(df, clustering)
        cluster_sil <- silhouette(clustering, dist_mat)
      } else {
        cluster_wss <- NA
        cluster_sil <- NA
      }
    }

    # K-medoids with optimization from Schubert and Rousseeuw (2019)---------

    else if (cluster_method == "pam") {
      if (is.null(init_clus_centers)) {
        clustering <- clustering <- pam(
          x = dist_mat,
          k = kk,
          diss = TRUE,
          pamonce = 3
        )
      } else if (length(init_clus_centers) == length(k)) {
        cluster_centers_index <- which(sapply(init_clus_centers, nrow) == kk)
        clustering <- clustering <- pam(
          x = dist_mat,
          k = kk,
          diss = TRUE,
          pamonce = 3,
          medoids = which(row.names(df) %in%
            row.names(init_clus_centers[[cluster_centers_index]]))
        )
      } else {
        stop("pam_medoids length", length(init_clus_centers), "!= k length", length(k))
      }

      euclidean_wss <- function(x, clustering) {
        spl <- split(x, clustering)

        wss <- sum(sapply(spl, function(d) sum(scale(d, scale = FALSE)^2)))
        return(wss)
      }
      if (int_metr) {
        cluster_wss <- euclidean_wss(df, clustering$cluster)
        cluster_sil <- silhouette(clustering$cluster, dist_mat)
      } else {
        cluster_wss <- NA
        cluster_sil <- NA
      }

      clustering <- clustering$cluster
    }

    cluster_res$clustering[[paste("k =", kk)]] <- clustering
    wss_vec[kk] <- cluster_wss
    sil_vec[kk] <- ifelse(int_metr, mean(cluster_sil[, 3]), NA)
    sil_list[[paste("clustering k =", kk)]] <- cluster_sil
  }

  cluster_res$WSS <- wss_vec
  cluster_res$Avg.silhouette.width <- sil_vec
  cluster_res$Silhouette <- sil_list

  return(cluster_res)
}


# Plots ggplot visualization of silhouette information from cluster::silhouette() data.frame
# Args:
#   df - as.data.frame(silhouette()[, 1:3]) object;
#   names_just - x axis justification for cluster names


plot.silhouette <- function(df, names_just = 0.5, name_angle = 45, base_p_size = 11) {
  library(ggplot2)
  library(dplyr)
  library(scales)

  # order by cluster and by sil_width
  df <- df[order(df$cluster, -df$sil_width), ]
  if (!is.null(rownames(df))) {
    df$name <- factor(rownames(df), levels = rownames(df))
  } else {
    df$name <- as.factor(1:nrow(df))
  }
  df$name <- df$name
  df$cluster <- as.factor(df$cluster)
  df_cl_mean <- aggregate(df$sil_width, by = list("cluster" = df$cluster), mean) %>%
    rename(ystart = x) %>%
    mutate(yend = ystart)

  df_clusters <- split(df$name, df$cluster)

  df_cl_mean <- df_clusters %>%
    lapply(function(x) x[c(1, length(x))]) %>%
    do.call(rbind.data.frame, args = .) %>%
    `names<-`(c("xstart", "xend")) %>%
    cbind.data.frame(df_cl_mean)

  middle_x_values <- df_clusters %>%
    lapply(function(x) x[length(x) %/% 2]) %>%
    do.call(rbind.data.frame, args = .) %>%
    `names<-`("x") %>%
    mutate("cluster" = levels(df$cluster), "y" = rep(names_just, nrow(.)))

  p <- ggplot(df, aes(x = name, y = sil_width)) +
    geom_bar(
      stat = "identity",
      fill = "gray",
      colour = "gray"
    ) +
    labs(
      y = "Silhouette Width", x = "Company in a Particular Year"
    ) +
    geom_segment(aes(x = xstart, xend = xend, y = ystart, yend = yend),
      df_cl_mean,
      linetype = "dashed", color = "black"
    ) +
    geom_text(aes(x = x, 
                  y = y, 
                  label = cluster),
      middle_x_values,
      angle = name_angle,
      size = base_p_size / ggplot2:::.pt
    ) +
    theme_classic(base_size = base_p_size) +
    scale_y_continuous(
      limits = c(
        min(df$sil_width),
        max(df$sil_width)
      ),
      breaks = pretty_breaks(10)
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_blank()
    )

  return(p)
}
