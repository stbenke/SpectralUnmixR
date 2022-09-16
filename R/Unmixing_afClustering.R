#' Select Autofluorescence Clusters
#'
#' For a clustering of autofluorescence populations, find similar clusters based on the cosine similarity.
#'
#' @param dat data read in using GetData
#' @param clusters vector with cluster assignment for each entry in dat
#' @param channels character vector, channel selection to include in the plot (default: NULL, i.e. all)
#' @param count_min integer, minimum number of events needed for a cluster to be considered in the selection (default: NULL, i.e. all)
#' @param csim_max numeric, maximum cosine similarity value. Clusters with a higher cosine similarity value will be combined. (default: NULL, i.e. no combining of any clusters)
#' @param selection either "max" or "median". If set to "max", the cluster with the highest intensity is chosen from a group of similar clusters as spectrum. If set to "median", the median of all events within a group of similar clusters is taken as spectrum.
#'
#' @return List of outputs. The spectra element contains the dataframe with the selected spectra.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
ClusterSelection <- function(dat, clusters,
                             channels = NULL,
                             count_min = NULL, csim_max = NULL,
                             selection = "max") {

  channels_used <- ChannelCheck(dat, channels)


  cluster_counts <- as.numeric(table(clusters))

  # add clustering result, calculate median signals
  dat_clust <- dat %>%
    dplyr::select(dplyr::all_of(channels_used)) %>%
    dplyr::mutate(cluster = clusters) %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::summarise_at(dplyr::all_of(channels_used), function(x) stats::median(x)) %>%
    dplyr::ungroup() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(max = max(dplyr::c_across(dplyr::all_of(channels_used)))) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data$cluster) %>%
    dplyr::mutate(cluster_count = cluster_counts)

  # remove small clusters
  if (!is.null(count_min)) {
    dat_clust <- dplyr::filter(dat_clust,
                               .data$cluster_count > count_min)

  }

  # find similar clusters
  if (!is.null(csim_max)) {

    CosSim <- function(m) {
      sim <- m / sqrt(rowSums(m * m))
      sim %*% t(sim)
    }

    sim <- CosSim(dat_clust %>%
                    dplyr::select(dplyr::all_of(channels_used)) %>%
                    as.matrix)
    diag(sim) <- 0

    inds <- which(sim > csim_max, arr.ind = T)

    groups <- purrr::map(1:nrow(inds), function(i) inds[i,])

    j <- 1
    while (j <= length(groups)) {
      i <- j+1
      while (i <= length(groups)) {
        if (any(groups[[j]] %in% groups[[i]])) {
          groups[[j]] <- unique(c(groups[[j]], groups[[i]]))
          groups <- groups[-i]
          break
        } else {
          i <- i+1
        }
      }
      if (i > length(groups)) {
        j <- j+1
      }
    }

    # convert to original cluster names
    groups <- purrr::map(groups, function(x) {
      dat_clust$cluster[x]
    })

    for (clust in dat_clust$cluster) {
      if (!clust %in% unlist(groups)) {
        groups <- c(groups, list(clust))
      }
    }

    if (selection == "max") {
      clust_max <- dat_clust$max
      names(clust_max) <- dat_clust$cluster

      cluster_sel <- purrr::map_dbl(groups, function(x) {
        x[which.max(clust_max[x])]
      })
    }

    cluster_groups <- purrr::map_dfr(seq_along(groups), function(i) {
      tibble::tibble(group = i, cluster = groups[[i]])
    })
  }

  normFunc <- function(x, max) x/max

  list(
    spectra = if (!is.null(csim_max)) {
      if (selection == "max") {
        dat_clust %>%
          dplyr::select(dplyr::all_of(channels_used)) %>%
          as.matrix %>%
          apply(1, function(x) x/max(x)) %>%
          t() %>%
          tibble::as_tibble() %>%
          dplyr::mutate(cluster = dat_clust$cluster) %>%
          dplyr::filter(.data$cluster %in% cluster_sel) %>%
          dplyr::rename(file = .data$cluster) %>%
          dplyr::select(dplyr::all_of(c(channels_used, "file")))
      } else if (selection == "median") {
        dat %>%
          dplyr::select(dplyr::all_of(channels_used)) %>%
          dplyr::mutate(cluster = clusters) %>%
          dplyr::filter(if (!is.null(count_min)) {
            .data$cluster %in% unique(dat_clust$cluster)
          } else {
            TRUE
          }) %>%
          dplyr::left_join(cluster_groups, by = "cluster") %>%
          dplyr::group_by(.data$group) %>%
          dplyr::summarise_at(dplyr::all_of(channels_used),
                              function(x) stats::median(x)) %>%
          dplyr::ungroup() %>%
          dplyr::arrange(.data$group) %>%
          dplyr::select(dplyr::all_of(channels_used)) %>%
          as.matrix %>%
          apply(1, function(x) x/max(x)) %>%
          t() %>%
          tibble::as_tibble() %>%
          dplyr::mutate(file = sort(unique(cluster_groups$group)))
      }
    } else {
      dat_clust %>%
        dplyr::select(dplyr::all_of(channels_used)) %>%
        as.matrix %>%
        apply(1, function(x) x/max(x)) %>%
        t() %>%
        tibble::as_tibble() %>%
        dplyr::mutate(cluster = dat_clust$cluster) %>%
        dplyr::rename(file = .data$cluster) %>%
        dplyr::select(dplyr::all_of(c(channels_used, "file")))
    },
    cluster_selection = if (!is.null(csim_max) && selection == "max") cluster_sel else NULL,
    cluster_groups = if (!is.null(csim_max)) cluster_groups else NULL,
    count_min = count_min,
    csim_max = csim_max,
    selection = selection,
    data_cluster = dat_clust,
    data = dat,
    channels_used = channels_used,
    clusters = clusters)
}


#' Plot clustered AF spectra
#'
#' @param input result of ClusterSelection
#' @param transformation numeric, cofactor used for the asinh transformation (default: 0.1) Note: The spectra are max normalized to 1 so your usual cofactors will not work here.
#' @param bins numeric, number of bins to use in the channel heatmaps (default: 50)
#' @param ncol integer, number of columns of the facet plot (default: NULL, i.e. automatic)
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return ggplot object
#' @export
#'
PlotSpectraClusters <- function(input, transformation = 0.1,
                                bins = 50, ncol = 3) {

  if (!is.null(input$csim_max)) {
    dat_cluster <- input$data_cluster %>%
      dplyr::left_join(input$cluster_groups, by = "cluster") %>%
      dplyr::rename(file = .data$group)
  } else {
    dat_cluster <- input$data_cluster %>%
      dplyr::mutate(file = .data$cluster)
  }

  spectra_prep_all <- dat_cluster %>%
    tidyr::pivot_longer(cols = dplyr::all_of(input$channels_used),
                        values_to = "mean",
                        names_to = "channel") %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::mutate(mean = .data$mean / max(.data$mean)) %>%
    dplyr::mutate(channel = factor(.data$channel, input$channels_used),
                  mean = asinh(mean/transformation))

  if (!is.null(input$csim_max)) {
    if (input$selection == "max") {
      spectra_prep_sel <- input$spectra %>%
        dplyr::left_join(input$cluster_groups, by = c("file" = "cluster")) %>%
        dplyr::mutate(file = .data$group)
    } else if (input$selection == "median") {
      spectra_prep_sel <- input$spectra
    }
    spectra_prep_sel <- spectra_prep_sel %>%
      dplyr::mutate(cluster = (max(input$clusters)+1):(max(input$clusters)+nrow(input$spectra))) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(input$channels_used),
                          values_to = "mean",
                          names_to = "channel") %>%
      dplyr::mutate(channel = factor(.data$channel, input$channels_used),
                    mean = asinh(mean/transformation))
  } else {
    spectra_prep_sel <- input$spectra %>%
      dplyr::mutate(cluster = .data$file) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(input$channels_used),
                          values_to = "mean",
                          names_to = "channel") %>%
      dplyr::mutate(channel = factor(.data$channel, input$channels_used),
                    mean = asinh(mean/transformation))
  }

  dat_prep <- input$data %>%
    dplyr::mutate(cluster = input$clusters)

  if (!is.null(input$csim_max)) {
    dat_prep <- dat_prep %>%
      dplyr::left_join(input$cluster_groups, by = "cluster") %>%
      dplyr::select(-file) %>%
      dplyr::rename(file = .data$group)
  }

  if (!is.null(input$count_min)) {
    cluster_count <- table(input$clusters)
    sel <- names(cluster_count[cluster_count > input$count_min])
    dat_prep <- dplyr::filter(dat_prep,
                              .data$cluster %in% sel)
  }

  dat_prep <- purrr::map_dfr(unique(dat_prep$file), function(x) {
    dat_prep %>%
      dplyr::filter(file == x) %>%
      SampleCorr(channels = input$channels_used)
  })

  PlotSpectra(dat_prep,
              channels = input$channels_used,
              bins = bins,
              transformation = transformation,
              ncol = ncol) +
    ggplot2::geom_line(data = dplyr::filter(spectra_prep_all),
                       ggplot2::aes(group = .data$cluster),
                       color = "black", size = 0.5) +
    ggplot2::geom_line(data = dplyr::filter(spectra_prep_sel),
                       ggplot2::aes(group = .data$cluster),
                       color = "red", size = 0.75)
}


#' Plot scatter signals of AF clusters
#'
#' @param input result of ClusterSelection
#' @param scatter_channels character vector, names of the two scatter channels to plot (default: c("FSC-A", "SSC-A"))
#' @param transformation numeric, cofactor used for the asinh transformation (default: NULL, i.e., none)
#' @param bins numeric, number of bins to use in the 2D histograms (default: 100)
#' @param alpha numeric, transparency parameter of the 2D histograms (default: 0.5)
#' @param ncol integer, number of columns of the facet plot (default: NULL, i.e. automatic)
#'
#' @importFrom rlang .data
#'
#' @return ggplot object, showing the scatter signals for each group of clusters (colored by cluster) found by ClusterSelection. If a group labelled NA is shown in the plot, clusters were excluded from the analysis in ClusterSelection for having too few events.
#' @export
#'
PlotScatterClusters <- function(input,
                                scatter_channels = c("FSC-A", "SSC-A"),
                                transformation = NULL,
                                bins = 100, alpha = 0.5, ncol = 3) {
  if (any(!c("data", "clusters", "cluster_groups") %in% names(input))) {
    stop("input must be result of ClusterSelection")
  }

  dat <- input$data

  dat$cluster <- input$clusters

  if (!is.null(input$cluster_groups)) {
    dat <- dplyr::left_join(dat, input$cluster_groups, by = "cluster")
  } else {
    dat <- dplyr::mutate(dat, group = 1)
  }

  if (!is.null(transformation)) {
    dat <- dplyr::mutate_at(dplyr::all_of(scatter_channels, function(x) {
      asinh(x/transformation)
    }))
  }

  ggplot2::ggplot(dat,
                  ggplot2::aes(.data[[scatter_channels[1]]],
                               .data[[scatter_channels[2]]],
                               fill = as.factor(.data$cluster))) +
    ggplot2::geom_hex(bins = bins, alpha = alpha) +
    ggplot2::facet_wrap(ggplot2::vars(.data$group), ncol = ncol) +
    ggplot2::theme_bw() +
    ggplot2::labs(fill = "Cluster")
}
