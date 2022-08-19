
# File handling -----------------------------------------------------------


#' Read data from FCS file into tibble or convert flowFrame into tibble
#'
#' @param f Path to FCS file(s), flowFrame or flowSet
#' @param pattern Regular expression to extract pattern from file name (default: NULL)
#' @param file_name File name(s) to be stored in tibble (default: NULL)
#' @param ... Arguments passed to read.FCS
#'
#' @return Tibble with the exprs matrix of the FCS file and a column containing the file name or extracted pattern. If both pattern and file_name are NA, the literal filename will be used.
#' @export
#' @importFrom magrittr %>%
#'
GetData <- function(f, pattern = NULL, file_name = NULL, ...) {

  if (methods::is(f, "flowFrame")) {
    dat <- flowCore::exprs(f)
    file <- flowCore::keyword(f, "$FIL")[[1]]

    out <- dat %>%
      tibble::as_tibble() %>%
      dplyr::mutate(file = if (!is.null(pattern)) {
        stringr::str_extract(file, pattern)
      } else if (!is.null(file_name)) {
        file_name
      } else {
        file
      }
      )
  } else if (methods::is(f, "flowSet")) {
    purrr::map_dfr(seq_along(f), function(i) {
      dat <- flowCore::exprs(f[[i]])
      file <- flowCore::keyword(f[[i]], "$FIL")[[1]]

      out <- dat %>%
        tibble::as_tibble() %>%
        dplyr::mutate(file = if (!is.null(pattern)) {
          stringr::str_extract(file, pattern)
        } else if (!is.null(file_name)) {
          if (length(file_name) == length(f)) {
            file_name[i]
          } else {
            stop("Provided file names must be as many as files in flowSet.")
          }
        } else {
          file
        }
        )
    })
  }

  if (all(is.character(f))) {
    if (length(f) == 1) {
      dat <- flowCore::exprs(flowCore::read.FCS(f, ...))
      file <- f

      out <- dat %>%
        tibble::as_tibble() %>%
        dplyr::mutate(file = if (!is.null(pattern)) {
          stringr::str_extract(f, pattern)
        } else if (!is.null(file_name)) {
          file_name
        } else {
          f
        }
        )
    } else {
      out <- purrr::map_dfr(seq_along(f), function(i) {
        flowCore::exprs(flowCore::read.FCS(f[i], ...)) %>%
          tibble::as_tibble() %>%
          dplyr::mutate(file = if (!is.null(pattern)) {
            stringr::str_extract(f[i], pattern)
          } else if (!is.null(file_name)) {
            if (length(file_name) == length(f)) {
              file_name[i]
            } else {
              stop("Provided file names must be as many as files.")
            }
          } else {
            f[i]
          }
          )
      })
    }
  } else {
    stop("f must be either a flowFrame, flowSet or path to one or multiple FCS file(s).")
  }

  if (any(is.na(out$file))) {
    warning("NAs produced for filenames. Most likely, there is a problem with the pattern argument.")
  }

  out
}

# read in FCS file and get median values for each channel.
# only apply to properly gated files!
#' Get channel medians for FCS file or flowFrame
#'
#' @param f Path to FCS file(s), flowFrame or flowSet, or output of GetData
#' @param pattern Regular expression to extract pattern from file name (default: NA)
#' @param file_name File name to be stored in tibble
#' @param ... Arguments passed to read.FCS
#'
#' @return Tibble with the channel medians in the columns with additional file column.
#' @export
#' @importFrom magrittr %>%
#'
GetMedianData <- function(f, pattern = NULL, file_name = NULL, ...) {

  if (methods::is(f, "flowFrame")) {
    dat <- flowCore::exprs(f)
    file <- flowCore::keyword(f, "$FIL")[[1]]

    out <- dat %>%
      tibble::as_tibble() %>%
      dplyr::summarise_all(function(x) stats::median(x)) %>%
      dplyr::mutate(file = if (!is.null(pattern)) {
        stringr::str_extract(file, pattern)
      } else if (!is.null(file_name)) {
        file_name
      } else {
        file
      }
      )
  } else if (methods::is(f, "flowSet")) {
    purrr::map_dfr(seq_along(f), function(i) {
      dat <- flowCore::exprs(f[[i]])
      file <- flowCore::keyword(f[[i]], "$FIL")[[1]]

      out <- dat %>%
        tibble::as_tibble() %>%
        dplyr::summarise_all(function(x) stats::median(x)) %>%
        dplyr::mutate(file = if (!is.null(pattern)) {
          stringr::str_extract(file, pattern)
        } else if (!is.null(file_name)) {
          if (length(file_name) == length(f)) {
            file_name[i]
          } else {
            stop("Provided file names must be as many as files in flowSet.")
          }
        } else {
          file
        }
        )
    })
  } else if (all(is.character(f))) {
    if (length(f) == 1) {
      dat <- flowCore::exprs(flowCore::read.FCS(f, ...))
      file <- f

      out <- dat %>%
        tibble::as_tibble() %>%
        dplyr::summarise_all(function(x) stats::median(x)) %>%
        dplyr::mutate(file = if (!is.null(pattern)) {
          stringr::str_extract(f, pattern)
        } else if (!is.null(file_name)) {
          file_name
        } else {
          f
        }
        )
    } else {
      out <- purrr::map_dfr(seq_along(f), function(i) {
        flowCore::exprs(flowCore::read.FCS(f[i], ...)) %>%
          tibble::as_tibble() %>%
          dplyr::summarise_all(function(x) stats::median(x)) %>%
          dplyr::mutate(file = if (!is.null(pattern)) {
            stringr::str_extract(f[i], pattern)
          } else if (!is.null(file_name)) {
            if (length(file_name) == length(f)) {
              file_name[i]
            } else {
              stop("Provided file names must be as many as files.")
            }
          } else {
            f[i]
          }
          )
      })
    }
  } else if (is.data.frame(f)) {
    if ("file" %in% colnames(f)) {
      out <- f %>%
        dplyr::group_by(file) %>%
        dplyr::summarise_all(function(x) stats::median(x))
    } else {
      stop("file column not in dataframe. Dataframe must be output of GetData.")
    }
  } else {
    stop("f must be either a flowFrame, flowSet or path to one or multiple FCS file(s).")
  }

  if (any(is.na(out$file))) {
    warning("NAs produced for filenames. Most likely, there is a problem with the pattern argument.")
  }

  out
}



#  Spectra calculation ----------------------------------------------------


#' Check if channel list is given and present in data
#'
#' Internal function
#'
#' @param dat tibble, output of GetData or GetMedianData
#' @param channels character vector, channels to be used.
#'
#' @return list of channels to be used
ChannelCheck <- function(dat, channels) {
  if (is.null(channels)) {
    channels_used <- colnames(dat)
    channels_used[channels_used != "file"]
  } else {
    if (all(channels %in% colnames(dat))) {
      channels
    } else {
      stop("The following channels were not found in the data:\n",
           stringr::str_c(channels[!channels %in% colnames(dat)],
                          collapse = "\n"))
    }
  }
}


#' Calculate reference spectra by median subtraction
#'
#' Classical way of calculating reference spectra by subtracting the medians of a negative population from the medians of a positive population and then normalizing to the highest channel reading.
#' Multiple dyes can be processed in one go, as long as they share the same negative population.
#'
#' @param pos Positive population(s) medians supplied as dataframe containing the channel medians in the columns. Can contain multiple rows.
#' @param neg Negative population medians supplied as dataframe containing the channel medians in the columns. Can contain only one row. If no negative population is supplied, only channel normalization is done (e.g., for autofluorescence populations).
#' @param channels Optional, the channels (column names) to include in the calculation.
#' @param force_non_negative TRUE (default) / FALSE. If TRUE, sets all negative values to 0. Careful: negative values (if not very small) usually indicate a mismatch between the positive and negative population autofluorescence!
#'
#' @return A tibble containing the reference spectra (channel values in columns) together with file column taken from pos.
#' @export
#' @importFrom magrittr %>%
#'
MedianSpectra <- function(pos, neg = NULL, channels = NULL,
                          force_non_negative = T) {

  # check pos and neg having same channels
  if (!is.null(neg)) {
    if (!all(colnames(pos) %in% colnames(neg))) {
      stop("Channel names in pos and neg do not match.")
    }
  }

  # check channels
  channels_used <- ChannelCheck(pos, channels)

  # take required channels and convert to matrix
  tmp_pos <- pos %>%
    dplyr::select(tidyselect::all_of(channels_used)) %>%
    as.matrix

  if (!is.null(neg)) {
    if (nrow(neg) > 1) {
      stop("Only one negative population allowed per call of MedianSpectra.")
    }
    tmp_neg <- neg %>%
      dplyr::select(tidyselect::all_of(channels_used)) %>%
      unlist
  }

  Correct <- function(row) {
    # only correct if negative is given, otherwise only normalize
    if (is.null(neg)) {
      pos_corr <- row
    } else {
      pos_corr <- row - tmp_neg
    }

    pos_corr / max(pos_corr)
  }

  # apply correction and normalization
  res <- apply(tmp_pos, 1, Correct) %>%
    t %>%
    magrittr::set_colnames(channels_used)

  # check for negative values
  for (i in 1:nrow(res)) {
    if (any(res[i,] < 0)) {
      warning(stringr::str_c("Warning: ", sum(res[i,] < 0),
                             " spectra values below 0 for ", pos$file[i]))

      if (force_non_negative) {
        # if negative values result, force to 0
        res[i, res[i,] < 0] <- 0
      }
    }
  }

  # return tibble
  res %>%
    tibble::as_tibble() %>%
    dplyr::mutate(file = pos$file)
}


#' Calculate reference spectra by regression
#'
#' Calculates reference spectra based on regression between channels (like AutoSpill but without the additional optimization). Does not require pre-gating on positive and negative populations, but positive and negative events should ideally be present to similar amounts.
#'
#' @param pos Positive population(s) medians supplied as dataframe containing the channel medians in the columns. Can contain multiple rows.
#' @param channels Optional, the channels (column names) to include in the calculation.
#' @param neg Optional, negative population medians supplied as dataframe containing the channel medians in the columns. Can contain only one row. Needed if pos contains to few negative events.
#' @param force_non_negative TRUE (default) / FALSE. If TRUE, sets all negative values to 0. Careful: negative values (if not very small) usually indicate a mismatch between the positive and negative population autofluorescence!
#' @param robust TRUE (default) / FALSE. If TRUE, rebust regression is done via MASS::rlm, otherwise standard regression is done using lm.
#'
#' @return A tibble containing the reference spectra (channel values in columns) together with file column taken from pos.
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
RegSpectra <- function(pos, channels = NULL, neg = NULL,
                       force_non_negative = T, robust = T) {

  # check pos and neg having same channels
  if (!is.null(neg)) {
    if (!all(colnames(pos) %in% colnames(neg))) {
      stop("Channel names in pos and neg do not match.")
    }
  }

  # if no channel list is given, use all columns of pos
  channels_used <- ChannelCheck(pos, channels)

  ch.fix <- make.names(channels_used)

  purrr::map_dfr(unique(pos$file), function(f) {

    # combine positive and negative events if necessary
    if (!is.null(neg)) {
      dat_tmp <- dplyr::bind_rows(dplyr::filter(.data$pos, file == f),
                                  neg) %>%
        dplyr::select(tidyselect::all_of(channels_used)) %>%
        magrittr::set_colnames(ch.fix)
    } else {
      dat_tmp <- pos %>%
        dplyr::filter(.data$file == f) %>%
        dplyr::select(tidyselect::all_of(channels_used)) %>%
        magrittr::set_colnames(ch.fix)
    }

    # find channels with strongest signal (reference set to 1)
    indices <- dat_tmp %>%
      dplyr::summarise_all(function(x) stats::quantile(x, 0.99)) %>%
      unlist %>%
      which.max

    ch_max <- ch.fix[indices]

    # calculate spillover
    res <- purrr::map_dbl(ch.fix, function(ch) {
      if (ch != ch_max) {
        if (robust) {
          MASS::rlm(stats::as.formula(stringr::str_c(ch, "~", ch_max)),
                    dat_tmp, maxit=100)$coefficients[2]
        } else {
          stats::lm(stats::as.formula(stringr::str_c(ch, "~", ch_max)),
                    data=dat_tmp)$coefficients[2]
        }
      } else {
        1
      }
    })

    # set negative values to 0
    if (any(res < 0)) {
      warning(stringr::str_c("Warning: ", sum(res < 0),
                             " spectra values below 0 for ", f))
    }
    if (force_non_negative) {
      res[res < 0] <- 0
    }

    # output
    names(res) <- channels_used
    res
  }) %>%
    dplyr::mutate(file = unique(pos$file))
}


# Single stain variability plots ------------------------------------------


#' Draw a sample and correct for background
#'
#' Used in SingleStainVariability
#'
#' @param pos dataframe containing raw events
#' @param channels optional, channels to use
#' @param neg_median dataframe containing the autofluorescence to subtract
#' @param sample_size size of the sample to draw from input data
#'
#' @importFrom magrittr %>%
#'
#' @return tibble with the corrected sample
SampleCorr <- function(pos, channels = NULL,
                       neg_median = NULL, sample_size = 1000) {

  # if no channel list is given, use all columns of pos
  channels_used <- ChannelCheck(pos, channels)

  tmp1 <- pos %>%
    dplyr::group_by(file) %>%
    dplyr::slice_sample(n=sample_size) %>%
    dplyr::ungroup()

  if (!is.null(neg_median)) {
    tmp2 <- tmp1 %>%
      dplyr::select(tidyselect::all_of(channels_used)) %>%
      as.matrix %>%
      sweep(2, unlist(dplyr::select(neg_median, tidyselect::all_of(channels_used))))
  } else {
    tmp2 <- tmp1 %>%
      dplyr::select(tidyselect::all_of(channels_used)) %>%
      as.matrix
  }

  # normalize to 1
  tmp2 <- tmp2 / apply(tmp2, 1, max)

  tmp2 %>%
    tibble::as_tibble() %>%
    dplyr::mutate(file = tmp1$file)
}


#' Calculate cosine similarity matrix
#'
#' https://stats.stackexchange.com/questions/31565/compute-a-cosine-dissimilarity-matrix-in-r
#'
#' @param m Matrix
#'
#' @return similarity matrix
CosSim <- function(m) {
  sim <- m / sqrt(rowSums(m * m))
  sim %*% t(sim)
}


#' Plot the variability of the single stains
#'
#' @param input_list Input list of lists in the form of list(list(pos_raw, neg_median), ...). Each sublist must contain pos_raw, the dataframe of the single stain raw data and neg_median, the channel medians of the corresponding negative control.
#' @param spectra Tibble of the calculated reference spectra (output of MedianSpectra or RegSpectra).
#' @param channels Optional, channels to be included in the calculation in calculation
#' @param sample_size Number of events to sample from each pos_raw
#' @param transformation optional asinh transform of the spectra. Default: NA. If you want to transform the normalized spectra, pass the cofactor of the asinh transformation. Note: The spectra are max normalized to 1 so your usual cofactors will not work here. Try cofactors around 0.001.
#'
#' @return Plot of single stain variability, faceted by single stain. The red line shows the calculated reference spectra, the blue lines the (autofluorescence subtracted) single stain events. The blue shade indicates the average similarity of that event to all other events (light blue: more similar, dark blue: less dissimilar).
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
SingleStainVariability <- function(input_list,
                                   spectra,
                                   channels = NULL,
                                   sample_size = 1000,
                                   transformation = NULL) {

  channels_used <- ChannelCheck(spectra, channels)

  # sample and normalize data
  dat_sample_norm <- purrr::map_dfr(input_list, function(li) {
    SampleCorr(li[[1]],
               channels_used,
               if (length(li) == 1) {NULL} else {li[[2]]},
               sample_size)
  })

  # get means cosine similarities
  sim_mean <- purrr::map(unique(dat_sample_norm$file), function(r) {
    tmp <- dat_sample_norm %>%
      dplyr::filter(file == r) %>%
      dplyr::select(-file) %>%
      as.matrix %>%
      CosSim
    # mean cosine similarity
    (rowSums(tmp) - 1) / (ncol(tmp)-1)
  }) %>%
    unlist

  # prepare spectra
  spectra_prep <- spectra %>%
    dplyr::filter(file %in% unique(dat_sample_norm$file)) %>%
    dplyr::select(tidyselect::all_of(c(channels_used, "file"))) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    tidyr::pivot_longer(cols = tidyselect::all_of(channels_used),
                        names_to = "channel",
                        values_to = "signal") %>%
    dplyr::mutate(channel = factor(.data$channel, channels_used),
                  signal = if (!is.null(transformation) & is.numeric(transformation)) {
                    asinh(.data$signal/transformation)
                  } else {
                    .data$signal
                  })

  # generate plots
  dat_sample_norm %>%
    dplyr::mutate(mSim = sim_mean,
                  id = 1:dplyr::n()) %>%
    tidyr::pivot_longer(cols = tidyselect::all_of(channels_used),
                        names_to = "channel",
                        values_to = "signal") %>%
    dplyr::mutate(channel = factor(.data$channel, channels_used),
                  signal = if (!is.null(transformation) & is.numeric(transformation)) {
                    asinh(.data$signal/transformation)
                  } else {
                    .data$signal
                  }) %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$channel,
                                 y = .data$signal,
                                 group = .data$id,
                                 color = .data$mSim)) +
    ggplot2::geom_line(alpha = 0.5) +
    ggplot2::geom_line(data = spectra_prep, color = "red") +
    ggplot2::facet_wrap(ggplot2::vars(file), ncol = 3, scales = "free_y") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                       vjust = 0.5,
                                                       hjust = 1))
}



# Spectra plots -----------------------------------------------------------

#' Plot Spectra
#'
#' @param dat tibble, output of GetData
#' @param channels character vector, channel selection to include in the plot (default: NULL, i.e. all)
#' @param bins numeric, number of bins to use in the channel heatmaps (default: 50)
#' @param transformation numeric, cofactor used for the asinh transformation (default: 1000)
#' @param ncol integer, number of columns of the facet plot (default: NULL, i.e. automatic)
#'
#' @return ggplot object
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
PlotSpectra <- function(dat,
                        channels = NULL,
                        bins = 50,
                        transformation = 1000,
                        ncol = NULL) {

  channels_used <- ChannelCheck(dat, channels)

  files <- unique(dat$file)

  dat <- dat %>%
    dplyr::mutate_at(dplyr::all_of(channels_used), function(x) asinh(x/transformation))

  dat_binned <- purrr::map_dfr(files, function(f) {
    dat_tmp <- dplyr::filter(dat, .data$file == f) %>%
      dplyr::select(dplyr::all_of(channels_used))
    purrr::map_dfr(channels_used, function(ch) {
      bin_counts <- table(cut(dat_tmp[[ch]], breaks = bins))

      names(bin_counts) %>%
        stringr::str_split(",") %>%
        purrr::map_dfr(function(x) {
          borders <- as.numeric(stringr::str_extract(x, "((?<=\\()(.*))|((.*)(?=\\]))"))
          tibble::tibble(height = diff(borders), mean = mean(borders))
        }) %>%
        dplyr::mutate(freq = as.numeric(bin_counts)/sum(bin_counts),
                      channel = ch,
                      file = f)
    })
  }) %>%
    dplyr::mutate(channel = factor(channel, channels_used))

  dat_binned %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$channel,
                                 y = .data$mean,
                                 height = .data$height,
                                 fill = .data$freq)) +
    ggplot2::geom_tile(width = 2, show.legend = F) +
    ggplot2::scale_fill_gradientn(colours = c("#FFFFFF",
                                              rev(grDevices::rainbow(100,
                                                                     start = 0,
                                                                     end = 0.6)))) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                       vjust = 0.5,
                                                       hjust = 1)) +
    ggplot2::scale_y_continuous(breaks = asinh(c(-10^seq(4,3),0,10^seq(3,7))/transformation),
                                labels = c(-10^seq(4,3),0,10^seq(3,7))) +
    ggplot2::labs(x = "Channel", y = "Signal") +
    ggplot2::facet_wrap(ggplot2::vars(file), ncol = ncol)
}


#' Plot reference spectra
#'
#' @param spectra spectra dataframe as produced by MedianSpectra or RegSpectra
#' @param transformation optional asinh transform of the spectra. Default: NA. If you want to transform the normalized spectra, pass the cofactor of the asinh transformation. Note: The spectra are max normalized to 1 so your usual cofactors will not work here. Try cofactors around 0.001.
#'
#' @return ggplot object
#'
#' @importFrom rlang .data
#'
#' @export
PlotRefSpectra <- function(spectra, transformation = NA) {
  channels <- colnames(spectra)

  spectra %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    tidyr::pivot_longer(cols = tidyr::all_of(channels[channels != "file"]),
                        names_to = "channel", values_to = "signal") %>%
    dplyr::mutate(channel = factor(.data$channel, channels),
                  signal = if (!is.null(transformation) & is.numeric(transformation)) {
                    asinh(.data$signal/transformation)
                  } else {
                    .data$signal
                  }) %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$channel,
                                 y = .data$signal,
                                 color = .data$file,
                                 group = .data$id)) +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
}


