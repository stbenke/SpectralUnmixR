

# Unmixing functions ------------------------------------------------------

#' Wrapper to apply unmixing function to each row in data matrix
#'
#' @param dat Dataframe to unmix
#' @param UnmixFun unmixing function to apply
#' @param ... params passed to unmixing function
#'
#' @return abundance matrix
UnmixMatrix <- function(dat, UnmixFun, ...) {
  t(apply(dat, 1, match.fun(UnmixFun), ...))
}


#' Ordinary least squares unmixing
#'
#' @param x event row
#' @param U (M^TM)*M^T
#'
#' @return abundance matrix
UnmixFunc.ols <- function(x, U) {
  U %*% x
}
#' Ordinary least squares unmixing
#'
#' https://en.wikipedia.org/wiki/Ordinary_least_squares
#'
#' @param dat dataframe, data to be unmixed
#' @param M unmixing matrix
#' @param bg vector of constant background to be subtracted from each event
#'
#' @return dataframe with component abundances
Unmix.ols <- function(dat, M, bg = NULL) {
  # calculate unmixing matrix to speed things up
  U <- solve(t(M) %*% M) %*% t(M)
  if (is.null(bg)) {
    tmp <- dat
  } else {
    tmp <- sweep(dat, 2, bg)
  }
  UnmixMatrix(tmp, UnmixFunc.ols, U)
}


#' Weighted least squares
#'
#' @param x event row
#' @param M unmixing matrix
#' @param bg vector of constant background to be subtracted from each event
#'
#' @return vector with component abundances
UnmixFunc.wls.bg <- function(x, M, bg) {
  solve(t(M) %*% diag(1/abs(x)) %*% M) %*% t(M) %*% diag(1/abs(x)) %*% (x-bg)
}
#' Weighted least squares
#'
#' @param x event row
#' @param M unmixing matrix
#'
#' @return vector with component abundances
UnmixFunc.wls.nobg <- function(x, M) {
  solve(t(M) %*% diag(1/abs(x)) %*% M) %*% t(M) %*% diag(1/abs(x)) %*% x
}
#' Weighted least squares
#'
#' https://en.wikipedia.org/wiki/Weighted_least_squares
#' question of how to treat negative values and whether to subtract background for the weights.
#'
#' @param dat dataframe, data to be unmixed
#' @param M unmixing matrix
#' @param bg vector of constant background to be subtracted from each event
#'
#' @return dataframe with component abundances
Unmix.wls <- function(dat, M, bg = NULL) {
  if (is.null(bg)) {
    UnmixMatrix(dat, UnmixFunc.wls.nobg, M)
  } else {
    UnmixMatrix(dat, UnmixFunc.wls.bg, M, bg)
  }
}


#' Non-negative least squares unmixing
#'
#' @param x event row
#' @param M unmixing matrix
#'
#' @return vector with component abundances
UnmixFunc.nnls <- function(x, M) {
  nnls::nnls(M, x)$x
}
#' Non-negative least squares unmixing
#'
#' @param dat dataframe, data to be unmixed
#' @param M unmixing matrix
#' @param bg vector of constant background to be subtracted from each event
#'
#' @return dataframe with component abundances
Unmix.nnls <- function(dat, M, bg = NULL) {
  if (is.null(bg)) {
    tmp <- dat
  } else {
    tmp <- sweep(dat, 2, bg)
  }
  UnmixMatrix(tmp, UnmixFunc.nnls, M)
}


#' Weighted NNLS with background
#'
#' @param x event row
#' @param M unmixing matrix
#' @param bg vector of constant background to be subtracted from each event
#'
#' @return vector with component abundances
UnmixFunc.wnnls.bg <- function(x, M, bg) {
  nnls::nnls(diag(sqrt(1/abs(x))) %*% M, sqrt(1/abs(x)) * (x-bg))$x
}
#' Weighted NNLS without background
#'
#' @param x event row
#' @param M unmixing matrix
#'
#' @return vector with component abundances
UnmixFunc.wnnls.nobg <- function(x, M) {
  nnls::nnls(diag(sqrt(1/abs(x))) %*% M, sqrt(1/abs(x)) * x)$x
}
#' Weighted NNLS unmixing
#'
#' https://stackoverflow.com/questions/47888996/weighted-nonnegative-least-squares-in-r
#'
#' @param dat dataframe, data to be unmixed
#' @param M unmixing matrix
#' @param bg vector of constant background to be subtracted from each event
#'
#' @return dataframe with component abundances
Unmix.wnnls <- function(dat, M, bg = NULL) {
  if (is.null(bg)) {
    UnmixMatrix(dat, UnmixFunc.wnnls.nobg, M)
  } else {
    UnmixMatrix(dat, UnmixFunc.wnnls.bg, M, bg)
  }
}


#' MAPE unmixing with backround
#'
#' Use in wrapper
#' question of how to treat negative values and whether to subtract background for the weights (probably not).
#'
#' @param x event row
#' @param M unmixing matrix
#' @param bg vector of constant background to be subtracted from each event
#'
#' @return vector with component abundances
UnmixFunc.mape.bg <- function(x, M, bg) {
  solve(t(M) %*% diag(1/abs(x^2)) %*% M) %*% t(M) %*% diag(1/abs(x^2)) %*% (x-bg)
}
#' MAPE unmixing without backround
#'
#' Use in wrapper
#' question of how to treat negative values and whether to subtract background for the weights (probably not).
#'
#' @param x event row
#' @param M unmixing matrix
#'
#' @return vector with component abundances
UnmixFunc.mape.nobg <- function(x, M) {
  solve(t(M) %*% diag(1/abs(x^2)) %*% M) %*% t(M) %*% diag(1/abs(x^2)) %*% x
}

#' MAPE unmixing wrapper
#'
#' @param dat dataframe, data to be unmixed
#' @param M unmixing matrix
#' @param bg vector of constant background to be subtracted from each event (default: NULL)
#'
#' @return dataframe with component abundances
Unmix.mape <- function(dat, M, bg = NULL) {
  if (is.null(bg)) {
    UnmixMatrix(dat, UnmixFunc.mape.nobg, M)
  } else {
    UnmixMatrix(dat, UnmixFunc.mape.bg, M, bg)
  }
}


#' Wrapper for all defined unmixing methods
#'
#' @param dat dataframe, data to be unmixed
#' @param M unmixing matrix
#' @param bg vector of constant background to be subtracted from each event (default: NULL)
#' @param unmixFunc string, unmixing method. OLS, WLS, NNLS, WNNLS or MAPE. (default: "OLS")
#'
#' @return matrix with component abundances
Unmix.wrapper <- function(dat, M, bg = NULL, unmixFunc = "OLS") {

  # Test if channels of dat and M match
  if (ncol(dat) != nrow(M)) {
    stop("Number of columns of the input data must match the number of rows of the spectra matrix.")
  } else if (!any(colnames(dat) == rownames(M))) {
    stop("Column names of input data must match the rownames of the spectra matrix in name and order.")
  }

  # Apply chosen unmixing function
  if (unmixFunc == "OLS") {
    dat %>%
      Unmix.ols(M, bg) %>%
      magrittr::set_colnames(colnames(M))
  } else if (unmixFunc == "NNLS") {
    dat %>%
      Unmix.nnls(M, bg) %>%
      magrittr::set_colnames(colnames(M))
  } else if (unmixFunc == "WLS") {
    # put in offset for 0 values to make weighted unmixing possible
    dat[dat==0] <- 0.05

    dat %>%
      Unmix.wls(M, bg) %>%
      magrittr::set_colnames(colnames(M))
  } else if (unmixFunc == "WNNLS") {
    # put in offset for 0 values to make weighted unmixing possible
    dat[dat==0] <- 0.05

    dat %>%
      Unmix.wnnls(M, bg) %>%
      magrittr::set_colnames(colnames(M))
  } else if (unmixFunc == "MAPE") {
    # put in offset for 0 values to make weighted unmixing possible
    dat[dat==0] <- 0.05

    dat %>%
      Unmix.mape(M, bg) %>%
      magrittr::set_colnames(colnames(M))
  } else {
    stop(stringr::str_c(unmixFunc, " is not a known unmixing function."))
  }
}


# Unmixing wrappers -------------------------------------------------------

#' Unmix a dataset
#'
#' Use inside Unmix wrapper.
#'
#' @param data dataframe to unmix
#' @param M unmixing matrix
#' @param bg dataframe containing background (optional)
#' @param unmixing unmixing function to use (string)
#' @param write_FCS write output into FCS files (default: FALSE)
#' @param FCS_suffix suffix added to the FCS files (default: "unmixed")
#' @param FCS_dir directory to export the FCS files to (default: ".")
#' @param fullOutput Add the reconstruction and unmixing matrix to the output (default: FALSE). Needed for the interactive visualization
#'
#' @return list(unmixed, reconstruction, M)
#' @export
#'
UnmixFile <- function(data, M, bg = NULL,
                      unmixing = "OLS",
                      write_FCS = FALSE, FCS_suffix = "unmixed",  FCS_dir = ".",
                      fullOutput = FALSE) {

  if (!is.data.frame(data)) {
    stop("Input data must be passed as dataframe or tibble.")
  }
  if (!is.matrix(M)) {
    stop("Spectral matrix must be passed as matrix.")
  }

  # prep data
  measurement <- data %>%
    dplyr::select(dplyr::all_of(rownames(M))) %>%
    as.matrix

  # unmixing
  abundances <- Unmix.wrapper(measurement, M, bg, unmixing)

  # reconstruction
  if (is.null(bg)) {
    reconstruction <- abundances %*% t(M)
  } else {
    reconstruction <- sweep(abundances %*% t(M), 2, bg, FUN = "+")
  }

  # get unmixing error
  # todo: implement other errors
  total_squared_error <- rowSums((measurement - reconstruction)^2)

  # add not unmixed columns and error
  abundances_out <- abundances %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(error = total_squared_error) %>%
    dplyr::bind_cols(dplyr::select(data, -c(rownames(M), "file"))) # add scatter

  # write FCS files
  if (write_FCS) {
    dir.create(FCS_dir, showWarnings = FALSE)

    # if multiple files in input, split into different FCS files
    files <- unique(data$file)

    for (f in files) {
      abundances_out %>%
        as.matrix %>%
        flowCore::flowFrame() %>%
        flowCore::write.FCS(file.path(FCS_dir,
                                      stringr::str_c(f,
                                                     "_", FCS_suffix,
                                                     "_", unmixing, ".fcs"))
        )
    }
  }

  # return output
  list(
    unmixed = abundances_out %>%
      dplyr::bind_cols(dplyr::select(data, "file")),
    reconstruction = if (fullOutput) reconstruction else NA,
    M = if (fullOutput) M else NA
  )
}



#' Unmixing with multiple AF
#'
#' Use multiple unmixing matrices, each with different AF and take best result for each cell.
#' AF columns in M list must have unique names (e.g. AF1, AF2, ...)!
#'
#' Intended to be used inside Unmix wrapper.
#'
#' @param data dataframe to unmix
#' @param M list of unmixing matrices
#' @param bg dataframe containing background (optional)
#' @param unmixing unmixing function to use (string)
#' @param write_FCS write output into FCS files (default: FALSE)
#' @param FCS_suffix suffix added to the FCS files (default: "unmixed")
#' @param FCS_dir directory to export the FCS files to (default: ".")
#' @param fullOutput Add the reconstruction and unmixing matrix to the output (default: FALSE). Needed for the interactive visualization
#'
#' @return list(unmixed, reconstruction, M)
#'
#' @export
#'
UnmixFile_MultiAF <- function(data, M, bg = NULL,
                              unmixing = "OLS",
                              write_FCS = FALSE, FCS_suffix = "unmixed",  FCS_dir = ".",
                              fullOutput = FALSE) {

  if (!is.data.frame(data)) {
    stop("Input data must be passed as dataframe or tibble.")
  }
  if (!is.list(M)) {
    stop("The collection of spectral matrices must be passed as list.")
  } else if (any(!purrr::map_lgl(M, is.matrix))) {
    stop("Spectral matrices must be passed as matrix.")
  }

  # autofluorescence must be last column in each M

  # prep data
  # check if all matrices have same rownames?
  measurement <- data %>%
    dplyr::select(dplyr::all_of(rownames(M[[1]]))) %>%
    as.matrix

  # unmixing
  # currently only supports one type of background
  abundances <- purrr::map(M, ~Unmix.wrapper(measurement, ., bg, unmixing))

  # reconstruction
  if (is.null(bg)) {
    reconstruction <- purrr::map2(abundances, M, function(a, m) {
      a %*% t(m)
    })
  } else {
    reconstruction <- purrr::map2(abundances, M, function(a, m) {
      sweep(a %*% t(m), 2, bg, FUN = "+")
    })
  }

  # unmixing error
  total_squared_error <- do.call(cbind, purrr::map(reconstruction, function(r) {
    rowSums((measurement - r)^2)
  }))

  # best unmixing for each event
  ind <- apply(total_squared_error, 1, which.min)

  abundances_out <- purrr::map_dfr(seq_along(M), function(i) {
    tmp <- which(ind == i)
    out <- abundances[[i]][tmp,] %>%
      dplyr::as_tibble() %>%
      dplyr::bind_cols(dplyr::select(data[tmp,],
                                     !c(colnames(measurement), "file"))) %>%
      dplyr::mutate(error = total_squared_error[tmp,i],
                    ind = which(ind == i))
    for (j in seq_along(M)) {
      if (j != i) {
        out[[colnames(M[[j]])[ncol(M[[j]])]]] = 0
      }
    }
    out
  }) %>%
    dplyr::arrange(ind) %>% # put into original order again
    dplyr::select(-ind)

  # write FCS
  if (write_FCS) {
    dir.create(FCS_dir, showWarnings = FALSE)

    # if multiple files in input, split into different FCS files
    files <- unique(data$file)

    for (f in files) {
      abundances_out %>%
        as.matrix %>%
        flowCore::flowFrame() %>%
        flowCore::write.FCS(file.path(FCS_dir,
                                      stringr::str_c(f,
                                                     "_", FCS_suffix,
                                                     "_", unmixing, ".fcs"))
        )
    }
  }

  # return output
  if (fullOutput) {
    M_out <- do.call(cbind, purrr::map(seq_along(M), function(i) {
      if (i == 1) {
        M[[i]]
      } else {
        M[[i]][,ncol(M[[i]])]
      }
    })) %>%
      magrittr::set_colnames(c(colnames(M[[1]]),
                               purrr::map_chr(2:length(M), ~colnames(M[[.]])[ncol(M[[.]])])))

    if (is.null(bg)) {
      reconstruction_out <- as.matrix(dplyr::select(abundances_out, colnames(M_out))) %*% t(M_out)
    } else {
      reconstruction_out <- sweep(as.matrix(dplyr::select(abundances_out, colnames(M_out))) %*% t(M_out), 2, bg, FUN = "+")
    }

    list(
      unmixed = abundances_out %>%
        dplyr::bind_cols(dplyr::select(data, "file")),
      reconstruction = reconstruction_out,
      M = M_out
    )
  } else {
    list(
      unmixed = abundances_out %>%
        dplyr::bind_cols(dplyr::select(data, "file")),
      reconstruction = NA,
      M = NA
    )
  }
}



#' Spectral unmixing
#'
#' Spectral unmixing using OLS, WLS, NNLS, WNNLS or MAPE with one or multiple autofluorescence populations.
#' For all weighted algorithms: For the weight of negative values, the absolute is taken and for 0 values, an offset of 0.05 is added.
#'
#' @param data Raw data to be unmixed. Can be passed as dataframe, tibble, flowframe or flowSet.
#' @param M Matrix containing the component spectra or list of such matrices.
#' @param bg Constant background or autofluorescence to be subtracted of events before unmixing (dataframe). Default: NULL
#' @param unmixing Unmixing method to use passed as string. Currently supported methods are: ordinary least squares ("OLS"),
#'  weighted least squares ("WLS"), non-negative least squares ("NNLS"), weighted non-negative least squares ("WNNLS") and mean
#'  absolute percentage error ("MAPE").
#' @param write_FCS Write unmixed file to FCS file (TRUE / FALSE (default))?
#' @param FCS_suffix Suffix to be added to the unmixed FCS file name. Default: "unmixed".
#' @param FCS_dir Directory to write the unmixed FCS file to. Default: ".".
#' @param fullOutput If TRUE, returns a list containing the reconstructed data and unmixing matrix in addition to the abundances.
#' Necessary for the interactive visualization.
#'
#' @return List containing the abundances (unmixed). If fullOutput = TRUE, the list also contains the reconstructed data (reconstruction),
#' and the unmixing matrix. If write_FCS = TRUE, the function also writes an FCS file containing the abundances.
#' @export
#' @importFrom magrittr %>%
#'
Unmix <- function(data, M, bg = NULL,
                  unmixing = "OLS",
                  write_FCS = FALSE, FCS_suffix = "unmixed",  FCS_dir = ".",
                  fullOutput = FALSE) {

  if (!is.null(bg)) {
    .bg <- unlist(bg[,rownames(M)])
  } else {
    .bg <- NULL
  }

  if (is.matrix(M)) {
    if (is.data.frame(data)) {
      UnmixFile(data, M, .bg,
                unmixing,
                write_FCS, FCS_suffix, FCS_dir,
                fullOutput)
    } else if (methods::is(data, "flowFrame")) {
      flowCore::exprs(data) %>%
        tibble::as_tibble() %>%
        dplyr::mutate(file = flowCore::keyword(data, "$FIL")[[1]]) %>%
        UnmixFile(M, .bg,
                  unmixing,
                  write_FCS, FCS_suffix, FCS_dir,
                  fullOutput)
    } else if (methods::is(data, "flowSet")) {
      purrr::map(flowCore::flowSet_to_list(data), function(dat) {
        flowCore::exprs(dat) %>%
          tibble::as_tibble() %>%
          dplyr::mutate(file = flowCore::keyword(data, "$FIL")[[1]]) %>%
          UnmixFile(M, .bg,
                    unmixing,
                    write_FCS, FCS_suffix, FCS_dir,
                    fullOutput)
      })
    } else {
      stop("Input data must be passed as dataframe, tibble, flowframe or flowset.")
    }
  } else if (is.list(M)) {
    if (all(purrr::map_lgl(M, is.matrix))) {
      if (is.data.frame(data)) {
        UnmixFile_MultiAF(data, M, .bg,
                          unmixing,
                          write_FCS, FCS_suffix, FCS_dir,
                          fullOutput)
      } else if (methods::is(data, "flowFrame")) {
        flowCore::exprs(data) %>%
          tibble::as_tibble() %>%
          dplyr::mutate(file = flowCore::keyword(data, "$FIL")[[1]]) %>%
          UnmixFile_MultiAF(M, .bg,
                            unmixing,
                            write_FCS, FCS_suffix, FCS_dir,
                            fullOutput)
      } else if (methods::is(data, "flowSet")) {
        purrr::map(flowCore::flowSet_to_list(data), function(dat) {
          flowCore::exprs(dat) %>%
            tibble::as_tibble() %>%
            dplyr::mutate(file = flowCore::keyword(data, "$FIL")[[1]]) %>%
            UnmixFile_MultiAF(M, .bg,
                              unmixing,
                              write_FCS, FCS_suffix, FCS_dir,
                              fullOutput)
        })
      } else {
        stop("Input data must be passed as dataframe, tibble, flowframe or flowset.")
      }
    } else {
      stop("Spectral matrix must be passed as matrix or list of matrices.")
    }
  } else {
    stop("Spectral matrix must be passed as matrix or list of matrices.")
  }
}

#' Create unmixing matrix
#'
#' @param dye_spectra Dataframe containing the dye spectra in the rows.
#' @param channels Optional, vector containing the channel names to be used. Default (NA) uses all channels.
#'
#' @return Unmixing matrix
#' @export
#'
MakeM <- function(dye_spectra, channels = NA) {
  if (is.na(channels)) {
    channels_used <- colnames(dye_spectra)[colnames(dye_spectra) != "file"]
  } else {
    channels_used <- channels
  }

  dye_spectra %>%
    dplyr::select(dplyr::all_of(channels_used)) %>%
    as.matrix %>%
    t %>%
    magrittr::set_colnames(dye_spectra$file)
}

#' Create list of unmixing matrices
#'
#' Creates a list of unmixing matrices to be used when there are multiple autofluorescence populations and each event should be unmixed using only one autofluorescence population. Each unmixing matrix contains the dye spectra plus one autofluorescence spectra.
#'
#' @param AF Dataframe containing the autofluorescence spectra in the rows.
#' @param dye_spectra Dataframe containing the dye spectra in the rows.
#' @param channels Optional, vector containing the channel names to be used. Default (NA) uses all channels.
#'
#' @return List of unmixing matrices to used in Unmixing
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
MakeMList <- function(AF, dye_spectra, channels = NULL) {
  if (is.null(channels)) {
    channels_used <- colnames(dye_spectra)[colnames(dye_spectra) != "file"]
  } else {
    channels_used <- channels
  }

  purrr::map(1:nrow(AF), function(i) {
    dplyr::bind_rows(dye_spectra, dplyr::mutate(AF[i,], file = as.character(.data$file))) %>%
      dplyr::select(dplyr::all_of(channels_used)) %>%
      as.matrix %>%
      t %>%
      magrittr::set_colnames(c(dye_spectra$file, stringr::str_c("AF", i)))
  })
}


