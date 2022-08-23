#' Visualize unmixing results
#'
#' The top plot shows the abundances resulting from the unmixing in a 2D dot plot. The user can choose which parameters to display and how to transform them. The color indicates the unmixing error.
#' By clicking on a data point, details for that point are shown in the two plots below. In the middle plot, the fitted spectral contributions for the selected data point are shown. The bottom plot shows the position of the selected data point in FSC/SSC.
#'
#' @param unmixing unmixing result (output of Unmix with full_output = TRUE)
#' @param original original raw data (output of GetData)
#' @param n_sample number of data points to sample for plots
#'
#' @return opens interactive visualization
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
UnmixingViz <- function(unmixing, original,
                        n_sample = 1000) {

  sel <- sample(1:nrow(unmixing$unmixed), n_sample)
  dat_plot <- unmixing$unmixed %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    dplyr::slice(sel)
  M <- unmixing$M
  abundances <- as.matrix(dplyr::select(unmixing$unmixed,
                                        dplyr::all_of(colnames(M))))
  reconstruction <- unmixing$reconstruction
  measurement <- dplyr::select(original,
                               dplyr::all_of(rownames(M)))
  markers <- colnames(M)


  ui <- miniUI::miniPage(
    miniUI::gadgetTitleBar("My Gadget"),
    miniUI::miniContentPanel(
      shiny::fluidRow(
        shiny::div(
          shiny::div(style="display: inline-block; vertical-align: top;",
                     shiny::uiOutput("marker_x")),
          shiny::div(style="display: inline-block; vertical-align: top;",
                     shiny::radioButtons("marker_x_scaling",
                                         "x-axis scaling",
                                         choices = list("asinh", "linear"),
                                         selected = "asinh",
                                         inline = T)),
          shiny::div(style="display: inline-block; vertical-align: top;",
                     shiny::numericInput("marker_x_cofactor",
                                         label = "x-axis scaling cofactor",
                                         value = 1000,
                                         min = 1))
        ),
        shiny::div(
          shiny::div(style="display: inline-block; vertical-align: top;",
                     shiny::uiOutput("marker_y")),
          shiny::div(style="display: inline-block; vertical-align: top;",
                     shiny::radioButtons("marker_y_scaling",
                                         "y-axis scaling",
                                         choices = list("asinh", "linear"),
                                         selected = "asinh",
                                         inline = T)),
          shiny::div(style="display: inline-block; vertical-align: top;",
                     shiny::numericInput("marker_y_cofactor",
                                         label = "y-axis scaling cofactor",
                                         value = 1000,
                                         min = 1))
        ),
        shiny::div(
          shiny::div(style="display: inline-block; vertical-align: top;",
                     shiny::uiOutput("color_error")),
          shiny::div(style="display: inline-block; vertical-align: top;",
                     shiny::radioButtons("error_scaling",
                                         "error (color) scaling",
                                         choices = list("log", "linear"),
                                         selected = "log",
                                         inline = T))
        )
      ),
      shiny::fluidRow(
        plotly::plotlyOutput("marker_plot_scatter")
      ),
      shiny::div(
        shiny::div(style="display: inline-block; vertical-align: top;",
                   shiny::radioButtons("signal_scaling",
                                       "signal scaling",
                                       choices = list("asinh", "linear"),
                                       selected = "asinh",
                                       inline = T)),
        shiny::div(style="display: inline-block; vertical-align: top;",
                   shiny::numericInput("signal_cofactor",
                                       label = "signal scaling cofactor",
                                       value = 1000,
                                       min = 1))
      ),
      shiny::fluidRow(
        plotly::plotlyOutput("spectra_plot")
      ),
      shiny::fluidRow(
        plotly::plotlyOutput("scatter_plot")
      )
    )
  )

  server <- function(input, output, session) {
    # Define reactive expressions, outputs, etc.

    output$marker_x <- shiny::renderUI({
      shiny::selectInput("marker_x_selection",
                         label = "x-axis",
                         choices = markers,
                         selected = markers[1])
    })
    output$marker_y <- shiny::renderUI({
      shiny::selectInput("marker_y_selection",
                         label = "y-axis",
                         choices = markers,
                         selected = markers[1])
    })

    output$marker_plot_scatter <- plotly::renderPlotly({
      if (!is.null(input$marker_x_selection)) {
        if (!is.null(input$marker_y_selection)) {
          if (input$marker_x_scaling == "asinh") {
            x_vals <- asinh(dat_plot[[input$marker_x_selection]] / input$marker_x_cofactor)
          } else {
            x_vals <- dat_plot[[input$marker_x_selection]]
          }
          if (input$marker_y_scaling == "asinh") {
            y_vals <- asinh(dat_plot[[input$marker_y_selection]] / input$marker_y_cofactor)
          } else {
            y_vals <- dat_plot[[input$marker_y_selection]]
          }
          if (input$error_scaling == "log") {
            color_vals <- log(dat_plot$error)
          } else {
            color_vals <- dat_plot$error
          }
          dat_plot %>%
            plotly::plot_ly(x = x_vals, y = y_vals,
                            text = dat_plot$id, color = color_vals,
                            customdata = dat_plot$id,
                            type = "scatter", mode = "markers",
                            showlegend = T, hoverinfo = "text",
                            source = "marker_scatter")
        }
      }
    })

    output$scatter_plot <- plotly::renderPlotly({
      suppressWarnings(
        sel <- as.numeric(plotly::event_data("plotly_click",
                                             source = "marker_scatter")$customdata)
      )
      if (length(sel) == 0) {
        sel <- 1
      }
      dat_plot %>%
        plotly::plot_ly(x = ~ `FSC-A`, y = ~ `SSC-A`,
                        type = "histogram2dcontour") %>%
        plotly::add_trace(data = dplyr::filter(dat_plot, .data$id %in% sel),
                          x = ~ `FSC-A`, y = ~ `SSC-A`,
                          type = "scatter", mode = "markers",
                          text = ~ id, hoverinfo = "text")
    })

    output$spectra_plot <- plotly::renderPlotly({
      suppressWarnings(
        sel <- as.numeric(plotly::event_data("plotly_click",
                                             source = "marker_scatter")$customdata)
      )
      if (length(sel) == 0) {
        sel <- 1
      }
      components <- do.call(rbind, c(
        purrr::map(seq_along(abundances[sel,]), function(j) {
          t(M)[j,] * abundances[sel,j]
        }),
        list(measurement[sel,], reconstruction[sel,])
      )
      ) %>%
        dplyr::as_tibble() %>%
        dplyr::mutate(component = c(colnames(M),
                                    "original",
                                    "reconstruction")) %>%
        tidyr::pivot_longer(cols = rownames(M),
                            names_to = "channel",
                            values_to = "signal") %>%
        dplyr::mutate(channel = factor(.data$channel, rownames(M)))

      if (input$signal_scaling == "asinh") {
        components <- components %>%
          dplyr::mutate(signal = asinh(.data$signal/input$signal_cofactor))
      }

      plotly::plot_ly(data = components,
                      x = ~ channel, y = ~ signal,
                      color = ~ component,
                      text = ~ component, hoverinfo = "text",
                      type = "scatter", mode = "lines")
    })

    # When the Done button is clicked, return a value
    shiny::observeEvent(input$done, {
      returnValue <- NULL
      shiny::stopApp(returnValue)
    })
  }

  shiny::runGadget(ui, server)
}



#' Compare two unmixing results side by side
#'
#' The top plot shows the abundances resulting from the unmixing in a 2D dot plot. The user can choose which parameters to display and how to transform them. The color indicates the unmixing error.
#' By clicking on a data point, details for that point are shown in the plots below. In the middle plots, the fitted spectral contributions for the selected data point are shown. The bottom plot shows the position of the selected data point in FSC/SSC.
#'
#' @param unmixing1 first unmixing result (output of Unmix with full_output = TRUE)
#' @param unmixing2 second unmixing result (output of Unmix with full_output = TRUE)
#' @param original original raw data (output of GetData)
#' @param n_sample number of data points to sample for plots
#'
#' @return opens interactive visualization
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
UnmixingCompViz <- function(unmixing1, unmixing2,
                            original, n_sample = 1000) {

  sel <- sample(1:nrow(original), n_sample)

  dat_plot1 <- unmixing1$unmixed %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    dplyr::slice(sel)
  dat_plot2 <- unmixing2$unmixed %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    dplyr::slice(sel)

  M1 <- unmixing1$M
  M2 <- unmixing2$M

  abundances1 <- as.matrix(dplyr::select(unmixing1$unmixed,
                                         dplyr::all_of(colnames(M1))))
  abundances2 <- as.matrix(dplyr::select(unmixing2$unmixed,
                                         dplyr::all_of(colnames(M2))))

  reconstruction1 <- unmixing1$reconstruction
  reconstruction2 <- unmixing2$reconstruction

  measurement <- dplyr::select(original, dplyr::all_of(rownames(M1)))

  markers <- colnames(M1)

  ui <- miniUI::miniPage(
    miniUI::gadgetTitleBar("My Gadget"),
    miniUI::miniContentPanel(
      shiny::fluidRow(
        shiny::div(
          shiny::div(style="display: inline-block; vertical-align: top;",
                     shiny::uiOutput("marker_x")),
          shiny::div(style="display: inline-block; vertical-align: top;",
                     shiny::radioButtons("marker_x_scaling",
                                         "x-axis scaling",
                                         choices = list("asinh", "linear"),
                                         selected = "asinh",
                                         inline = T)),
          shiny::div(style="display: inline-block; vertical-align: top;",
                     shiny::numericInput("marker_x_cofactor",
                                         label = "x-axis scaling cofactor",
                                         value = 1000,
                                         min = 1))
        ),
        shiny::div(
          shiny::div(style="display: inline-block; vertical-align: top;",
                     shiny::uiOutput("marker_y")),
          shiny::div(style="display: inline-block; vertical-align: top;",
                     shiny::radioButtons("marker_y_scaling",
                                         "y-axis scaling",
                                         choices = list("asinh", "linear"),
                                         selected = "asinh",
                                         inline = T)),
          shiny::div(style="display: inline-block; vertical-align: top;",
                     shiny::numericInput("marker_y_cofactor",
                                         label = "y-axis scaling cofactor",
                                         value = 1000,
                                         min = 1))
        ),
        shiny::div(
          shiny::div(style="display: inline-block; vertical-align: top;",
                     shiny::radioButtons("error_scaling",
                                         "error scaling",
                                         choices = list("log", "linear"),
                                         selected = "log",
                                         inline = T))
        )
      ),
      shiny::fluidRow(
        shiny::div(style="display: inline-block; vertical-align: top;",
                   plotly::plotlyOutput("marker_plot_scatter1")
        ),
        shiny::div(style="display: inline-block; vertical-align: top;",
                   plotly::plotlyOutput("marker_plot_scatter2")
        )
      ),
      shiny::div(
        shiny::div(style="display: inline-block; vertical-align: top;",
                   shiny::radioButtons("signal_scaling",
                                       "signal scaling",
                                       choices = list("asinh", "linear"),
                                       selected = "asinh",
                                       inline = T)),
        shiny::div(style="display: inline-block; vertical-align: top;",
                   shiny::numericInput("signal_cofactor",
                                       label = "signal scaling cofactor",
                                       value = 1000,
                                       min = 1))
      ),
      shiny::fluidRow(
        shiny::div(style="display: inline-block; vertical-align: top;",
                   plotly::plotlyOutput("spectra_plot1")
        ),
        shiny::div(style="display: inline-block; vertical-align: top;",
                   plotly::plotlyOutput("spectra_plot2")
        )
      ),
      shiny::fluidRow(
        plotly::plotlyOutput("scatter_plot")
      )
    )
  )

  server <- function(input, output, session) {
    # Define reactive expressions, outputs, etc.

    output$marker_x <- shiny::renderUI({
      shiny::selectInput("marker_x_selection",
                         label = "x-axis",
                         choices = markers)
    })
    output$marker_y <- shiny::renderUI({
      shiny::selectInput("marker_y_selection",
                         label = "y-axis",
                         choices = markers)
    })

    selection <- shiny::reactiveVal(1)
    shiny::observeEvent(plotly::event_data("plotly_click",
                                                  source = "marker_scatter1",
                                                  priority = "event"), {
      selection(
        as.numeric(plotly::event_data("plotly_click",
                                      source = "marker_scatter1",
                                      priority = "event")$customdata)
      )
    })
    shiny::observeEvent(plotly::event_data("plotly_click",
                                                  source = "marker_scatter2",
                                                  priority = "event"), {
      selection(
        as.numeric(plotly::event_data("plotly_click",
                                      source = "marker_scatter2",
                                      priority = "event")$customdata)
      )
    })

    output$marker_plot_scatter1 <- plotly::renderPlotly({
      if (!is.null(input$marker_x_selection)) {
        if (!is.null(input$marker_y_selection)) {
          if (input$marker_x_scaling == "asinh") {
            x_vals <- asinh(dat_plot1[[input$marker_x_selection]] / input$marker_x_cofactor)
          } else {
            x_vals <- dat_plot1[[input$marker_x_selection]]
          }
          if (input$marker_y_scaling == "asinh") {
            y_vals <- asinh(dat_plot1[[input$marker_y_selection]] / input$marker_y_cofactor)
          } else {
            y_vals <- dat_plot1[[input$marker_y_selection]]
          }
          if (input$error_scaling == "log") {
            color_vals <- log(dat_plot1$error)
          } else {
            color_vals <- dat_plot1$error
          }

          size_vals <- rep(1, length(dat_plot2$id))
          size_vals[dat_plot2$id == selection()] <- 5

          dat_plot1 %>%
            plotly::plot_ly(x = x_vals, y = y_vals,
                            text = dat_plot1$id, color = color_vals,
                            customdata = dat_plot1$id,
                            size = size_vals,
                            type = "scatter", mode = "markers",
                            showlegend = T, hoverinfo = "text",
                            source = "marker_scatter1") %>%
            plotly::event_register('plotly_click')
        }
      }
    })

    output$marker_plot_scatter2 <- plotly::renderPlotly({

      if (!is.null(input$marker_x_selection)) {
        if (!is.null(input$marker_y_selection)) {
          if (input$marker_x_scaling == "asinh") {
            x_vals <- asinh(dat_plot2[[input$marker_x_selection]] / input$marker_x_cofactor)
          } else {
            x_vals <- dat_plot2[[input$marker_x_selection]]
          }
          if (input$marker_y_scaling == "asinh") {
            y_vals <- asinh(dat_plot2[[input$marker_y_selection]] / input$marker_y_cofactor)
          } else {
            y_vals <- dat_plot2[[input$marker_y_selection]]
          }
          if (input$error_scaling == "log") {
            color_vals <- log(dat_plot2$error)
          } else {
            color_vals <- dat_plot2$error
          }

          size_vals <- rep(1, length(dat_plot2$id))
          size_vals[dat_plot2$id == selection()] <- 5

          dat_plot2 %>%
            plotly::plot_ly(x = x_vals, y = y_vals,
                            text = dat_plot2$id, color = color_vals,
                            customdata = dat_plot2$id,
                            size = size_vals,
                            type = "scatter", mode = "markers",
                            showlegend = T, hoverinfo = "text",
                            source = "marker_scatter2") %>%
            plotly::event_register('plotly_click')
        }
      }
    })

    output$scatter_plot <- plotly::renderPlotly({

      dat_plot1 %>%
        plotly::plot_ly(x = ~ `FSC-A`, y = ~ `SSC-A`,
                        type = "histogram2dcontour") %>%
        plotly::add_trace(data = dplyr::filter(dat_plot1, .data$id %in% selection()),
                          x = ~ `FSC-A`, y = ~ `SSC-A`,
                          type = "scatter", mode = "markers",
                          text = ~ id, hoverinfo = "text")
    })

    output$spectra_plot1 <- plotly::renderPlotly({

      components <- do.call(rbind, c(
        purrr::map(seq_along(abundances1[selection(),]), function(j) {
          t(M1)[j,] * abundances1[selection(),j]
        }),
        list(measurement[selection(),], reconstruction1[selection(),])
      )
      ) %>%
        dplyr::as_tibble() %>%
        dplyr::mutate(component = c(names(abundances1[selection(),]),
                                    "original",
                                    "reconstruction")) %>%
        tidyr::pivot_longer(cols = rownames(M1),
                            names_to = "channel",
                            values_to = "signal") %>%
        dplyr::mutate(channel = factor(.data$channel, rownames(M1)))

      if (input$signal_scaling == "asinh") {
        components <- components %>%
          dplyr::mutate(signal = asinh(.data$signal/input$signal_cofactor))
      }

      components %>%
        plotly::plot_ly(x = ~ channel, y = ~ signal,
                        color = ~ component,
                        text = ~ component, hoverinfo = "text",
                        type = "scatter", mode = "lines")
    })

    output$spectra_plot2 <- plotly::renderPlotly({

      components <- do.call(rbind, c(
        purrr::map(seq_along(abundances2[selection(),]), function(j) {
          t(M2)[j,] * abundances2[selection(),j]
        }),
        list(measurement[selection(),], reconstruction2[selection(),])
      )
      ) %>%
        dplyr::as_tibble() %>%
        dplyr::mutate(component = c(names(abundances2[selection(),]),
                                    "original",
                                    "reconstruction")) %>%
        tidyr::pivot_longer(cols = rownames(M2),
                            names_to = "channel",
                            values_to = "signal") %>%
        dplyr::mutate(channel = factor(.data$channel, rownames(M2)))

      if (input$signal_scaling == "asinh") {
        components <- components %>%
          dplyr::mutate(signal = asinh(.data$signal/input$signal_cofactor))
      }

      components %>%
        plotly::plot_ly(x = ~ channel, y = ~ signal,
                        color = ~ component,
                        text = ~ component, hoverinfo = "text",
                        type = "scatter", mode = "lines")
    })

    # When the Done button is clicked, return a value
    shiny::observeEvent(input$done, {
      returnValue <- NULL
      shiny::stopApp(returnValue)
    })
  }

  shiny::runGadget(ui, server)
}


#' Plot one marker against all others
#'
#' Controls the unmixing quality by plotting one marker against all others.
#'
#' Best done on data from single stained cells that was unmixed with a given unmixing matrix. If the plots show any curved distributions, there is a problem with the unmixing (e.g., a single stain needs to be recorded on cells instead of beads).
#'
#' @param marker string, the marker to be plotted agains all other markers
#' @param dat_unmixed dataframe, the unmixed data to be plotted. Can contain data from multiple single stains.
#' @param markers vector, all markers to be included in the plot
#' @param output either "screen", "file" or "both". Whether to return the plot ("screen"), write it to a png ("file") or do both
#' @param output_path string, where to save the plot if output is set to either "file" or "both". Default: ".", the script directory.
#'
#' @return either a ggplot2 plot, the plot as a png file or both
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
PlotMxN <- function(marker, dat_unmixed,
                    markers = NULL,
                    output = "screen", output_path = ".") {

  # take dat_unmixed directly from Unmix output
  if ("unmixed" %in% names(dat_unmixed)) {
    dat_unmixed <- dat_unmixed$unmixed
  }

  if (!marker %in% unique(dat_unmixed$file)) {
    stop("No data found for ", marker, " in the input data.")
  }

  if (is.null(markers)) {
    markers <- colnames(dat_unmixed)
    markers <- markers[!stringr::str_detect(markers, "error|Time|SSC|FSC|file")]
  }

  if (!all(markers %in% colnames(dat_unmixed))) {
    stop("The following markers were not found in the data: ",
         stringr::str_c(markers[!markers %in% colnames(dat_unmixed)],
                        collapse = " "))
  } else if (!marker %in% markers) {
    stop(marker, " not found in data columns.")
  }

  plot <- dat_unmixed %>%
    dplyr::filter(file == marker) %>%
    dplyr::select(dplyr::all_of(markers)) %>%
    tidyr::pivot_longer(cols = markers[markers != marker],
                 names_to = "channel", values_to = "signal") %>%
    dplyr::mutate(channel = factor(.data$channel, markers)) %>%
    ggplot2::ggplot(ggplot2::aes(.data[[!!rlang::sym(marker)]], .data$signal)) +
    ggplot2::geom_hex(bins = 300, show.legend = F) +
    ggplot2::scale_fill_gradientn(colours = rev(grDevices::rainbow(20,
                                                                   start = 0,
                                                                   end = 0.7))) +
    ggplot2::facet_wrap(ggplot2::vars(.data$channel), ncol = 5) +
    ggcyto::scale_x_logicle(t=2^22) +
    ggcyto::scale_y_logicle(t=2^22) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(marker)

  if (output %in% c("file", "both")) {
    dir.create(output_path, showWarnings = FALSE)

    ggplot2::ggsave(file.path(output_path, stringr::str_c(marker, ".png")),
                    plot,
                    width = 7.5,
                    height = ceiling((length(markers)-1)/5) * 1.5 + 1)
  }

  if (output %in% c("screen", "both")) {
    plot
  }
}


