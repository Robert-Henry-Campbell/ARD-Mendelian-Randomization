# Internal helpers for attaching and exporting plot data frames

#' @keywords internal
.ardmr_sanitize_component <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- "component"
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- sub("^_+", "", x)
  x <- sub("_+$", "", x)
  x <- ifelse(!nzchar(x), "component", x)
  x
}

#' @keywords internal
.ardmr_coerce_data_list <- function(x, default_prefix = "data") {
  out <- list()
  counter <- 1L

  add_item <- function(item, nm = NULL) {
    if (is.null(item)) return()

    if (is.data.frame(item)) {
      name_use <- nm
      if (is.null(name_use) || !nzchar(name_use)) {
        name_use <- sprintf("%s%d", default_prefix, counter)
        counter <<- counter + 1L
      }
      out[[name_use]] <<- tibble::as_tibble(item)
      return()
    }

    if (is.list(item)) {
      nms <- names(item)
      if (length(item) == 0) return()
      for (i in seq_along(item)) {
        nm_i <- nms[i]
        if (!is.null(nm) && nzchar(nm)) {
          if (is.null(nm_i) || !nzchar(nm_i)) {
            nm_i <- nm
          } else {
            nm_i <- paste(nm, nm_i, sep = "__")
          }
        }
        add_item(item[[i]], nm_i)
      }
      return()
    }
  }

  if (is.list(x) && !is.data.frame(x)) {
    nms <- names(x)
    if (length(x) == 0) return(out)
    for (i in seq_along(x)) {
      nm_i <- nms[i]
      add_item(x[[i]], nm_i)
    }
  } else {
    add_item(x, NULL)
  }

  out
}

#' @keywords internal
.ardmr_attach_plot_data <- function(plot, ...) {
  data_list <- list(...)
  if (length(data_list) == 1L && is.list(data_list[[1]]) &&
      (is.null(names(data_list)) || !nzchar(names(data_list)[1]))) {
    data_list <- data_list[[1]]
  }
  data_list <- .ardmr_coerce_data_list(data_list)
  attr(plot, "ardmr_plot_data") <- data_list
  plot
}

#' @keywords internal
.ardmr_collect_plot_data <- function(plot, fallback = NULL) {
  data_list <- attr(plot, "ardmr_plot_data", exact = TRUE)
  data_list <- .ardmr_coerce_data_list(data_list)
  if (length(data_list)) return(data_list)

  fallback_list <- .ardmr_coerce_data_list(fallback)
  if (length(fallback_list)) return(fallback_list)

  if (is.data.frame(plot$data)) {
    return(.ardmr_coerce_data_list(plot$data))
  }

  built <- tryCatch(ggplot2::ggplot_build(plot), error = function(e) NULL)
  if (!is.null(built) && length(built$data)) {
    built_list <- .ardmr_coerce_data_list(built$data)
    if (length(built_list)) return(built_list)
  }

  list()
}

#' @keywords internal
.ardmr_write_plot_data <- function(plot, dir_path, base_name, fallback = NULL) {
  data_list <- .ardmr_collect_plot_data(plot, fallback)
  if (!length(data_list)) return(invisible(NULL))

  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }

  nms <- names(data_list)
  if (is.null(nms)) nms <- rep("", length(data_list))
  if (!length(data_list)) return(invisible(NULL))

  # ensure unique, sanitized names
  if (length(data_list) == 1L) {
    dataset_names <- "data"
  } else {
    dataset_names <- ifelse(nzchar(nms), nms, paste0("dataset", seq_along(data_list)))
    dataset_names <- make.unique(dataset_names, sep = "_")
  }

  for (i in seq_along(data_list)) {
    df <- data_list[[i]]
    if (!is.data.frame(df)) next
    nm <- if (length(data_list) == 1L) {
      "data"
    } else {
      dataset_names[i]
    }
    file_stub <- if (length(data_list) == 1L) {
      paste0(base_name, "__data")
    } else {
      paste0(base_name, "__", .ardmr_sanitize_component(nm))
    }
    file_path <- file.path(dir_path, paste0(file_stub, ".csv"))
    if (requireNamespace("readr", quietly = TRUE)) {
      readr::write_csv(df, file_path)
    } else {
      utils::write.csv(df, file_path, row.names = FALSE)
    }
  }

  invisible(NULL)
}
