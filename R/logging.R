#' Set up logging for ardmr
#'
#' Configures the [logger] package to emit messages to the console and,
#' optionally, to a file. The log level controls verbosity.
#'
#' @param logfile Path to a log file. If empty, logs only to the console.
#' @param level Logging level threshold passed to [logger::log_threshold()].
#' @return Invisibly returns the logfile path.
#' @export
setup_logging <- function(logfile, level = "INFO") {
  if (nzchar(logfile)) {
    logger::log_appender(logger::appender_tee(logfile))
  } else {
    logger::log_appender(logger::appender_console())
  }
  logger::log_layout(logger::layout_glue_colors)
  logger::log_threshold(level)
  logger::log_info("Logging initialised at level {level}")
  invisible(logfile)
}
