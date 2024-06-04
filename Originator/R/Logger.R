#' Logger
#' @description
#' A logger object
#'
#' @section Functions:
#' - \code{\link[=Logger.log]{Logger$log}}
#' - \code{\link[=Logger.info]{Logger$info}}
#' - \code{\link[=Logger.debug]{Logger$debug}}
#' - \code{\link[=Logger.warn]{Logger$warn}}
#' - \code{\link[=Logger.error]{Logger$error}}
#'
#' @export
Logger <- new.env()

#' Logger$log
#' @name Logger.log
#' @description
#' Log a UTC timestamped message with specified level
#'
#' @param text Text message
#' @param level Log level, default INFO
Logger$log <- function(text = "", level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC")
  message(sprintf("[%s] [%s] %s", timestamp, level, text))
}

#' Logger$info
#' @name Logger.info
#' @description
#' Log a INFO level message
#'
#' @param text Text message
Logger$info <- function(text = "") {
  Logger$log(text)
}

#' Logger$debug
#' @name Logger.debug
#' @description
#' Log a DEBUG level message
#'
#' @param text Text message
Logger$debug <- function(text = "") {
  Logger$log(text, "DEBUG")
}

#' Logger$warn
#' @name Logger.warn
#' @description
#' Log a WARN level message
#'
#' @param text Text message
Logger$warn <- function(text = "") {
  Logger$log(text, "WARN")
}

#' Logger$error
#' @name Logger.error
#' @description
#' Log a ERROR level message
#'
#' @param text Text message
Logger$error <- function(text = "") {
  Logger$log(text, "ERROR")
}
