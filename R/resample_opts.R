resample_opts <- function(
  dom,
  radius,
  method = c("inverse_distance", "nearest", "mean")
) {

  if (missing(dom)) {
    stop("`dom` must be passed to resample_opts.", call. = FALSE)
  }

  if (missing(radius)) {
    stop("`radius` must be passed to resample_opts.", call. = FALSE)
  }

  if (!(meteogrid::is.geofield(dom) || meteogrid::is.geodomain(dom))) {
    stop("`dom` must be a geofield or geodomain.", call. = FALSE)
  }

  if (!(is.numeric(radius) && length(radius) == 1)) {
    stop("`radius` must be a numeric scalar.", call. = FALSE)
  }

  method <- match.arg(method)

  list(dom = dom, radius = radius, method = method)
}
