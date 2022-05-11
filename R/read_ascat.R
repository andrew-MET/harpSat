#' Title
#'
#' @param file_name
#' @param parameter
#' @param rounding
#' @param bounding_date_times
#' @param bounding_domain
#' @param resample_domain
#'
#' @return
#' @export
#'
#' @examples
read_ascat_nc <- function(
  file_name,
  parameter,
  rounding            = c(NA, "secs", "mins", "hours"),
  bounding_date_times = NA,
  bounding_domain     = NA,
  resample_options    = NA,
  lonlat_to_id        = FALSE
) {

  rounding <- match.arg(rounding)
  if (is.list(resample_options)) {
    resample_options <- do.call(resample_opts, resample_options)
  }

  has_bbox <- TRUE

  if (!is.na(bounding_date_times)) {
    if (length(bounding_date_times) != 2) {
      stop("`bounding_date_times` must be a vector of length 2")
    }
    bounding_date_times <- sort(harpIO::str_datetime_to_datetime(bounding_date_times))
  }

  bounding_domain_try <- try(meteogrid::as.geodomain(bounding_domain), silent = TRUE)
  if (inherits(bounding_domain_try, "try-error")) {
    if (!is.na(bounding_domain)) {
      stop("`bounding_domain` must be a `geofield` or `geodomain`")
    }
    has_bbox <- FALSE
    bounding_domain_try <- NA
  }
  bounding_domain <- bounding_domain_try

  nc_id <- ncdf4::nc_open(file_name)

  nc_vars <- names(nc_id[["var"]])
  params  <- intersect(parameter, nc_vars)

  if (length(params) < 1) {
    warning("None of the requested paramters found in file.", immediate. = TRUE)
    ncdf4::nc_close(nc_id)
    return(NULL)
  }

  missing_params <- setdiff(sort(parameter), sort(params))

  if (length(missing_params) > 0) {
    warning("`", paste(missing_params, sep = "`, `"), "` not found in file.", immediate. = TRUE)
  }

  times      <- ncdf4::ncvar_get(nc_id, "time")
  time_units <- ncdf4::ncatt_get(nc_id, "time", "units")

  origin <- sub("^ ", "", gsub("[[:alpha:]]+ since", "", time_units[["value"]]))
  times  <- times * harpIO:::units_multiplier(substr(time_units[["value"]], 1, 1))
  times  <- as.POSIXct(times, tz = "UTC", origin = origin)

  if (!is.na(rounding)) {
    times <- round(times, rounding)
  }

  time_range <- range(times)

  if (!is.na(bounding_date_times)) {
    if (max(time_range) < min(bounding_date_times) | min(time_range) > max(bounding_date_times)) {
      warning(
        "Times in file are outside of `bounding_date_times`.\n",
        "`bounding_date_times`: ", paste(bounding_date_times, collapse = " - "), "\n",
        "Time range in file   : ", paste(time_range, collapse = " - "),
        immediate. = TRUE
      )
      ncdf4::nc_close(nc_id)
      return(NULL)
    }
  }

  lat <- ncdf4::ncvar_get(nc_id, "lat")
  lon <- ncdf4::ncvar_get(nc_id, "lon")

  if (has_bbox) {

    indices <- meteogrid::point.index(
      lon    = as.vector(lon),
      lat    = as.vector(lat),
      domain = bounding_domain
    )

    if (!any(complete.cases(indices))) {
      warning("Whole swath is outside of `bounding_domain`.", immediate. = TRUE)
      ncdf4::nc_close(nc_id)
      return(NULL)
    }

  }

  swath_df <- tibble::tibble(
    validdate = times[1:length(times)],
    lon       = lon[1:length(lon)],
    lat       = lat[1:length(lat)],
    SID       = seq_along(times)
  )

  if (lonlat_to_id) {
    swath_df[["SID"]] <- paste(
      round(swath_df[["lon"]], 5), round(swath_df[["lat"]], 5), sep = "_"
    )
  }

  for (param in params) {
    swath_df[[param]] <- as.vector(ncdf4::ncvar_get(nc_id, param))
  }

  ncdf4::nc_close(nc_id)

  if (is.list(resample_options)) {
    return(
      resample_swath(
        swath_df$lon, swath_df$lat, swath_df[[param]], resample_options
      )
    )
  }

  if (has_bbox) {
    swath_df <- cbind(swath_df, indices)
    swath_df <- dplyr::filter(
      swath_df,
      dplyr::if_all(dplyr::all_of(c("i", "j")), ~!is.na(.x))
    )
    swath_df <- cbind(
      swath_df,
      meteogrid::project(swath_df$lon, swath_df$lat, bounding_domain[["projection"]])
    )
    attr(swath_df, "domain") <- bounding_domain
  }

  tibble::as_tibble(
    dplyr::filter(
      swath_df, dplyr::if_any(dplyr::all_of(params), ~!is.na(.x))
    )
  )
}


