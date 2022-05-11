swath_to_geofield <- function(
  swath_df,
  data_col,
  dom,
  search_radius,
  group_col = NULL,
  method    = c("inverse_distance", "nearest", "mean"),
  output    = c("data_frame", "geofield"),
  ...
) {

  method <- match.arg(method)
  output <- match.arg(output)

  data_col_quo  <- rlang::enquo(data_col)
  group_col_quo <- rlang::enquo(group_col)

  data_col_name  <- rlang::as_label(data_col_quo)
  expected_cols <- c("lon", "lat", data_col_name)

  if (
    !rlang::quo_is_null(group_col_quo) &&
      !rlang::quo_is_missing(group_col_quo)
  ) {
    group_col_name <- rlang::as_label(group_col_quo)
    expected_cols <- c(expected_cols, group_col_name)
  }

  missing_cols  <- setdiff(expected_cols, colnames(swath_df))

  if (length(missing_cols) > 0) {
    stop("`swath_df` must contain columns: `", paste(expected_cols, collapse = "`, `"), "`")
  }

  dom <- try(meteogrid::as.geodomain(dom))
  if (inherits(dom, "try-error")) {
    stop("`dom` must be a `geofield` or `geodomain`.")
  }

  indices <- meteogrid::point.index(
    swath_df$lon,
    swath_df$lat,
    domain = dom
  )

  if (length(which(!is.na(indices[["i"]]) && !is.na(indices[["j"]]))) < 1) {
    stop("Whole swath is outside of domain.", immediate. = TRUE)
  }

  swath_raster <- harpIO::geofield_to_raster(
    meteogrid::as.geofield(NA_real_, dom)
  )

  proj4str <- meteogrid::proj4.list2str(dom$projection)

  coords   <- swath_df[c("lon", "lat")]
  swath_df <- swath_df[!colnames(swath_df) %in% c("lon", "lat")]

  swath_df <- sp::SpatialPointsDataFrame(
    coords      = data.frame(coords),
    data        = data.frame(swath_df),
    proj4string = sp::CRS("+proj=longlat")
  ) %>%
    sp::spTransform(proj4str)

  swath_raster <- raster::rasterize(
    swath_df, swath_raster, field = data_col_name#, fun = mean
  )

  harpIO::raster_to_geofield(swath_raster)

}
