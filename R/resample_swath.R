# Internal function called by swath_to_geofield so shouldn't
# need to check the inputs

resample_swath <- function(lon, lat, value, opts) {

  # reproject pixel lon and lat to projection of dom
  dom <- meteogrid::as.geodomain(opts[["dom"]])
  xy  <- meteogrid::project(lon, lat, dom[["projection"]])

  xy[["value"]] <- value

  # Remove swath pixels that are outside of the domain
  dom_extent <- meteogrid::DomainExtent(dom)
  xy <- dplyr::filter(
    xy,
    .data[["x"]] >= dom_extent[["x0"]],
    .data[["x"]] <= dom_extent[["x1"]],
    .data[["y"]] >= dom_extent[["y0"]],
    .data[["y"]] <= dom_extent[["y1"]]
  )

  if (nrow(xy) < 1) {
    stop("No pixels found inside `dom`.", call. = FALSE)
  }

  # Make a subdomain that includes the swath to reduce search time.
  # The subgrid should be the extent of the swath plus the radius
  # unless at the edge of the domain
  dom_x <- seq(dom_extent[["x0"]], dom_extent[["x1"]], dom_extent[["dx"]])
  dom_y <- seq(dom_extent[["y0"]], dom_extent[["x1"]], dom_extent[["dy"]])

  x0 <- max(c(
    1,
    which.min(abs(dom_x - (min(xy[["x"]]) - opts[["radius"]] * 2))) - 1
  ))
  x1 <- min(c(
    dom_extent[["nx"]],
    which.min(abs(dom_x - (max(xy[["x"]]) + opts[["radius"]] * 2))) + 1
  ))
  y0 <- max(c(
    1,
    which.min(abs(dom_y - (min(xy[["y"]]) - opts[["radius"]] * 2))) - 1
  ))
  y1 <- min(c(
    dom_extent[["ny"]],
    which.min(abs(dom_y - (max(xy[["y"]]) + opts[["radius"]] * 2))) + 1
  ))

  sub_dom        <- meteogrid::subgrid(dom, x0, x1, y0, y1)
  sub_dom_extent <- meteogrid::DomainExtent(sub_dom)

  # Use kd tree to find nearest neighbours within radius
  grid_points <- expand.grid(
    x = seq(sub_dom_extent[["x0"]], sub_dom_extent[["x1"]], sub_dom_extent[["dx"]]),
    y = seq(sub_dom_extent[["y0"]], sub_dom_extent[["y1"]], sub_dom_extent[["dy"]])
  )

  resample_func <- get(paste0("resample_", opts[["method"]]))
  resample_func(xy, grid_points, opts[["radius"]])

}

resample_nearest <- function(xy, grid_points, radius) {

  nn <- RANN::nn2(
    xy[c("x", "y")],
    grid_points,
    k          = 1,
    searchtype = "radius",
    radius     = radius
  )

  nn[["nn.dists"]][nn[["nn.dists"]] > 1e100] <- NA

  idx         <- which(nn[["nn.idx"]] != 0)
  values      <- rep(NA, length(nn[["nn.idx"]]))
  values[idx] <- xy[["value"]][nn[["nn.idx"]][idx]]

  cbind(cbind(grid_points, value = values), distance = nn[["nn.dists"]])

}

resample_inverse_distance <- function(xy, grid_points, radius) {

  nn <- RANN::nn2(
    xy[c("x", "y")],
    grid_points,
    searchtype = "radius",
    radius     = radius
  )

  nn[["nn.dists"]][nn[["nn.dists"]] > 1e100] <- NA

  idx         <- which(nn[["nn.idx"]] != 0)
  values      <- array(0, dim(nn[["nn.idx"]]))
  values[idx] <- xy[["value"]][nn[["nn.idx"]][idx]]

  weights <- (1 / nn[["nn.dists"]]) /
    rowSums(1 / nn[["nn.dists"]], na.rm = TRUE)

  counts <- rowSums(values / values, na.rm = TRUE)
  values <- rowSums(values * weights, na.rm = TRUE)
  values[counts == 0] <- NA

  cbind(cbind(grid_points, value = values), count = counts)

}

resample_mean <- function(xy, grid_points, radius) {

  nn <- RANN::nn2(
    xy[c("x", "y")],
    grid_points,
    searchtype = "radius",
    radius     = radius
  )

  nn[["nn.dists"]][nn[["nn.dists"]] > 1e100] <- NA

  idx         <- which(nn[["nn.idx"]] != 0)
  values      <- array(NA, dim(nn[["nn.idx"]]))
  values[idx] <- xy[["value"]][nn[["nn.idx"]][idx]]

  counts <- rowSums(values / values, na.rm = TRUE)
  values <- rowMeans(values, na.rm = TRUE)
  values[counts == 0] <- NA

  cbind(cbind(grid_points, value = values), count = counts)

}
