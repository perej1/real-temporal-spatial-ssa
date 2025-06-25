#' Compute negative square root of a positive definite matrix
#'
#' More precisely, function computes matrix lambda
#' s.t. sigma^{-1} = lambda %*% lambda.
#'
#' @param sigma A Double matrix, positive definite matrix.
#'
#' @return Double matrix, square root of the matrix.
sqrtmat_inv <- function(sigma) {
  eigenval <- eigen(sigma)$values
  if (any(eigenval <= 0) || any(sigma != t(sigma))) {
    rlang::abort("`sigma` must be a symmetric positive definite matrix.")
  }
  eigenvec <- eigen(sigma)$vectors
  eigenvec %*% diag(eigenval^(-0.5)) %*% t(eigenvec)
}


#' Compute segments for spatio-temporal coordinates
#'
#' @param data An sftime object with locations and time
#' @param x_prop,y_prop,time_prop Division of the corresponding axis in terms of
#'   percentages. Percentages are given as a double vector.
#'
#' @returns Tibble of segments for different axes (3 columns).
compute_segments <- function(data, x_prop, y_prop, time_prop) {
  if (sum(x_prop) != 100 || sum(y_prop) != 100 || sum(time_prop) != 100) {
    rlang::abort("Sum of proportions must be 100.")
  }
  
  time <- sftime::st_time(data) %>%
    unique()
  
  geom <- sf::st_coordinates(data) %>%
    unique() %>%
    as_tibble() %>%
    rename_with(tolower)
  
  xlen <- max(geom$x) - min(geom$x)
  x_cuts <- min(geom$x) + c(0, cumsum(x_prop / 100)) * xlen
  
  ylen <- max(geom$y) - min(geom$y)
  y_cuts <- min(geom$y) + c(0, cumsum(y_prop / 100)) * ylen
  
  time_cuts <- quantile(1:max(time), c(0, cumsum(time_prop / 100)))
  
  data %>%
    mutate(x_segment = cut(sf::st_coordinates(.)[, "X"], x_cuts,
                           include.lowest = TRUE),
           y_segment = cut(sf::st_coordinates(.)[, "Y"], y_cuts,
                           include.lowest = TRUE),
           time_segment = cut(time, time_cuts, include.lowest = TRUE)) %>%
    sftime::st_drop_time() %>%
    sf::st_drop_geometry() %>%
    select(x_segment, y_segment, time_segment)
}
