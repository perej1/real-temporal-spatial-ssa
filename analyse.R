library(tidyverse)
library(sftime)
library(optparse)
#library(SpaceTimeBSS)

option_list <- list(
  make_option("--xblocks", type = "character", default = "33:33:34",
              help = stringr::str_c("Segmentation of x coord, string gives ",
                                    "proportions of the segment lengths")),
  make_option("--yblocks", type = "character", default = "33:33:34",
              help = stringr::str_c("Segmentation of y coord, string gives ",
                                    "proportions of the segment lengths")),
  make_option("--tblocks", type = "character", default = "33:33:34",
              help = stringr::str_c("Segmentation of time, string gives ",
                                    "proportions of the segment lengths")),
              make_option("--dim_nonstationary", type = "integer", default = 1,
                          help = "Number of nonstationary components")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

source("functions.R")
load("data/veneto.RData")

# Create sftime object with the desired variables
veneto <- file_veneto %>%
  select(X, Y, codt, ET0, tmax, hmax, wind) %>%
  rename(x = X, y = Y, time = codt, evapotranspiration = ET0, temp_max = tmax,
         humid_max = hmax) %>%
  as_tibble() %>%
  sf::st_as_sf(coords = c("x", "y"), sf_column_name = "geometry",
               crs = 3003) %>%
  sftime::st_sftime(time_column_name = "time")

# Data without locations and time
observed <- veneto %>%
  sftime::st_drop_time() %>%
  sf::st_drop_geometry()

# Coordinates
coords <- veneto %>%
  select(geometry, time)

# Compute whitened field
dim <- ncol(observed)
mean_p <- colMeans(observed)
cov_p_inv_sqrt <- sqrtmat_inv(cov(observed))

whitened <- observed %>%
  sweep(2, mean_p, "-") %>%
  apply(1, function(x) cov_p_inv_sqrt %*% matrix(x, ncol = 1)) %>%
  t() %>%
  as_tibble(.name_repair = "minimal")
colnames(whitened) <- stringr::str_c("f", 1:dim)

# Parse cut points
x_prop <- stringr::str_split(opt$xblocks, ":", simplify = TRUE) %>%
  as.integer()
y_prop <- stringr::str_split(opt$yblocks, ":", simplify = TRUE) %>%
  as.integer()
time_prop <- stringr::str_split(opt$tblocks, ":", simplify = TRUE) %>%
  as.integer()

whitened_with_seg <- compute_segments(veneto, x_prop, y_prop, time_prop) %>%
  cbind(whitened) %>%
  group_by(x_segment, y_segment, time_segment)

# Compute means and proportional sample sizes in segments
means_segment <- whitened_with_seg %>%
  summarise(across(starts_with("f"), ~ mean(.x)), .groups = "drop") %>%
  mutate(n_prop = group_size(whitened_with_seg) / (nrow(veneto)))

# Compute variance of the segment means
outer_prods <- means_segment %>%
  select(starts_with("f")) %>%
  as.matrix() %>%
  apply(1, function(x) x %o% x, simplify = FALSE)

mean_var <- purrr::pmap(list(arg1 = means_segment$n_prop, arg2 = outer_prods),
                        \(arg1, arg2) arg1 * arg2) %>%
  purrr::reduce(`+`)


# Plot scree plot
eigen_val <- eigen(mean_var)$values
plot(1:4, eigen_val)

scree_name <- stringr::str_c("xblocks_", opt$xblocks,
                             "_yblocks_", opt$yblocks,
                             "_tblocks_", opt$tblocks, ".pdf")

scree <- tibble(x = 1:4, y = eigen_val) %>%
  ggplot(aes(x, y)) +
  geom_point(size = 3) +
  geom_line() +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        panel.grid.major.y = element_line(colour = alpha("black", 0.2)),
        panel.grid.major.x = element_line(colour = alpha("black", 0.2))) +
  xlab(expression(lambda[i])) +
  ylab("Value") +
  ylim(0, 0.35)
  ggsave(stringr::str_c("plots/scree/", scree_name), scree)

# Compute unmixing matrix
v_transpose <- t(eigen(mean_var)$vectors)
w <- v_transpose %*% cov_p_inv_sqrt
w_nonstationary <- w[1:opt$dim_nonstationary, , drop = FALSE]
w_stationary <- w[(opt$dim_nonstationary + 1):dim, , drop = FALSE]
