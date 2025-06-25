library(tidyverse)
library(sftime)
library(optparse)
library(tmap)

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
  make_option("--seed", type = "integer", default = 500,
              help = stringr::str_c("Seed for plotting the field on a",
                                    "fixed time/spatial point"))
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

filename <- stringr::str_c("xblocks_", opt$xblocks,
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
ggsave(stringr::str_c("plots/scree/", filename), scree)

# Compute unmixing matrix and latent components
v_transpose <- t(eigen(mean_var)$vectors)
w <- v_transpose %*% cov_p_inv_sqrt

latent <- t(apply(observed, 1, function(x) w %*% x)) %>%
  as_tibble(.name_repair = "universal")
colnames(latent) <- stringr::str_c("f", 1:dim)

# Add coordinates to the latent components
latent <- latent %>%
  mutate(time = veneto$time, geometry = veneto$geometry) %>%
  sftime::st_sftime(time_column_name = "time", sf_column_name = "geometry")

# Plot first and last component at a random location for all time points
set.seed(opt$seed) # 500
location <- latent$geometry[sample(1:nrow(latent), 1)]
latent_location <- latent %>%
  filter(geometry == location) %>%
  select(f1, f4, time) %>%
  sf::st_drop_geometry()


fix_location <- latent_location %>%
  gather("field", "value", -time) %>%
  ggplot(aes(x = time, y = value, linetype = factor(field))) +
  geom_line() +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        panel.grid.major.y = element_line(colour = alpha("black", 0.2)),
        panel.grid.major.x = element_line(colour = alpha("black", 0.2))) +
  xlab("Time") +
  ylab("Value") +
  scale_linetype_manual(values = c("f1" = "solid", "f4" = "dashed"),
                        name = "Fields",
                        labels = c("f1" = "Nonstationary", "f4" = "Stationary"))
ggsave(stringr::str_c("plots/paths/", "fix_location_", filename),
       fix_location)

# Plot first and last component at a random time point for all time locations
t <- latent$time[sample(1:nrow(latent), 1)]
latent_time <- latent %>%
  filter(time == t) %>%
  select(f1, f4, time) %>%
  sftime::st_drop_time() %>%
  mutate(f1 = f1  - mean(f1), f4 = f4  - mean(f4))

latent_time_f1 <- latent_time %>%
  select(-f4)

latent_time_f4 <- latent_time %>%
  select(-f1)


veneto_box <- sf::st_buffer(veneto$geometry, dist = 50 * 1000) %>%
  st_bbox()

sf_italy <-
  rnaturalearth::ne_countries(returnclass = "sf", scale = "medium") %>%
  dplyr::filter(name == "Italy") %>%
  st_transform(crs = 3003)

# Show Veneto area in Italy
pdf("plots/veneto.pdf")
tm_shape(sf_italy) +
  tm_polygons() +
  tm_shape(st_as_sfc(veneto_box)) +
  tm_polygons(col = "red", fill_alpha = 0)
dev.off()


n_break <- 10
minbreak <- min(c(range(latent_time_f1$f1), range(latent_time_f4$f4)))
maxbreak <- max(c(range(latent_time_f1$f1), range(latent_time_f4$f4)))
breaks <- seq(minbreak, maxbreak, length.out = n_break + 1)


# Show nonstationary field at a fixed time point
pdf(stringr::str_c("plots/paths/", "fix_time_nonstationary_", filename))
tm_shape(sf_italy, bbox = veneto_box) +
  tm_polygons() +
  tm_shape(latent_time_f1) +
  tm_bubbles(size = 1, fill = "f1",
             fill.scale = tm_scale_intervals(n = n_break, style = "fixed",
                                             breaks = breaks,
                                             midpoint = NA,
                                             values = "-brewer.spectral"),
             fill.legend = tm_legend("Nonstationary"))
dev.off()

# Show stationary field at a fixed time point
pdf(stringr::str_c("plots/paths/", "fix_time_stationary_", filename))
tm_shape(sf_italy, bbox = veneto_box) +
  tm_polygons() +
  tm_shape(latent_time_f4) +
  tm_bubbles(size = 1, fill = "f4",
             fill.scale = tm_scale_intervals(n = n_break, style = "fixed",
                                             breaks = breaks,
                                             midpoint = NA,
                                             values = "-brewer.spectral"),
             fill.legend = tm_legend("Stationary"))
dev.off()
