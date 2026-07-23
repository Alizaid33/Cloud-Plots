
# R code to Construct the Cloud plot for spherical data
# with spherical coordinates (longitude,colatitude) = (\phi ,\theta)

# Required libraries
library(sphunif)
library(rgl)
library(geometry)
library(Directional)
library(readxl)

set.seed(123)

# Helper Functions
geodesic_distance <- function(p1, p2) {
  dp <- pmin(pmax(sum(p1 * p2), -1), 1)
  acos(dp)
}
cdf_geodesic_distance <- function(d, kappa) {
  if (kappa == 0) return(d / pi)
  (exp(kappa) - exp(kappa * cos(d))) / (2 * sinh(kappa))
}
find_percentile <- function(target_probability, kappa) {
  objective_function <- function(d) cdf_geodesic_distance(d, kappa) - target_probability
  if (target_probability == 0) return(0)
  if (target_probability == 1) return(pi)
  if (kappa == 0) return(target_probability * pi)
  tryCatch(
    uniroot(objective_function, c(.Machine$double.eps, pi - .Machine$double.eps))$root,
    error = function(e) NA
  )
}
estimate_kappa_3D <- function(data) {
  R_bar <- sqrt(sum(colMeans(data)^2))
  if (is.na(R_bar) || R_bar < 0 || R_bar >= 1 - .Machine$double.eps) {
    if (R_bar >= 1 - .Machine$double.eps) return(Inf)
    if (R_bar <= .Machine$double.eps) return(0)
    return(NA)
  }
  (3 * R_bar - R_bar^3) / (1 - R_bar^2)
}
cross_product <- function(u, v) {
  c(u[2]*v[3] - u[3]*v[2], u[3]*v[1] - u[1]*v[3], u[1]*v[2] - u[2]*v[1])
}
generate_spherical_region <- function(center, inner_radius, outer_radius, nr=20, nphi=40) {
  if (sqrt(center[1]^2 + center[2]^2) < 1e-6) {
    e1 <- c(1, 0, 0)
  } else {
    e1 <- c(-center[2], center[1], 0)
    e1 <- e1 / sqrt(sum(e1^2))
  }
  e2 <- cross_product(center, e1); e2 <- e2 / sqrt(sum(e2^2))
  rs <- seq(inner_radius, outer_radius, length.out = nr + 1)
  phis <- seq(0, 2*pi, length.out = nphi + 1)
  vertices <- matrix(NA, nrow = (nr+1)*(nphi+1), ncol = 3)
  idx <- 1
  for (r in rs) {
    for (phi in phis) {
      p <- cos(r)*center + sin(r)*(cos(phi)*e1 + sin(phi)*e2)
      vertices[idx, ] <- p / sqrt(sum(p^2))
      idx <- idx + 1
    }
  }
  faces <- matrix(NA, nrow = 2 * nr * nphi, ncol = 3)
  idx <- 1
  for (i in 1:nr) {
    for (j in 1:nphi) {
      a <- (i-1)*(nphi+1) + j
      b <- a + 1
      c <- a + (nphi+1)
      d <- c + 1
      faces[idx, ] <- c(a, b, c)
      faces[idx+1, ] <- c(b, d, c)
      idx <- idx + 2
    }
  }
  list(vertices = vertices, faces = faces)
}
plot_spherical_region <- function(center, inner_radius, outer_radius, col_shade, alpha_val=0.3, nr=20, nphi=40) {
  reg <- generate_spherical_region(center, inner_radius, outer_radius, nr, nphi)
  mesh <- tmesh3d(vertices = t(reg$vertices), indices = t(reg$faces), homogeneous = FALSE)
  shade3d(mesh, color = col_shade, alpha = alpha_val)
}

add_spherical_grid <- function(nlat = 5, nlon = 8) {
  # Latitude lines (constant colatitude theta)
  thetas <- seq(0, pi, length.out = nlat + 1)
  for (th in thetas[-c(1, length(thetas))]) {
    phi <- seq(0, 2*pi, length.out = 100)
    x <- sin(th) * cos(phi)
    y <- sin(th) * sin(phi)
    z <- cos(th) #* rep(1, 100)
    lines3d(x, y, z, col = "gray", alpha = 0.5)
  }
  # Longitude lines (meridians, constant phi)
  phis <- seq(0, 2*pi, length.out = nlon + 1)
  for (ph in phis) {
    th <- seq(0, pi, length.out = 100)
    x <- sin(th) * cos(ph)
    y <- sin(th) * sin(ph)
    z <- cos(th)
    lines3d(x, y, z, col = "gray", alpha = 0.5)
  }
}

# Main function
analyze_and_visualize_spherical_data <- function(data_spherical) {
  if (!is.data.frame(data_spherical) || !all(c("theta", "phi") %in% colnames(data_spherical))) {
    stop("Input must be a data.frame with 'theta' and 'phi' columns.")
  }
  
  # Convert to Cartesian
  data_cartesian <- as.matrix(data.frame(
    x = sin(data_spherical$theta) * cos(data_spherical$phi),
    y = sin(data_spherical$theta) * sin(data_spherical$phi),
    z = cos(data_spherical$theta)
  ))
  n <- nrow(data_cartesian)
  if (n < 2) stop("Need at least 2 data points.")
  
  # Summary statistics
  cat("\n--- Summary Statistics ---\n")
  mean_dir <- colMeans(data_cartesian)
  R_bar <- sqrt(sum(mean_dir^2))
  kappa_hat <- estimate_kappa_3D(data_cartesian)
  mean_theta <- if (R_bar > .Machine$double.eps) acos(mean_dir[3]/R_bar) else NA
  mean_phi <- if (R_bar > .Machine$double.eps) atan2(mean_dir[2], mean_dir[1]) else NA
  cat("Mean Direction (Cartesian):", round(mean_dir, 4), "\n")
  cat("Mean Resultant Length (RL):", round(R_bar, 4), "\n")
  cat("Estimated Concentration (N:):", round(kappa_hat, 4), "\n")
  if (!is.na(mean_theta)) {
    cat(sprintf("Mean Direction (rad): N8 = %.4f, O = %.4f\n", mean_theta, mean_phi))
    cat(sprintf("Mean Direction (deg): N8 = %.2fB0, O = %.2fB0\n", mean_theta * 180/pi, mean_phi * 180/pi))
  } else {
    cat("Mean Direction (Spherical): Undefined (uniform data)\n")
  }

  s.dist <- as.matrix(1 - data_cartesian %*% t(data_cartesian))
  diag(s.dist) <- NA
  Q1_vals <- apply(s.dist, 1, min, na.rm = TRUE)
  E_vals <- sapply(1:n, function(i) {
    if (n - 1 < 1) return(NA)
    Dm <- data_cartesian[-i, , drop = FALSE]
    Rm <- sqrt(sum(colMeans(Dm)^2)) * (n - 1)
    Rg <- R_bar * n
    denom <- (n - 1 - Rm)
    if (denom == 0) return(NA)
    (n - 2) * (1 + Rm - Rg) / denom
  })
  R_minus <- sapply(1:n, function(i) {
    Dm <- if (n - 1 >= 1) data_cartesian[-i, , drop = FALSE] else NULL
    if (is.null(Dm) || nrow(Dm) < 1) return(NA)
    sqrt(sum(colMeans(Dm)^2)) * (n - 1)
  })
  C_vals <- if (!all(is.na(R_minus))) abs((R_minus / (n - 1) - R_bar)) / R_bar else rep(NA, n)
  
  cloud_outliers <- integer(0)
  if (!is.na(kappa_hat) && kappa_hat >= 0 && !is.infinite(kappa_hat)) {
    p50 <- find_percentile(0.50, kappa_hat)
    p75 <- find_percentile(0.75, kappa_hat)
    p99 <- find_percentile(0.99, kappa_hat)
    if (!anyNA(c(p50, p75, p99)) && p50 > 0) {
      estimated_k_cloud <- (p99 - p75) / p50
      med <- Directional::mediandir(data_cartesian)
      dists <- sapply(1:n, function(i) geodesic_distance(data_cartesian[i, ], as.numeric(med)))
   
      iqr_val <- quantile(dists, 0.50, na.rm = TRUE)
      Q3 <- quantile(dists, 0.75, na.rm = TRUE)
      whisker_upper <- Q3 + estimated_k_cloud * iqr_val
      cloud_outliers <- which(dists > whisker_upper)
      cat("cloud Method: Outliers at indices:", if (length(cloud_outliers)) paste(cloud_outliers, collapse = ", ") else "None", "\n")
      cat(sprintf("Estimated N: (cloud): %.4f, Whisker upper: %.4f, SMAD: %.4f, Q3: %.4f\n",
                  estimated_k_cloud, whisker_upper, iqr_val, Q3))
    }
  }
  cat("Spherical Median:", med, "\n")
  
  # Outlier Detection
  cat("\n--- Outlier Detection ---\n")
  cat("Q1 Method: Max stat at", which.max(Q1_vals), "(stat =", round(max(Q1_vals, na.rm = TRUE), 4), ")\n")
  if (!all(is.na(E_vals))) cat("E Method: Max stat at", which.max(E_vals), "(stat =", round(max(E_vals, na.rm = TRUE), 4), ")\n")
  if (!all(is.na(C_vals))) cat("C Method: Max stat at", which.max(C_vals), "(stat =", round(max(C_vals, na.rm = TRUE), 4), ")\n")
 
  # 3D Visualization
  cat("\n--- 3D Visualization ---\n")
  if (rgl.cur() != 0) close3d()
  open3d()
  spheres3d(0, 0, 0, radius = 1, col = "white", alpha = 0.1)
  add_spherical_grid()
  text3d(1.1, 0, 0, texts = "0", col = "blue", cex = 1)
  text3d(-1.1, 0, 0, texts = expression(pi), col = "blue", cex = 1)
  text3d(0, 1.1, 0, texts = expression(pi/2), col = "blue", cex = 1)
  text3d(0, -1.1, 0, texts = expression(3*pi/2), col = "blue", cex = 1)
  text3d(0, 0, 1.1, texts = "N", col = "blue", cex = 1.1)
  text3d(0, 0, -1.1, texts = "S", col = "blue", cex = 1.1)
  points3d(data_cartesian, col = "black", size = 10)
  
  for (i in 1:n) {
    pts <- data_cartesian[i, ]
    
    # small offset outward from the sphere
  
    label_pos <- pts + 0.05
    
    text3d(label_pos[1], label_pos[2], label_pos[3],
           texts = i, col = "black", cex = 1.1)
  }
  
  if (!is.na(R_bar)) {
    med <- Directional::mediandir(data_cartesian)
    points3d(matrix(med, nrow = 1), col = "yellow", size = 15)
  }
  
  if (length(cloud_outliers) > 0) {
    points3d(data_cartesian[cloud_outliers, , drop = FALSE], col = "red", size = 15)
  }
  
  if (!is.na(R_bar) && exists("whisker_upper") && exists("iqr_val")) {
    if (iqr_val > 0 && !anyNA(med)) plot_spherical_region(med, 0, iqr_val, col_shade = "#1f78b4", alpha_val = 0.4)
    if (whisker_upper > iqr_val && !anyNA(med)) plot_spherical_region(med, iqr_val, whisker_upper, col_shade = "#FFB6C1", alpha_val = 0.4)
  }
  
  for (i in 1:n) {
    pts <- data_cartesian[i, ]
    text3d(pts[1]- 0.5, pts[2], pts[3] , texts = i, col = "black", cex = 1.1)
  }
}



# A Graphical Goodness-of-Fit Analysis of Spherical Data to the von Mises-Fisher Distribution

n <- length(Datatheta)

# 2. Calculation of Sample Mean Direction (theta_bar, phi_bar)
L <- sum(sin(Datatheta) * cos(Dataphi))
M <- sum(sin(Datatheta) * sin(Dataphi))
N <- sum(cos(Datatheta))
R <- sqrt(L^2 + M^2 + N^2)

theta_bar <- acos(N / R)
phi_bar <- atan2(M, L)

# First Rotational Transformation to (theta_prime, phi_prime)
theta_prime <- numeric(n)
phi_prime <- numeric(n)
for (i in 1:n) {
  cos_theta_prime <- sin(Datatheta[i]) * sin(theta_bar) * cos(Dataphi[i] - phi_bar) + cos(Datatheta[i]) * cos(theta_bar)
  sin_theta_prime_cos_phi_prime <- sin(Datatheta[i]) * cos(theta_bar) * cos(Dataphi[i] - phi_bar) - cos(Datatheta[i]) * sin(theta_bar)
  sin_theta_prime_sin_phi_prime <- sin(Datatheta[i]) * sin(Dataphi[i] - phi_bar)
  theta_prime[i] <- acos(cos_theta_prime)
  phi_prime[i] <- atan2(sin_theta_prime_sin_phi_prime, sin_theta_prime_cos_phi_prime)
}

# 4. Second Rotational Transformation to (theta_prime_prime, phi_prime_prime)
theta_prime_prime <- numeric(n)
phi_prime_prime <- numeric(n)

for (i in 1:n) {
  cos_theta_pp <- sin(Datatheta[i]) * cos(theta_bar) * cos(Dataphi[i] - phi_bar) - cos(Datatheta[i]) * sin(theta_bar)
  sin_theta_pp_cos_phi_pp <- sin(Datatheta[i]) * sin(theta_bar) * cos(Dataphi[i] - phi_bar) + cos(Datatheta[i]) * cos(theta_bar)
  sin_theta_pp_sin_phi_pp <- -sin(Datatheta[i]) * sin(Dataphi[i] - phi_bar)
  
  theta_prime_prime[i] <- acos(cos_theta_pp)
  phi_prime_prime[i] <- atan2(sin_theta_pp_sin_phi_pp, sin_theta_pp_cos_phi_pp)
}
par(mfrow = c(1, 3))
# 5. Colatitude Plot Procedure (c' vs. Exponential Quantiles)
c_prime <- 1 - cos(theta_prime)
c_prime_sorted <- sort(c_prime)
exp_quantiles <- log(n / (n - (1:n) + 0.5))

plot(exp_quantiles, c_prime_sorted,
     main="(a) Colatitude Plot",
     xlab=expression(paste("Theoretical Exponential Quantiles, ", E[i,n])),
     ylab=expression(paste("Ordered Sample Values, ", c[i]^{"'"})))
lm_colat <- lm(c_prime_sorted ~ 0 + exp_quantiles)
abline(lm_colat, col="red")


# 6. Longitude Plot Procedure (phi' vs. Uniform Quantiles)
phi_prime_normalized <- phi_prime / (2 * pi)
phi_prime_normalized_sorted <- sort(phi_prime_normalized)
uniform_quantiles <- ((1:n) - 0.5) / n

plot(uniform_quantiles, phi_prime_normalized_sorted,
     main="(b) Longitude Plot",
     xlab=expression(paste("Theoretical Uniform Quantiles, ", U[i,n])),
     ylab=expression(paste("Ordered Normalized Sample Values, ", (phi[i]^{"'"} / 2*pi))))
abline(a=-0.5, b=1, col="red")

# 7. Two-Variable Plot Procedure (t vs. Normal Quantiles)
t_values <- phi_prime_prime * sqrt(sin(theta_prime_prime))
t_sorted <- sort(t_values)
normal_quantiles <- qnorm(((1:n) - 0.5) / n)

plot(normal_quantiles, t_sorted,
     main="(c) Two-Variable Plot",
     xlab=expression(paste("Theoretical Normal Quantiles, ", Q[i,n])),
     ylab=expression(paste("Ordered Sample Values, ", t[i])))
lm_two_var <- lm(t_sorted ~ 0 + normal_quantiles)
abline(lm_two_var, col="red")





## ===========================
## Choose dataset
## ===========================

dataset <- "eye"      # "eye" or "gait"

if (dataset == "eye") {
  
  library(readxl)
  
  my_data_df <- read_excel("eye_data.xlsx")
  
  
} else if (dataset == "gait") {
  
  library(fda)
  
  subject_data <- gait[1,,]
  
  my_data_df <- data.frame(
    theta = subject_data[,1] * pi/180,
    phi   = subject_data[,2] * pi/180
  )
}

analyze_and_visualize_spherical_data(my_data_df)
