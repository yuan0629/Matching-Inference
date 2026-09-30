library(ggplot2)

#setwd("")

Iteration <- as.integer(Sys.getenv("SIM_ITERATIONS", unset = "500"))
d1 <- 100 # number of rows
K <- 5
d2 <- 150*K # number of columns
rreal <- 2 # real rank
Time <- 200*K # Time horizon (half); total sample size is 2*Time
rreal <- r <- 2 # rank used in practice
m <- 20
sigma <- 0.2
p0 <- 0.8 # probability used to construct the one-to-many target Q

eta <- 0.5 # step size
nu <- 1/d2
N0 <- Time/m/2

if (Iteration < 1L || N0 != as.integer(N0)) {
  stop("Iteration must be positive and Time must be divisible by 2*m.")
}
N0 <- as.integer(N0)
rows <- seq_len(d1)

set.seed(1124)

# low rank matrix
Mraw <- matrix(runif(d1*d2, -20, 20), nrow = d1, ncol = d2)
Mraw_svd <- svd(Mraw, nu = rreal, nv = rreal)
U0 <- Mraw_svd$u[, seq_len(rreal), drop = FALSE] %*%
  diag(sqrt(Mraw_svd$d[seq_len(rreal)]), nrow = rreal)
V0 <- Mraw_svd$v[, seq_len(rreal), drop = FALSE] %*%
  diag(sqrt(Mraw_svd$d[seq_len(rreal)]), nrow = rreal)
M <- U0 %*% t(V0)
lmin <- Mraw_svd$d[rreal]

# Store only the observed column and response for each row and time point.
# This is equivalent to the original dense X/Y lists, whose unobserved entries
# were all zero, but it avoids allocating thousands of d1-by-d2 matrices.
generate_oto_data <- function() {
  allcols <- matrix(0L, nrow = d1, ncol = Time)
  ally <- matrix(0, nrow = d1, ncol = Time)

  for (t in seq_len(Time)) {
    cols <- sample.int(d2, size = d1, replace = FALSE)
    allcols[, t] <- cols
    ally[, t] <- M[cbind(rows, cols)] + rnorm(d1, 0, sigma)
  }

  list(cols = allcols, y = ally)
}

observed_index <- function(cols) {
  rows + (cols - 1L)*d1
}

spectral_initialization <- function(allcols, ally, batch) {
  M0 <- matrix(0, nrow = d1, ncol = d2)
  scale <- 1/(N0*nu)

  for (i in batch) {
    index <- observed_index(allcols[, i])
    M0[index] <- M0[index] + scale*ally[, i]
  }

  M0_svd <- svd(M0, nu = r, nv = r)
  list(
    U = M0_svd$u[, seq_len(r), drop = FALSE],
    V = M0_svd$v[, seq_len(r), drop = FALSE]
  )
}

# Solve the same least-squares problem for G as the original observation-wise
# A/B accumulation, using a compact design matrix and BLAS cross-products.
fit_G <- function(U, V, allcols, ally, batch) {
  batch_cols <- as.vector(allcols[, batch, drop = FALSE])
  batch_y <- as.vector(ally[, batch, drop = FALSE])
  batch_rows <- rep.int(rows, times = length(batch))

  Ubatch <- U[batch_rows, , drop = FALSE]
  Vbatch <- V[batch_cols, , drop = FALSE]
  design <- matrix(0, nrow = length(batch_y), ncol = r^2)

  for (b in seq_len(r)) {
    block <- ((b - 1L)*r + 1L):(b*r)
    design[, block] <- Ubatch*Vbatch[, b]
  }

  Gvec <- solve(crossprod(design), crossprod(design, batch_y))
  matrix(Gvec, nrow = r, ncol = r)
}

# Accumulate the gradient only at observed entries. The original dense
# residual matrices were zero everywhere else.
gradient_sum <- function(Mcurrent, allcols, ally, batch) {
  sumXY <- matrix(0, nrow = d1, ncol = d2)

  for (i in batch) {
    index <- observed_index(allcols[, i])
    sumXY[index] <- sumXY[index] + Mcurrent[index] - ally[, i]
  }

  sumXY
}

estimate_initial <- function(allcols, ally) {
  error <- rep(0, m)

  init <- spectral_initialization(allcols, ally, seq_len(N0))
  U <- init$U
  V <- init$V
  G <- fit_G(U, V, allcols, ally, (N0 + 1L):(2L*N0))

  G_svd <- svd(G, nu = r, nv = r)
  LG <- G_svd$u
  Lambda_inv <- diag(1/G_svd$d, nrow = r, ncol = r)
  RG <- G_svd$v

  for (p in seq_len(m - 1L)) {
    gradient_batch <- (N0*(2L*p) + 1L):(N0*(2L*p + 1L))
    Mcurrent <- U %*% G %*% t(V)
    sumXY <- gradient_sum(Mcurrent, allcols, ally, gradient_batch)

    Um <- U %*% LG - eta/nu/N0*sumXY %*% V %*% RG %*% Lambda_inv
    Vm <- V %*% RG - eta/nu/N0*t(sumXY) %*% U %*% LG %*% Lambda_inv

    U <- svd(Um, nu = r, nv = 0)$u[, seq_len(r), drop = FALSE]
    V <- svd(Vm, nu = r, nv = 0)$u[, seq_len(r), drop = FALSE]

    regression_batch <- (N0*(2L*p + 1L) + 1L):(N0*2L*(p + 1L))
    G <- fit_G(U, V, allcols, ally, regression_batch)

    G_svd <- svd(G, nu = r, nv = r)
    LG <- G_svd$u
    Lambda_inv <- diag(1/G_svd$d, nrow = r, ncol = r)
    RG <- G_svd$v

    Minit <- U %*% G %*% t(V)
    error[p + 1L] <- log(sum((Minit - M)^2)/sum(M^2))
  }

  list(error = error, Minit = U %*% G %*% t(V))
}

gradient_Grassmannians <- function(data1, data2) {
  # Cross fitting: Minit1 is estimated on sample 2 and Minit2 on sample 1.
  result1 <- estimate_initial(data2$cols, data2$y)
  result2 <- estimate_initial(data1$cols, data1$y)

  list(
    error1 = result1$error,
    error2 = result2$error,
    Minit1 = result1$Minit,
    Minit2 = result2$Minit
  )
}

debias <- function(Minit, allcols, ally) {
  Munbs <- Minit
  scale <- 1/(Time*nu)

  for (t in seq_len(Time)) {
    index <- observed_index(allcols[, t])
    Munbs[index] <- Munbs[index] + scale*(ally[, t] - Minit[index])
  }

  Munbs
}

rank_r_projection <- function(Munbs) {
  Munbs_svd <- svd(Munbs, nu = r, nv = r)
  L <- Munbs_svd$u[, seq_len(r), drop = FALSE]
  R <- Munbs_svd$v[, seq_len(r), drop = FALSE]
  L %*% (t(L) %*% Munbs %*% R) %*% t(R)
}

Qlist <- list()

set.seed(1124)
# single entry
Q <- matrix(0, nrow = d1, ncol = d2)
Q[1, 1] <- 1
Qlist[[1]] <- Q

# one-to-one matching
Q <- matrix(0, nrow = d1, ncol = d2)
cols <- sample.int(d2, size = d1, replace = FALSE)
Q[cbind(rows, cols)] <- 1
Qlist[[2]] <- Q

# difference between two one-to-one matchings
Q1 <- matrix(0, nrow = d1, ncol = d2)
cols <- sample.int(d2, size = d1, replace = FALSE)
Q1[cbind(rows, cols)] <- 1

Q2 <- matrix(0, nrow = d1, ncol = d2)
cols <- sample.int(d2, size = d1, replace = FALSE)
Q2[cbind(rows, cols)] <- 1

Qlist[[3]] <- Q1 - Q2

# one-to-many target
num <- rbinom(d1, K, p0)
cols <- sample.int(d2, size = sum(num), replace = FALSE)
Q <- matrix(0, nrow = d1, ncol = d2)
Q[cbind(rep(rows, num), cols)] <- 1
Qlist[[4]] <- Q

value <- matrix(0, nrow = length(Qlist), ncol = Iteration)

for (iter in seq_len(Iteration)) {
  print(iter)

  data1 <- generate_oto_data()
  data2 <- generate_oto_data()

  result <- gradient_Grassmannians(data1, data2)
  Minit1 <- result$Minit1
  Minit2 <- result$Minit2

  Munbs1 <- debias(Minit1, data1$cols, data1$y)
  Munbs2 <- debias(Minit2, data2$cols, data2$y)

  Mproj1 <- rank_r_projection(Munbs1)
  Mproj2 <- rank_r_projection(Munbs2)
  Mfinal <- (Mproj1 + Mproj2)/2

  # Variance estimation. Each one-to-one observation contains exactly d1 entries.
  index1 <- observed_index(data1$cols)
  index2 <- observed_index(data2$cols)
  residual1 <- data1$y - matrix(Minit1[as.vector(index1)], nrow = d1, ncol = Time)
  residual2 <- data2$y - matrix(Minit2[as.vector(index2)], nrow = d1, ncol = Time)
  sigmahat <- (sum(residual1^2) + sum(residual2^2))/(2*Time*d1)

  Mfinal_svd <- svd(Mfinal, nu = r, nv = r)
  Lproj <- Mfinal_svd$u[, seq_len(r), drop = FALSE]
  Rproj <- Mfinal_svd$v[, seq_len(r), drop = FALSE]

  for (q in seq_along(Qlist)) {
    Q <- Qlist[[q]]
    S <- sum((t(Lproj) %*% Q)^2) +
      sum((Q %*% Rproj)^2) -
      sum((t(Lproj) %*% Q %*% Rproj)^2)

    value[q, iter] <- (sum(Mfinal*Q) - sum(M*Q))/
      sqrt(sigmahat*S/Time/2/nu)
  }
}

make_inference_plot <- function(x, trim_extreme = FALSE) {
  if (trim_extreme) {
    x <- x[abs(x) <= 5]
  }
  ggplot(data.frame(value = x), aes(x = value)) +
    geom_histogram(
      aes(y = after_stat(density)),
      fill = "#377eb8", color = "black", alpha = 0.6, stat = "bin", bins = 20
    ) +
    stat_function(fun = dnorm, colour = "red", linewidth = 1) +
    xlab("") +
    theme_classic()
}

inference_plots <- lapply(seq_along(Qlist), function(q) {
  make_inference_plot(value[q, ], trim_extreme = q == 1L)
})
fig1 <- inference_plots[[1]]
fig2 <- inference_plots[[2]]
fig3 <- inference_plots[[3]]
fig4 <- inference_plots[[4]]

for (q in seq_along(inference_plots)) {
  ggsave(
    filename = sprintf("density_oto_%d.png", q),
    plot = inference_plots[[q]], width = 5.5, height = 5, units = "in", dpi = 300
  )
}

value_oto <- value
save(value_oto, file = "value_oto.RData")



