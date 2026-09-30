library(ggplot2)

setwd("")

Iteration <- as.integer(Sys.getenv("SIM_ITERATIONS", unset = "500"))
d1 <- 100 # number of rows
K <- 5
d2 <- 150*K # number of columns
p0 <- 0.8 # probability used only to construct the fourth target Q
p1 <- 0.9
p2 <- 0.9
rreal <- 2 # real rank
Time <- 200*K # Time horizon (half); total sample size is 2*Time
rreal <- r <- 2 # rank used in practice
m <- 20
sigma <- 0.2
gamma <- 0.2
cr <- cs <- 0.8

eta <- 0.5 # step size
N0 <- Time/m/2

if (Iteration < 1L || N0 != as.integer(N0)) {
  stop("Iteration must be positive and Time must be divisible by 2*m.")
}
N0 <- as.integer(N0)
rows <- seq_len(d1)

set.seed(1124)

# Draw the same truncated-binomial matching size as in the original code.
draw_matching_size <- function() {
  repeat {
    repeat {
      Br <- rbinom(1L, d1, p1)
      if (Br >= d1*cr) break
    }
    repeat {
      Bs <- rbinom(1L, d2, p2)
      if (Bs >= d2*cs) break
    }
    if (Br >= (1 + gamma)*Bs || Bs >= (1 + gamma)*Br) break
  }
  min(Br, Bs)
}

# Sampling probability. Retaining this Monte Carlo calculation and its RNG
# order keeps the definition of nu identical to the original script.
nu <- 0
for (i in seq_len(1000L)) {
  nu <- nu + draw_matching_size()
}
nu <- nu/1000/d1/d2

# Low-rank matrix. One truncated SVD replaces the repeated full SVD calls.
Mraw <- matrix(runif(d1*d2, -20, 20), nrow = d1, ncol = d2)
Mraw_svd <- svd(Mraw, nu = rreal, nv = rreal)
U0 <- Mraw_svd$u[, seq_len(rreal), drop = FALSE] %*%
  diag(sqrt(Mraw_svd$d[seq_len(rreal)]), nrow = rreal)
V0 <- Mraw_svd$v[, seq_len(rreal), drop = FALSE] %*%
  diag(sqrt(Mraw_svd$d[seq_len(rreal)]), nrow = rreal)
M <- U0 %*% t(V0)
lmin <- Mraw_svd$d[rreal]

# Store only observed coordinates and responses instead of 2*Time dense
# d1-by-d2 X/Y matrices. Each observation is still generated in the same order.
generate_tside_data <- function() {
  allrows <- vector("list", Time)
  allcols <- vector("list", Time)
  ally <- vector("list", Time)

  for (t in seq_len(Time)) {
    Bmin <- draw_matching_size()
    observed_rows <- sample.int(d1, size = Bmin, replace = FALSE)
    cols <- sample.int(d2, size = Bmin, replace = FALSE)

    allrows[[t]] <- observed_rows
    allcols[[t]] <- cols
    ally[[t]] <- M[cbind(observed_rows, cols)] + rnorm(Bmin, 0, sigma)
  }

  list(rows = allrows, cols = allcols, y = ally)
}

observed_index <- function(observed_rows, cols) {
  observed_rows + (cols - 1L)*d1
}

spectral_initialization <- function(data, batch) {
  M0 <- matrix(0, nrow = d1, ncol = d2)
  scale <- 1/(N0*nu)

  for (i in batch) {
    index <- observed_index(data$rows[[i]], data$cols[[i]])
    M0[index] <- M0[index] + scale*data$y[[i]]
  }

  M0_svd <- svd(M0, nu = r, nv = r)
  list(
    U = M0_svd$u[, seq_len(r), drop = FALSE],
    V = M0_svd$v[, seq_len(r), drop = FALSE]
  )
}

# Solve the original observation-wise least-squares problem for G using one
# compact design matrix and BLAS cross-products.
fit_G <- function(U, V, data, batch) {
  batch_rows <- unlist(data$rows[batch], use.names = FALSE)
  batch_cols <- unlist(data$cols[batch], use.names = FALSE)
  batch_y <- unlist(data$y[batch], use.names = FALSE)

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

# The original residual matrices are zero outside observed coordinates.
gradient_sum <- function(Mcurrent, data, batch) {
  sumXY <- matrix(0, nrow = d1, ncol = d2)

  for (i in batch) {
    index <- observed_index(data$rows[[i]], data$cols[[i]])
    sumXY[index] <- sumXY[index] + Mcurrent[index] - data$y[[i]]
  }

  sumXY
}

estimate_initial <- function(data) {
  error <- rep(0, m)

  init <- spectral_initialization(data, seq_len(N0))
  U <- init$U
  V <- init$V
  G <- fit_G(U, V, data, (N0 + 1L):(2L*N0))

  G_svd <- svd(G, nu = r, nv = r)
  LG <- G_svd$u
  Lambda_inv <- diag(1/G_svd$d, nrow = r, ncol = r)
  RG <- G_svd$v

  for (p in seq_len(m - 1L)) {
    gradient_batch <- (N0*(2L*p) + 1L):(N0*(2L*p + 1L))
    Mcurrent <- U %*% G %*% t(V)
    sumXY <- gradient_sum(Mcurrent, data, gradient_batch)

    Um <- U %*% LG - eta/nu/N0*sumXY %*% V %*% RG %*% Lambda_inv
    Vm <- V %*% RG - eta/nu/N0*t(sumXY) %*% U %*% LG %*% Lambda_inv

    U <- svd(Um, nu = r, nv = 0)$u[, seq_len(r), drop = FALSE]
    V <- svd(Vm, nu = r, nv = 0)$u[, seq_len(r), drop = FALSE]

    regression_batch <- (N0*(2L*p + 1L) + 1L):(N0*2L*(p + 1L))
    G <- fit_G(U, V, data, regression_batch)

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
  result1 <- estimate_initial(data2)
  result2 <- estimate_initial(data1)

  list(
    error1 = result1$error,
    error2 = result2$error,
    Minit1 = result1$Minit,
    Minit2 = result2$Minit
  )
}

debias <- function(Minit, data) {
  Munbs <- Minit
  scale <- 1/(Time*nu)

  for (t in seq_len(Time)) {
    index <- observed_index(data$rows[[t]], data$cols[[t]])
    Munbs[index] <- Munbs[index] + scale*(data$y[[t]] - Minit[index])
  }

  Munbs
}

rank_r_projection <- function(Munbs) {
  Munbs_svd <- svd(Munbs, nu = r, nv = r)
  L <- Munbs_svd$u[, seq_len(r), drop = FALSE]
  R <- Munbs_svd$v[, seq_len(r), drop = FALSE]
  L %*% (t(L) %*% Munbs %*% R) %*% t(R)
}

estimate_sigma2 <- function(Minit1, Minit2, data1, data2) {
  sigmahat <- 0

  for (t in seq_len(Time)) {
    index <- observed_index(data1$rows[[t]], data1$cols[[t]])
    residual <- data1$y[[t]] - Minit1[index]
    sigmahat <- sigmahat + sum(residual^2)/(2*Time*length(residual))
  }
  for (t in seq_len(Time)) {
    index <- observed_index(data2$rows[[t]], data2$cols[[t]])
    residual <- data2$y[[t]] - Minit2[index]
    sigmahat <- sigmahat + sum(residual^2)/(2*Time*length(residual))
  }

  sigmahat
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
  if (iter == 1L || iter %% 10L == 0L) message("replication ", iter, "/", Iteration)

  data1 <- generate_tside_data()
  data2 <- generate_tside_data()

  result <- gradient_Grassmannians(data1, data2)
  Minit1 <- result$Minit1
  Minit2 <- result$Minit2

  Munbs1 <- debias(Minit1, data1)
  Munbs2 <- debias(Minit2, data2)

  Mproj1 <- rank_r_projection(Munbs1)
  Mproj2 <- rank_r_projection(Munbs2)
  Mfinal <- (Mproj1 + Mproj2)/2

  sigmahat <- estimate_sigma2(Minit1, Minit2, data1, data2)

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

value_tside <- value
save(value_tside, file = "value_tside.RData")

make_inference_plot <- function(x) {
  ggplot(data.frame(value = x), aes(x = value)) +
    geom_histogram(
      aes(y = after_stat(density)),
      fill = "#377eb8", color = "black", alpha = 0.6, stat = "bin", bins = 20
    ) +
    stat_function(fun = dnorm, colour = "red", linewidth = 1) +
    xlab("") +
    theme_classic()
}

# Keep the four Q statistics in their construction order and draw each one in
# its own file without a title.
inference_plots <- lapply(seq_along(Qlist), function(q) {
  make_inference_plot(value_tside[q, ])
})
fig1 <- inference_plots[[1]]
fig2 <- inference_plots[[2]]
fig3 <- inference_plots[[3]]
fig4 <- inference_plots[[4]]

for (q in seq_along(inference_plots)) {
  ggsave(
    filename = sprintf("density_tside_%d.png", q),
    plot = inference_plots[[q]], width = 5.5, height = 5, units = "in", dpi = 300
  )
}
