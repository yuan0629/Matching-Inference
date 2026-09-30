library(ggplot2)

# Reproduce Appendix A: Table 3 and Figures 5--7.
# Outputs: one CSV table and 12 PDF figures (3 sample sizes x 4 query matrices).
#setwd("")

# Simulation settings
Iteration <- 500
d1 <- 100
K <- 5
d2 <- 150*K
rreal <- r <- 2
Tlist <- c(1000, 2000, 5000)
m <- 20
sigma <- 1
p0 <- 0.8
eta <- 0.7
nu <- 1/d2

# Generate the fixed rank-two signal matrix.
set.seed(1124)
M <- matrix(runif(d1*d2, -20, 20), nrow = d1, ncol = d2)
U0 <- svd(M)$u[, 1:rreal] %*% diag(c(sqrt(svd(M)$d)[1:rreal]))
V0 <- svd(M)$v[, 1:rreal] %*% diag(c(sqrt(svd(M)$d)[1:rreal]))
M <- U0 %*% t(V0)

# Estimate the low-rank matrix from sequential data batches.
estimate_initial <- function(alldataX, alldataY, alldatacode) {
  batch <- split(
    seq_along(alldataY),
    cut(seq_along(alldataY), breaks = 2*m, labels = FALSE)
  )

  M0 <- 0
  for (i in batch[[1]]) {
    M0 <- M0 + alldataY[[i]]/length(batch[[1]])/nu
  }
  U <- svd(M0)$u[, 1:r]
  V <- svd(M0)$v[, 1:r]

  fit_G <- function(U, V, batch_index) {
    A <- matrix(0, nrow = r^2, ncol = r^2)
    B <- matrix(0, nrow = r^2, ncol = 1)
    for (i in batch_index) {
      for (j in seq_len(nrow(alldatacode[[i]]))) {
        row_index <- alldatacode[[i]][j, 1]
        col_index <- alldatacode[[i]][j, 2]
        b <- matrix(
          c(matrix(t(U)[, row_index], ncol = 1) %*%
              matrix(V[col_index, ], nrow = 1)),
          ncol = 1
        )
        A <- A + b %*% t(b)
        B <- B + alldataY[[i]][row_index, col_index]*b
      }
    }
    matrix(solve(A) %*% B, nrow = r, ncol = r)
  }

  G <- fit_G(U, V, batch[[2]])
  G_svd <- svd(G)
  LG <- G_svd$u
  Lambda <- diag(G_svd$d)
  RG <- G_svd$v

  for (p in seq_len(m - 1)) {
    sumXY <- 0
    for (i in batch[[2*p + 1]]) {
      sumXY <- sumXY + (U %*% G %*% t(V))*alldataX[[i]] - alldataY[[i]]
    }

    Um <- U %*% LG - eta/nu/length(batch[[2*p + 1]]) *
      sumXY %*% V %*% RG %*% solve(Lambda)
    Vm <- V %*% RG - eta/nu/length(batch[[2*p + 1]]) *
      t(sumXY) %*% U %*% LG %*% solve(Lambda)

    U <- svd(Um)$u[, 1:r]
    V <- svd(Vm)$u[, 1:r]
    G <- fit_G(U, V, batch[[2*p + 2]])

    G_svd <- svd(G)
    LG <- G_svd$u
    Lambda <- diag(G_svd$d)
    RG <- G_svd$v
  }

  U %*% G %*% t(V)
}

# Four query matrices used in each row block of Table 3 and each figure panel.
Qlist <- list()
set.seed(1124)

# Single entry
Q <- matrix(0, nrow = d1, ncol = d2)
Q[1, 1] <- 1
Qlist[[1]] <- Q

# One-to-one matching
Q <- matrix(0, nrow = d1, ncol = d2)
cols <- sample(d2, size = d1, replace = FALSE)
Q[cbind(seq_len(d1), cols)] <- 1
Qlist[[2]] <- Q

# Difference of two one-to-one matchings
Q1 <- matrix(0, nrow = d1, ncol = d2)
cols <- sample(d2, size = d1, replace = FALSE)
Q1[cbind(seq_len(d1), cols)] <- 1

Q2 <- matrix(0, nrow = d1, ncol = d2)
cols <- sample(d2, size = d1, replace = FALSE)
Q2[cbind(seq_len(d1), cols)] <- 1
Qlist[[3]] <- Q1 - Q2

# One-to-many matching
num <- rbinom(d1, K, p0)
cols <- sample(d2, size = sum(num), replace = FALSE)
Q <- matrix(0, nrow = d1, ncol = d2)
Q[cbind(rep(seq_len(d1), num), cols)] <- 1
Qlist[[4]] <- Q

Qlabels <- c(
  "$Q=e_1e_1^{\\top}$",
  "$Q\\in\\mathcal{M}^{\\mathrm{oto}}$",
  "$Q=Q_1-Q_2$",
  "$Q\\in\\mathcal{M}_{K,p_0}^{\\mathrm{otm}}$"
)

# Standardized statistics produce the figures; raw estimates produce Table 3.
array_dim <- c(length(Qlist), Iteration, length(Tlist))
value_split <- array(NA_real_, dim = array_dim)
value_nosplit <- array(NA_real_, dim = array_dim)
theta_split <- array(NA_real_, dim = array_dim)
theta_nosplit <- array(NA_real_, dim = array_dim)

# Use the same simulated dataset for the splitting and no-splitting estimators.
for (tt in seq_along(Tlist)) {
  T <- Tlist[tt]
  Time <- T/2

  for (iter in seq_len(Iteration)) {
    message("T=", T, ", replication=", iter, "/", Iteration)
    set.seed(as.integer(1124 + d1*10000000 + d2*10000 + T + iter))

    alldataX <- alldataY <- alldatacode <- vector("list", T)
    for (t in seq_len(T)) {
      cols <- sample(d2, size = d1, replace = FALSE)
      X <- Y <- matrix(0, nrow = d1, ncol = d2)
      X[cbind(seq_len(d1), cols)] <- 1
      Y[cbind(seq_len(d1), cols)] <-
        M[cbind(seq_len(d1), cols)] + rnorm(d1, 0, sigma)
      alldataX[[t]] <- X
      alldataY[[t]] <- Y
      alldatacode[[t]] <- cbind(seq_len(d1), cols)
    }

    first_half <- seq_len(Time)
    second_half <- (Time + 1):T

    # Sample splitting: estimate on one half and debias on the other half.
    Minit1 <- estimate_initial(
      alldataX[second_half], alldataY[second_half], alldatacode[second_half]
    )
    Minit2 <- estimate_initial(
      alldataX[first_half], alldataY[first_half], alldatacode[first_half]
    )

    Munbs1 <- Minit1
    for (t in first_half) {
      Munbs1 <- Munbs1 +
        (alldataY[[t]] - Minit1*alldataX[[t]])/Time/nu
    }
    Munbs2 <- Minit2
    for (t in second_half) {
      Munbs2 <- Munbs2 +
        (alldataY[[t]] - Minit2*alldataX[[t]])/Time/nu
    }

    L1 <- svd(Munbs1)$u[, 1:r]
    R1 <- svd(Munbs1)$v[, 1:r]
    Mproj1 <- L1 %*% t(L1) %*% Munbs1 %*% R1 %*% t(R1)

    L2 <- svd(Munbs2)$u[, 1:r]
    R2 <- svd(Munbs2)$v[, 1:r]
    Mproj2 <- L2 %*% t(L2) %*% Munbs2 %*% R2 %*% t(R2)
    Mfinal <- (Mproj1 + Mproj2)/2

    sigmahat <- 0
    for (t in first_half) {
      sigmahat <- sigmahat +
        sum((alldataY[[t]] - alldataX[[t]]*Minit1)^2)/Time/2/sum(alldataX[[t]])
    }
    for (t in second_half) {
      sigmahat <- sigmahat +
        sum((alldataY[[t]] - alldataX[[t]]*Minit2)^2)/Time/2/sum(alldataX[[t]])
    }

    Lproj <- svd(Mfinal)$u[, 1:r]
    Rproj <- svd(Mfinal)$v[, 1:r]
    for (q in seq_along(Qlist)) {
      Q <- Qlist[[q]]
      S <- sum((t(Lproj) %*% Q)^2) +
        sum((Q %*% Rproj)^2) -
        sum((t(Lproj) %*% Q %*% Rproj)^2)
      theta_split[q, iter, tt] <- sum(Mfinal*Q)
      value_split[q, iter, tt] <-
        (theta_split[q, iter, tt] - sum(M*Q))/sqrt(sigmahat*S/Time/2/nu)
    }

    # No sample splitting: estimate and debias using the full dataset.
    Minit <- estimate_initial(alldataX, alldataY, alldatacode)
    Munb <- Minit
    for (t in seq_len(T)) {
      Munb <- Munb + (alldataY[[t]] - Minit*alldataX[[t]])/T/nu
    }

    L <- svd(Munb)$u[, 1:r]
    R <- svd(Munb)$v[, 1:r]
    Mproj <- L %*% t(L) %*% Munb %*% R %*% t(R)

    sigmahat <- 0
    for (t in seq_len(T)) {
      sigmahat <- sigmahat +
        sum((alldataY[[t]] - alldataX[[t]]*Minit)^2)/T/sum(alldataX[[t]])
    }

    Lproj <- svd(Mproj)$u[, 1:r]
    Rproj <- svd(Mproj)$v[, 1:r]
    for (q in seq_along(Qlist)) {
      Q <- Qlist[[q]]
      S <- sum((t(Lproj) %*% Q)^2) +
        sum((Q %*% Rproj)^2) -
        sum((t(Lproj) %*% Q %*% Rproj)^2)
      theta_nosplit[q, iter, tt] <- sum(Mproj*Q)
      value_nosplit[q, iter, tt] <-
        (theta_nosplit[q, iter, tt] - sum(M*Q))/sqrt(sigmahat*S/T/nu)
    }
  }
}

# Save the five columns reported in Table 3.
format_sig5 <- function(x) {
  formatC(x, format = "fg", digits = 5, flag = "#")
}

table_lines <- "T,Q,Var. (splitting),Var. (no splitting),Variance ratio"
for (tt in seq_along(Tlist)) {
  for (q in seq_along(Qlist)) {
    variance_split <- var(theta_split[q, , tt])
    variance_nosplit <- var(theta_nosplit[q, , tt])
    q_label <- Qlabels[q]
    if (q == 4) q_label <- paste0('"', q_label, '"')
    table_lines <- c(
      table_lines,
      paste(
        Tlist[tt], q_label,
        format_sig5(variance_split),
        format_sig5(variance_nosplit),
        format_sig5(variance_split/variance_nosplit),
        sep = ","
      )
    )
  }
}
writeLines(table_lines, "app_variance_table.csv")

# Save the 12 title-free density figures.
for (tt in seq_along(Tlist)) {
  for (q in seq_along(Qlist)) {
    value1df <- data.frame(
      value1 = c(value_split[q, , tt], value_nosplit[q, , tt]),
      method = rep(
        c("With sample splitting", "Without sample splitting"),
        each = Iteration
      )
    )
    breaks <- seq(min(value1df$value1), max(value1df$value1), length.out = 21)
    figure <- ggplot(value1df, aes(x = value1)) +
      geom_histogram(
        aes(y = after_stat(density)),
        fill = "#377eb8", color = "black", alpha = 0.6,
        stat = "bin", breaks = breaks
      ) +
      stat_function(fun = dnorm, colour = "red", linewidth = 1) +
      xlab("") +
      theme_classic(base_size = 11) +
      facet_wrap(~method, nrow = 1, scales = "fixed") +
      theme(
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 11),
        strip.text = element_text(size = 11),
        strip.background = element_blank(),
        panel.spacing = grid::unit(12, "pt"),
        plot.margin = margin(5.5, 5.5, 5.5, 5.5)
      )

    ggsave(
      sprintf("app_T%d_Q%d.pdf", Tlist[tt], q),
      figure, units = "px", width = 2054, height = 930, device = "pdf"
    )
  }
}
