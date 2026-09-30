# Reproduce Table 1 in Section 5.2 of MatchingInference.pdf.
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(file_arg)) {
  setwd(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
}
if (!requireNamespace("clue", quietly = TRUE)) {
  stop("Package 'clue' is required.")
}

data <- read.csv(
  "clean data.csv", stringsAsFactors = FALSE,
  na.strings = c("", "NA"), check.names = FALSE
)
MLBdata <- data[data$sample == "analysis_2025", ]

# Convert the original zero-based identifiers to consecutive R indices.
row_levels <- sort(unique(MLBdata$row_index))
col_levels <- sort(unique(MLBdata$col_index))
MLBdata$row <- match(MLBdata$row_index, row_levels)
MLBdata$col <- match(MLBdata$col_index, col_levels)

d1 <- length(row_levels)
d2 <- length(col_levels)
T_total <- length(unique(MLBdata$time_index))
T_half <- T_total/2L
r <- 2L
m <- 3L
eta <- 0.2
N0 <- T_half/(2*m)
nu <- nrow(MLBdata)/(T_total*d1*d2)
split_seed <- 324L

stopifnot(
  nrow(MLBdata) == 2073L,
  d1 == 30L, d2 == 75L, T_total == 180L,
  T_half == 90L, N0 == 15L,
  !anyDuplicated(MLBdata[c("time_index", "row")]),
  !anyDuplicated(MLBdata[c("time_index", "col")])
)

# Construct daily observation matrices X_t and reward matrices Y_t.
build_matrices <- function(date_order) {
  X <- Y <- code <- vector("list", length(date_order))
  for (t in seq_along(date_order)) {
    daily <- MLBdata[MLBdata$time_index == date_order[t], ]
    rows <- daily$row
    cols <- daily$col
    if (anyDuplicated(rows) || anyDuplicated(cols)) {
      stop("A daily matrix violates the row/column at-most-one constraint.")
    }
    X[[t]] <- Y[[t]] <- matrix(0, d1, d2)
    X[[t]][cbind(rows, cols)] <- 1
    Y[[t]][cbind(rows, cols)] <- daily$accuracy_above_expected
    code[[t]] <- cbind(rows, cols)
  }
  list(X = X, Y = Y, code = code)
}

# Algorithm 1, using only the identification half D1.
estimate_initial <- function(D) {
  stopifnot(length(D$X) == 2*m*N0)

  M0 <- matrix(0, d1, d2)
  for (t in seq_len(N0)) M0 <- M0 + D$Y[[t]]/(N0*nu)
  U <- svd(M0)$u[, seq_len(r), drop = FALSE]
  V <- svd(M0)$v[, seq_len(r), drop = FALSE]

  fit_G <- function(batch, U, V) {
    A <- matrix(0, r^2, r^2)
    B <- matrix(0, r^2, 1)
    for (t in batch) {
      for (k in seq_len(nrow(D$code[[t]]))) {
        i <- D$code[[t]][k, 1]
        j <- D$code[[t]][k, 2]
        z <- matrix(c(
          matrix(t(U)[, i], ncol = 1) %*%
            matrix(V[j, ], nrow = 1)
        ), ncol = 1)
        A <- A + z %*% t(z)
        B <- B + D$Y[[t]][i, j]*z
      }
    }
    matrix(solve(A) %*% B, r, r)
  }

  G <- fit_G((N0 + 1):(2*N0), U, V)
  for (p in seq_len(m - 1L)) {
    G_svd <- svd(G)
    L_G <- G_svd$u
    Lambda <- diag(G_svd$d, r, r)
    R_G <- G_svd$v

    residual <- matrix(0, d1, d2)
    gradient_batch <- (N0*(2*p) + 1):(N0*(2*p + 1))
    for (t in gradient_batch) {
      residual <- residual +
        (U %*% G %*% t(V))*D$X[[t]] - D$Y[[t]]
    }

    U_update <- U %*% L_G -
      eta/(nu*N0)*residual %*% V %*% R_G %*% solve(Lambda)
    V_update <- V %*% R_G -
      eta/(nu*N0)*t(residual) %*% U %*% L_G %*% solve(Lambda)
    U <- svd(U_update)$u[, seq_len(r), drop = FALSE]
    V <- svd(V_update)$u[, seq_len(r), drop = FALSE]

    regression_batch <- (N0*(2*p + 1) + 1):(2*N0*(p + 1))
    G <- fit_G(regression_batch, U, V)
  }
  U %*% G %*% t(V)
}

make_matching <- function(rows, cols) {
  Q <- matrix(0, d1, d2)
  Q[cbind(rows, cols)] <- 1
  Q
}

optimal_matching <- function(M, rows) {
  score <- M[rows, , drop = FALSE]
  score <- score - min(score) + 1
  cols <- as.integer(clue::solve_LSAP(score, maximum = TRUE))
  make_matching(rows, cols)
}

# Fixed 90/90 split used in the paper.
set.seed(split_seed)
date_order <- sample(sort(unique(MLBdata$time_index)), T_total)
D1_dates <- date_order[seq_len(T_half)]
D2_dates <- date_order[T_half + seq_len(T_half)]
D1 <- build_matrices(D1_dates)
D2 <- build_matrices(D2_dates)
M_initial <- estimate_initial(D1)

# Debias with D2 and project to rank r.
M_debiased <- M_initial
for (t in seq_len(T_half)) {
  M_debiased <- M_debiased +
    (D2$Y[[t]] - M_initial*D2$X[[t]])/(T_half*nu)
}
s <- svd(M_debiased)
L <- s$u[, seq_len(r), drop = FALSE]
R <- s$v[, seq_len(r), drop = FALSE]
M_projected <- L %*% t(L) %*% M_debiased %*% R %*% t(R)

sigma_hat <- mean(vapply(seq_len(T_half), function(t) {
  sum((D2$Y[[t]] - D2$X[[t]]*M_initial)^2)/sum(D2$X[[t]])
}, numeric(1)))

projected_svd <- svd(M_projected)
L_projected <- projected_svd$u[, seq_len(r), drop = FALSE]
R_projected <- projected_svd$v[, seq_len(r), drop = FALSE]

# One-sided test of <M,Q_optimal-Q_observed> > 0.
test_contrast <- function(Q) {
  tangent_norm_sq <-
    sum((t(L_projected) %*% Q)^2) +
    sum((Q %*% R_projected)^2) -
    sum((t(L_projected) %*% Q %*% R_projected)^2)
  estimate <- sum(M_projected*Q)
  standard_error <- sqrt(sigma_hat*tangent_norm_sq/(T_half*nu))
  z_value <- estimate/standard_error
  data.frame(
    estimate = estimate,
    p_value = pnorm(z_value, lower.tail = FALSE)
  )
}

# All 10 dates in D1 with exactly 15 arrived home teams.
D1_counts <- table(MLBdata$time_index[MLBdata$time_index %in% D1_dates])
observed_days <- sort(as.integer(names(D1_counts)[D1_counts == 15L]))
stopifnot(length(observed_days) == 10L)

table1 <- do.call(rbind, lapply(observed_days, function(time) {
  daily <- MLBdata[MLBdata$time_index == time, ]
  rows <- daily$row
  Q_observed <- make_matching(rows, daily$col)
  Q_optimal <- optimal_matching(M_initial, rows)
  cbind(
    data.frame(date = as.character(unique(daily$date))),
    test_contrast(Q_optimal - Q_observed)
  )
}))
rownames(table1) <- NULL

# Verify and save the exact values reported in Table 1.
stopifnot(
  sum(table1$p_value < 0.05) == 8L,
  isTRUE(all.equal(
    round(table1$estimate, 3),
    c(21.085, 22.422, 21.652, 20.820, 17.702,
      10.517, 3.740, 0.203, 12.417, 11.380)
  )),
  isTRUE(all.equal(
    round(table1$p_value, 4),
    c(0.0001, 0.0002, 0.0002, 0.0008, 0.0001,
      0.0490, 0.2247, 0.4849, 0.0449, 0.0175)
  ))
)
write.csv(table1, "Table 1 - optimal versus observed.csv", row.names = FALSE)

print(transform(
  table1,
  estimate = sprintf("%.3f", estimate),
  p_value = sprintf("%.4f", p_value)
), row.names = FALSE)
