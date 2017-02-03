gram_schmidt <- function(X) {
  n <- nrow(X)
  n_vars <- ncol(X)
  Q <- matrix(0, nrow=n, ncol=n_vars)
  R <- matrix(0, nrow=n_vars, ncol=n_vars)
  v <- as.matrix(X[, 1])
  R[1, 1] <- norm(v, type="F")
  Q[, 1] <- v/R[1, 1]
  for (j in 2:n_vars) {
    v <- X[, j]
    for (i in 1:(j - 1)) {
      R[i, j] <- t(Q[, i]) %*% X[, j]
      v <- v - R[i, j] %*% Q[, i]
    }
    R[j, j] <- norm(v, type="F")
    Q[, j] <- v/R[j, j]
  }
  return(list(Q=Q, R=R))
}
