householder <- function(X) {
    n <- nrow(X)
    p <- ncol(X)
    U <- matrix(0, nrow=n, ncol=p)
    Q <- diag(nrow=n)
    for (k in 1:p) {
        w <-  matrix(X[k:n, k])
        w[1] <- w[1] - norm(w, type="F")
        u <-  w/norm(w, type="F")
        U[k:n, k] <- u
        X[k:n, k:p] <- X[k:n, k:p] - 2*u %*% (t(u) %*% as.matrix(X[k:n, k:p]))
        Q <- Q %*% (diag(n) - 2*U[, k] %*% t(U[, k]))
    }
    R <- matrix(0, nrow=p, ncol=p)
    R[which(upper.tri(X, diag=T), arr.ind=T)] <- X[which(upper.tri(X, diag=T),
                                                         arr.ind=T)]
    return(list(U=U, Q=Q[, 1:p], R=R))
}
