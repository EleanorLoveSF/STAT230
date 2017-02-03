householder <- function(X) {
    n <- nrow(X)
    p <- ncol(X)
    U <- matrix(0, nrow=n, ncol=p)
    Q <- diag(nrow=n)
    for (k in 1:p) {
        w = matrix(X[k:n, k])
        print(w)
        w[1] = w[1] - norm(w, type="F")
        u = w/norm(w, type="F")
        U[k:n, k] = u
        X[k:n, k:p] = X[k:n, k:p] - 2*u %*% (t(u) %*% X[k:n, k:p])
        Q <- Q %*% (diag(n) - 2*U[, k] %*% t(U[, k]))
    }
    R <- matrix(0, nrow=p, ncol=p)
    R[which(upper.tri(X, diag=T), arr.ind=T)] <- X[which(upper.tri(X, diag=T), arr.ind=T)]
    return(list(U=U, R=R))
}

## ####
# householder <- function(A) {
#     m <- nrow(A)
#     n <- ncol(A)
#     Q <- diag(nrow = m, ncol = m)
#     R <- A
#     for (j in 1:n) {
#         normx <- norm(matrix(R[j:m, j]))
#         s <- -sign(R[j, j])
#         u1 <- R[j, j] - s*normx
#         w <- R[j:m, j]/u1
#         w[1] <- 1
#         tau <- s*u1/normx
#     }
#     R[j:m, ] <- R[j:m, ] - (tau*w) %*% (t(w) %*% R[j:m, ])
#     print(Q[, j:m], w)
#     Q[, j:m] <- Q[, j:m] - (Q[, j:m] %*% w) %*% t(tau*w)
#     return(list(Q, R))
# }

## ####
X <- matrix(c(1:3, 5, 7, 8), nrow = 3)
QR <- householder(X)

## ####
y <- X %*% c(1, -1)

## ####
Qy <- y
for (j in 1:ncol(QR[[1]])) {
    u <- QR[[1]][ ,j]
    print(u)
    Qy <- (diag(nrow(X)) - 2*u %*% t(u)) %*% Qy
}
Qy
u <- QR[[1]][, 1]
H1 <- diag(length(u)) - 2*u %*% t(u)
H1 %*% H1

H2 <- diag(length(u)) - 2*QR[[1]][, 2] %*% t(QR[[1]][, 2])
p <- 2
solve(QR$R) %*% (H2 %*% H1 %*% y)[1:p]

