jacobi <- function(A, x_0, b, tol = 1e-6) {
    n <- nrow(A)
    x <- list(x_0)
    k <- 1
    terminate <- FALSE
    while (!terminate) {
        x[[k + 1]] <- x[[k]]
        for (i in 1:n) {
            r <- 0
            for (j in 1:n) {
                if (i != j)
                    r <- r + A[i, j]*x[[k]][j]
            }
            x[[k + 1]][i] <- 1/A[i, i]*(b[i] - r)
        }
        if (norm(matrix(x[[k + 1]]) - matrix(x[[k]]), type = "F") < tol | k > 1e3)
            terminate <- T
        k <- k + 1
    }
    return(x[[k]])
}
