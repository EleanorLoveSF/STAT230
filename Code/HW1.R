# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Dustin Pluta
# Assignment 1
# STAT 230: Winter 2017
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

setwd("~/Dropbox/Coursework/Winter2017/STAT230")
dat1 <- read.csv("Data/Assignment1.csv", row.names = 1)

X <- as.matrix(dat1[, 1:5])
y <- as.matrix(dat1$y)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## PROBLEM 1 ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# (a) Gram-Schmidt
source("Code/gram_schmidt.R")
QR <- gram_schmidt(X)
beta <- solve(QR$R) %*% t(QR$Q) %*% y
print(beta)

# (b) Householder
source("Code/householder.R")
QR <- householder(X)
beta <- solve(QR$R) %*% t(QR$Q) %*% y
print(beta)

# (c) Jacobi
source("Code/jacobi.R")
beta <- jacobi(t(X) %*% X, rep(1, 5), t(X) %*% y)
print(beta)

# Verify values using lm()
fit <- lm(y ~ . - 1, data = dat1)
print(fit$coefficients)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## PROBLEM 2 ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Set BATCH 1 data and compute coefficient
# estimates for using the jacobi method.
X_n <- X[1:80, ]
y_n <- y[1:80]
beta_batch1 <- jacobi(t(X_n) %*% X_n,
               rep(1, 5), t(X_n) %*% y_n)

# Set BATCH 2 data and update coefficient
# estimates from BATCH 1
# using formula from Lecture 3.
X_k <- X[81:100, ]
y_k <- y[81:100]
A <- t(X) %*% X
b <- A %*% beta_batch1 + t(X_k) %*% (y_k - X_k %*% beta_batch1)
beta <- jacobi(A, beta_batch1, b)


