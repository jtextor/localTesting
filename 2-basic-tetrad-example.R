# This example shows how to test tetrad constraints implied by a structural
# equation model using both a parametric test by Wishart and by bootstrapping
# standard errors. The code only shows how to test the first implied tetrad.
# We only use the standard R package 'boot' for bootstrapping

library(boot)
set.seed(1234)

# Sample size
N <- 1000

# Standardized coefficient for all paths
p <- 0.3

# Generate data according to two-factor latent variable model, where X,Y,Z are
# indicators of the first factor, and W is the single indicator of the second
# factor
U1 <- rnorm(N)
U2 <- p * U1 + rnorm(N, 0, sqrt(1 - p^2))
X <- p * U1 + rnorm(N, 0, sqrt(1 - p^2))
Y <- p * U1 + rnorm(N, 0, sqrt(1 - p^2))
Z <- p * U1 + rnorm(N, 0, sqrt(1 - p^2))
W <- p * U2 + rnorm(N, 0, sqrt(1 - p^2))

# Compute sample covariance matrix
d <- data.frame(X, Y, Z, W)
S <- cov(d)

# Print value of tetrad WYZX. Using the paper's notation, I={W,X} and J={Y,Z}
tetrad <- det(S[c("W", "X"), c("Y", "Z")])

# Determine standard error of tetrad WYZX using Wishart's formula
d.IJ <- det(S[c("W", "Y", "Z", "X"), c("W", "Y", "Z", "X")])
d.I <- det(S[c("W", "X"), c("W", "X")])
d.J <- det(S[c("Y", "Z"), c("Y", "Z")])

tetrad.se <- sqrt((d.I * d.J * (N + 1)/(N - 1) - d.IJ)/(N - 2))
(tetrad.pval <- 2*pnorm(abs(tetrad/tetrad.se), lower.tail = FALSE))

#test of close fit, here against arbitrary value of .05
tetrad.z <- tetrad / tetrad.se 
tol.z <- .05 / tetrad.se 
tetrad.pval.close <- pchisq( tetrad.z^2, 1, ncp=tol.z^2, lower.tail=FALSE )

# Bootstrap standard error of tetrad WYZX
boot.out <- boot(d, function(d, i) det(cov(d[i, ])[c("W", "X"), c("Y", "Z")]), 
    R = 1000)
print(boot.ci(boot.out, type = "norm"))

#Remaining tetrads are computed via dagitty
ex2model <- dagitty("dag {
                    U1 [latent] ; U2 [latent]
                     { X Y Z } <- U1 -> U2 -> {W}
                     }")

#Wishart's test
(localTests(ex2model, data.frame(X,Y,Z,W), "tetrads"))
#Bootstrapped
(localTests(ex2model, data.frame(X,Y,Z,W), "tetrads",R=1000))
#close fit
(localTests(ex2model, data.frame(X,Y,Z,W), "tetrads",tol=.05))

