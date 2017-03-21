# This example shows how to simulate data from a standardized structural
# equation model and test a d-separation constraint, without using any R
# packages

set.seed(1234)

# Sample size
N <- 1000

# Path coefficient for all paths
p <- 0.3

# The bi-directed arrow between a and b is replaced by a structure a <- L -> b
# for the simulation
L <- rnorm(N, 0, 1)

# We now generate data according to the model structure, setting the residual
# variance such that the expected variance of each variable is 1
A <- sqrt(p) * L + rnorm(N, 0, sqrt(1 - sqrt(p)^2))
B <- p * A + rnorm(N, 0, sqrt(1 - p^2))
C <- sqrt(p) * L + rnorm(N, 0, sqrt(1 - sqrt(p)^2))
D <- p * B + p * C + rnorm(N, 0, sqrt(1 - p^2 - p^2 - 2 * p^4))

# Now we test the implication Cov( B, C | A )=0. First, we regress both B and
# C on A, and compute the residuals
rB.A <- lm(B ~ A)$resid
rC.A <- lm(C ~ A)$resid

# Now, we test if there is any correlation between those residuals
cor.test(rB.A, rC.A)

#Then we perform a test of not-close fit, testing observed value against .05
z <- atanh(cor(rB.A, rC.A))
sigmafz <-  1/sqrt(length(A) - 3 -2)
pval <- pchisq((z/sigmafz)^2,1, ncp =  (atanh(.05)/sigmafz)^2,lower.tail = FALSE)

# Next we test the implication Cov( A, D | B, C )=0 in the same manner
rA.BC <- lm(A ~ B + C)$resid
rD.BC <- lm(D ~ B + C)$resid
cor.test(rA.BC, rD.BC)

#Then we perform a test of not-close fit, testing observed value against .05
z2 <- atanh(cor(rA.BC, rD.BC))
sigmafz2 <-  1/sqrt(length(A) - 3 -2)
pval2 <- pchisq((z2/sigmafz)^2,1, ncp =  (atanh(.05)/sigmafz)^2,lower.tail = FALSE)
