library(lavaan)
library(dagitty)
set.seed(123)
N <- 1000

# We generate non-linear data using an 'additive noise model'. The noise is
# still Gaussian, but the variables are non-linearly dependent
X <- 4 + rnorm(N)
Y <- X^2 + 0.2 * sd(X^2) * rnorm(N)
Z <- Y^2 + 0.2 * sd(Y^2) * rnorm(N)

# Scale data to variance 1 for numerical reasons
d <- as.data.frame(scale(cbind(X, Y, Z)))

# A standard full mediation model does not fit to this data. Modification
# indices suggest adding a direct effect X -> Z
m <- sem("Z~Y\nY~X", d)
print(summary(m, fit = TRUE, mod = TRUE))


# The (linear) local test also fails
print(localTests("dag{X->Y->Z}", type = "cis", data = d))

# However, the semi-parametric local test indicates that conditional
# independence does hold (the confidence interval includes 0). Thus, the
# misfit is due to data distribution rather than model structure
print(localTests("dag{X->Y->Z}", type = "cis.loess", data = d, R = 500))

# Visualize linear regressions & residual correlations
par(mfrow = c(2, 3))
lmX.Y <- lm(X ~ Y, d)
lmZ.Y <- lm(Z ~ Y, d)
with(d, plot(Y, X))
abline(lmX.Y, col = 2)
with(d, plot(Y, Z))
abline(lmZ.Y, col = 2)
scatter.smooth(lmX.Y$resid, lmZ.Y$resid, span = 5, lpars = list(col = 2))

# Visualize non-linear smoothing using loess
lsX.Y <- loess(X ~ Y, d)
lsZ.Y <- loess(Z ~ Y, d)
with(d, scatter.smooth(Y, X, lpars = list(col = 2)))
with(d, scatter.smooth(Y, Z, lpars = list(col = 2)))
scatter.smooth(lsX.Y$resid, lsZ.Y$resid, span = 5, lpars = list(col = 2))
