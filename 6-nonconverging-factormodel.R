library(lavaan)
library(dagitty)
set.seed(123)
N <- 100

# We define a single-factor model with residual correlations between the
# first two and the last two indicators
g.true <- dagitty("dag{ 
L [latent]
L -> { X1 <-> X2 [beta=.25] }
L -> { X3 <-> X4 [beta=.25] }
L -> X3 [beta=.1]
}")

#A simple graph of the model
coordinates(g.true) <- list(x=c(L=0,X1=-1,X2=-.5,X3=.5,X4=1),
                            y=c(L=0,X1=1,X2=1,X3=1,X4=1))
plot(g.true)

# Simulate data with all unspecified loadings set to .8
d <- simulateSEM(g.true, 0.8, N = N)

# We estimate the model by imposing an equality constraint between the
# loadings of X1 and X3. Without such a constraint, the model would not be
# identified
m.assumed <- sem("u=~l*X1+X2+l*X3+X4
X1~~X2
X3~~X4", d)

# The model does not converge. That could be due to the model structure or
# the equality constraint
print(summary(m.assumed, fit = TRUE))

# We can use local tests to test the structure in isolation. This shows that
# the structure is OK
print(localTests(g.true, d, type = "tetrads"))
print(localTests(g.true, d, type = "tetrads",tol=.05))

# So the problem must be the equality constraint. We try a different one:
m.assumed.2 <- sem("u=~l*X1+X2+X3+l*X4
X1~~X2
X3~~X4", d)

# This model does fit
print(summary(m.assumed.2, fit = TRUE))
