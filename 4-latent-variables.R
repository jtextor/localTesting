library(lavaan)
library(dagitty)
set.seed(1234)
N <- 1000

# Define the model in dagitty syntax. Variable and arrow attributes can be
# set in square brackets. We use this below to define which variables are
# latent and we also define some path coefficients that we will use to
# generate data
g.true <- dagitty("dag {
LX [latent]
LY [latent]
LX -> LY
LX -> { X1 X2 X3 }
LY -> { Y1 Y2 Y3 } 
X1 <-> X3 [beta=.25]
X1 <-> Y1 [beta=.25]
}")

#A simple graph of the model
coordinates(g.true) <- list(x=c(LX=1,LY=5,X1=0,X2=1,X3=2,Y1=4,Y2=5,Y3=6),
                            y=c(LX=-.5,LY=-.5,X1=0,X2=0,X3=0,Y1=0,Y2=0,Y3=0))
plot(g.true)

# Generate data, using .7 for all path coefficients not set in the syntax
d <- simulateSEM(g.true, 0.7, N = N)

# Lists all tetrad implications
vanishingTetrads(g.true)

# We assume a model that has more constraints (fewer arrows)
g.assumed <- dagitty("dag {
LX [latent] ; LY [latent]
{ X1 X2 X3 } <- LX -> LY -> {Y1 Y2 Y3}
}")

#A simple graph of the model
coordinates(g.assumed) <- list(x=c(LX=1,LY=5,X1=0,X2=1,X3=2,Y1=4,Y2=5,Y3=6),
                            y=c(LX=-.5,LY=-.5,X1=0,X2=0,X3=0,Y1=0,Y2=0,Y3=0))
plot(g.assumed)


# List all vanishing tetrads implied by the assumed model
print(vanishingTetrads(g.assumed))

# Convert the assumed model to lavaan syntax
m.assumed <- toString(g.assumed, "lavaan")

# Unrestricted model (EFA)
factanal(~X1+X2+X3+Y1+Y2+Y3,data=d,factors=2)

# Fit assumed model and request fit indices
fit <- sem(m.assumed, d, std.lv = TRUE)
print(summary(fit, fit = TRUE, mod = TRUE))
resid(fit,"standardized")

# Execute all tetrad tests and return p-values and confidence intervals
print(localTests(g.assumed, d, "tetrads"))

# Test of close fit for tetrads
print(localTests(g.assumed, d, "tetrads",tol=.05))

# Fit model in which largest modificiation index was added 
# (i.e., the covariance between X1 and Y1)
m.assumed2 <- paste(m.assumed,"\n","X1~~Y1")
fit2 <- sem(m.assumed2, d, std.lv = TRUE)
print(summary(fit2, fit = TRUE, mod = TRUE))
resid(fit2,"standardized")
