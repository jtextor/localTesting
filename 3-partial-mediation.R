# Load required packages.
library(lavaan)
library(dagitty)
# At least version 0.2.2 of the dagitty package is required for this and all
# following examples.
if(packageVersion("dagitty") < "0.2.2") {
    stop("Please update the dagitty package to run this example!")
}
set.seed(1234)
N <- 1000

# Below we define a graphical model in dagitty. Models are defined using
# arrow operators like ->, <-, or <->. Variables can be grouped using curly
# braces to shorten notation
g.true <- dagitty("dag {
U1 [latent]
X -> { U1 M2 <-> M3 } -> Y
U1 -> M1
}")

# A simple graph can be generated automatically after X and Y coordinates 
# for each variable have been supplied
coordinates(g.true) <- list(x=c(X=0,U1=1,M1=2,M2=2,M3=2,Y=3),
                            y=c(X=0,U1=-1,M1=-1,M2=0,M3=1,Y=0))
plot(g.true)

# Generate data according to the true model. All standardized path
# coefficients are set to 0.9
d <- simulateSEM(g.true, 0.4, N = N)

# Enumerate all d-separation constraints implied by the true model
print(impliedConditionalIndependencies(g.true))

# Define a second graphical model that is slightly misspecified
g.assumed <- dagitty("dag { X -> { M1 M2 M3 } -> Y }")

#Again a simple graph
coordinates(g.assumed) <- list(x=c(X=0,M1=2,M2=2,M3=2,Y=3),
                            y=c(X=0,M1=-1,M2=0,M3=1,Y=0))
plot(g.assumed)


# Enumerate all d-separation constraints implied by the assumed model
print(impliedConditionalIndependencies(g.assumed))

# Perform statistical tests of conditional independence and report results
print(localTests(g.assumed, d, "cis",tol=0))

#Also perform tests of close fit using .05 as tolerage
print(localTests(g.assumed, d, "cis",tol=.05))


# Convert the model to lavaan syntax
mymodel <- toString(g.assumed, "lavaan")

# Estimate the model in lavaan and request fit statistics and standardized 
# residuals
fit <- sem(mymodel, d)
print(summary(fit, fit = TRUE, modindices = TRUE))
resid(fit,type="standardized")


# Add a covariance between M2 and M3, as suggested by the modification index
g.assumed.2 <- dagitty("dag { X -> { M1 M2 <-> M3 } -> Y }")
mymodel.2 <- toString(g.assumed.2, "lavaan")
fit.2 <- sem(mymodel.2, d)
print(summary(fit.2, fit = TRUE, modindices = TRUE))
resid(fit.2,type="standardized")
