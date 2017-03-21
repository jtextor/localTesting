library(lavaan)
library(dagitty)
set.seed(123)
N <- 1000

# True model
g.true <- dagitty("dag{ 
J -> I [beta=-.8] 
I -> X [beta=.4]
X <-> M [beta=.25]
J -> X -> M -> Y
J -> Y
}")
d <- simulateSEM(g.true, 0.5, N = N)

#A simple graph of the model
coordinates(g.true) <- list(x=c(I=1,J=0,X=2,M=3,Y=2),
                            y=c(I=1,J=0,X=0,M=0,Y=-1))
plot(g.true)


# Assumed model; lacks the arrow I -> X
g.assumed <- dagitty("dag{ J -> Y ; J -> I -> X -> M -> Y ; X <-> M }")
m.assumed <- toString(g.assumed, "lavaan")

#A simple graph of the model
coordinates(g.assumed) <- list(x=c(I=1,J=0,X=2,M=3,Y=2),
                            y=c(I=1,J=0,X=0,M=0,Y=-1))
plot(g.assumed)



# The model fails to converge. Hence we cannot list modification indices
fit <- sem(m.assumed, d)
print(summary(fit))

# The most significant local test suggests adding an arrow I -> X to the
# model (an arrow X -> I would yield a cyclic model)
print(localTests(g.assumed, d,tol=0))
print(localTests(g.assumed, d,tol=.05))

# Adding this arrow gives the true model, which passes local tests ...
print(localTests(g.true, d))

# ... and also also converges and fits globally
m.true <- toString(g.true, "lavaan")
fit.true <- sem(m.true, d)
print(summary(fit.true))
