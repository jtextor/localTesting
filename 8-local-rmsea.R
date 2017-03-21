# This script performs local test evaluation of latent variable models by
# testing single-implication path models.
#
# The script was added in the revised version of the accompanying article
# and it depends on the newest version of dagitty, available from github.

library(dagitty)
library(lavaan)
set.seed(1234)

if(packageVersion("dagitty") < "0.2.3") {
    stop("Please update the dagitty package to run this example by running
      devtools::install_github('jtextor/dagitty/r')")
}

# This function takes a SEM g with latent variables, and returns another
# SEM in which the structural model is saturated.
saturateStructure <- function(g){
	g <- as.dagitty(g)
	gm <- measurementPart(g)
	gs <- completeDAG(latents(g))
	c(gs,gm)
}

# This function takes a conditional independence implication x (in dagitty
# format) and a vector of variable names v, and returns a graph in which the
# only implied d-separation implication is x.
singleImplicationGraph <- function( x, v ){
	upper.part <- x$Z
	lower.part <- setdiff(v,c(x$X,x$Y,x$Z))
	r.edges <- ""
	for( i in upper.part ){
		r.edges <- paste(r.edges,i,"->",x$X,i,"->",x$Y)
		for( j in lower.part ){
			r.edges <- paste(r.edges,i,"->",j)
		}
	}
	for( i in lower.part ){
		r.edges <- paste(r.edges,i,"<->",x$X,i,"<->",x$Y)
	}
	if( length(upper.part)>1 ){
		ux <- combn(upper.part,2)
		r.edges <- paste(r.edges,paste(ux[1,],"<->",ux[2,]),collapse="\n")
	}
	if( length(lower.part)>1 ){
		ux <- combn(lower.part,2)
		r.edges <- paste(r.edges,paste(ux[1,],"<->",ux[2,]),collapse="\n")
	}
	as.dagitty(
		paste("dag{",r.edges,"}")
	)
}

# This function takes a SEM g (in dagitty syntax) and a data frame d, and evaluates 
# - the overall fit of g to d
# - the fit of the structural model only, given the measurement model
# - the fit of each individual implication of the structural model, given the
#   measurement model
evaluatePathModel <- function( g, d ){
	g <- as.dagitty(g)
	N <- nrow(d)
	m <- toString(g,"lavaan")
	m.SS <- toString(saturateStructure(g),"lavaan")
	m.SN <- toString(measurementPart(g),"lavaan")
	tst.m <- attributes(lavaan( m, d, auto.var=T, std.lv=T ))$test[[1]]
	tst.SS <- attributes(lavaan( m.SS, d, auto.var=T, std.lv=T ))$test[[1]]
	tst.SN <- attributes(lavaan( m.SN, d, auto.var=T, std.lv=T ))$test[[1]]
	r <- as.data.frame(rbind(
		(c(stat=(tst.m$stat),df=(tst.m$df))),
		(c(stat=(tst.m$stat - tst.SS$stat),df=(tst.m$df-tst.SS$df)))
	))
	rownames(r) <- c("Total","Structural")

	gs <- structuralPart(tested.model)
	gm <- measurementPart(tested.model)

	latents(gs)<-list()
	dseps <- impliedConditionalIndependencies(gs)
	for( i in seq_along(dseps) ){
		g.i <- c(singleImplicationGraph(dseps[[i]],names(gs)),gm)
		m.imp <- toString(g.i,"lavaan")
		tst.imp <- attributes(lavaan( m.imp, d, auto.var=T, std.lv=T ))$test[[1]]
		r.i <- data.frame(stat=(tst.imp$stat-tst.SS$stat),df=1)
		rownames(r.i) <- dseps[[i]]
		r <- rbind(r,r.i)
	}
	r$p.value <- apply( r, 1, 
		function(x) pchisq(x[1],x[2],lower.tail=FALSE) )
	r$RMSEA <- apply( r, 1, 
		function(x) sqrt(max(0,(x[1]/x[2]-1)/(N-1))) )
	r$C9 <- apply( r, 1, 
		function(x) (tst.SN$stat-x[1]-tst.SS$stat)/(tst.SN$stat-tst.SS$stat) )
	r$C10 <- apply( r, 1, 
		function(x) (x[1])/(tst.SN$stat-tst.SS$stat) )
	# C9 and C10 are path-based indices and as such undefined for the whole model.
	r[1,c("C9","C10")] <- NA
	r
}

# This graph is the same as in our previous partial mediation example.
# All variables are now latent variables, measured by 3 indicators each.
real.model <- dagitty("dag{
X[u]Y[u]M1[u]M2[u]M3[u]U1[u]
X -> U1 [beta=.4]
X -> M2 [beta=.4]
X -> M3 [beta=.4]
U1 -> Y [beta=.4]
U1 -> M1 [beta=.4]
M2 -> Y [beta=.4]
M3 -> Y [beta=.4]
M2 <-> M3 [beta=.4]

X -> {x1 x2 x3}
Y -> {x4 x5 x6}
M1 -> {x7 x8 x9}
M2 -> {x10 x11 x12}
M3 -> {x13 x14 x15}
}")

plot( graphLayout( real.model ) )

# The tested model has a correct measurement part, but an incorrect 
# structural part.
tested.model <- c(
	measurementPart( real.model ), 
	"dag{X[u]Y[u]M1[u]M2[u]M3[u]
	X -> { M1 M2 M3 } -> Y}"
)

# Run the analysis and print the results.
results <- evaluatePathModel(tested.model,simulateSEM(real.model,.8,.8))
print(signif(results,3))
