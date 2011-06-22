sim.spacodi <-
function(species=100, plots=30, missing.prop=0.6, mean.abund=15, sim.tree=FALSE,...){
	require(geiger)
	l.data=species*plots
	nulls=round(l.data*missing.prop)
	non.nulls=l.data-nulls
	m.log=log(mean.abund)
	out <- as.data.frame(matrix(sample(c(round(rlnorm(n=non.nulls, meanlog=m.log)),rep(0,nulls))),species,plots))
	rownames(out) <- paste("sp",seq(1:species),sep="")
	names(out) <- paste("plot",seq(1:plots),sep="")
	if(sim.tree) {
		if(species<2)stop("Simulation called for a tree with fewer than two taxa.")
		phy=prune.last.split(bd.tree(taxa.stop=(1+species),return.all.extinct=FALSE,...))
		phy$tip.label=row.names(out)
		res=list(out,phy)
		names(res)=c("sp.plot","sp.tree")
		return(res)
	} else {return(out)}
}

