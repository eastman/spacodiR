r.plot <-
function(species=100, plots=30, missing.prop=0.6, mean.abund=15, sim.tree=FALSE){
	require(geiger)
	l.data=species*plots
	nulls=round(l.data*missing.prop)+1
	non.nulls=l.data-nulls
	m.log=log(mean.abund)
	out <- as.data.frame(matrix(sample(c(round(rlnorm(n=non.nulls, meanlog=m.log)),rep(0,nulls))),species,plots))
	rownames(out) <- paste("sp",seq(1:species),sep="")
	names(out) <- paste("plot",seq(1:plots),sep="")
	if(sim.tree) {
		phy=prune.last.split(bd.tree(b=0.02,d=0.0,taxa.stop=(1+species),return.all.extinct=FALSE))
		phy$tip.label=row.names(out)
		return(list(out,phy))
	} else {return(out)}
}

