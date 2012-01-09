r.sp.plot <-
function(obj,n.rep=10){
	sp.plot=obj
	zz=unlist(unname(obj))
	nulls=length(which(zz==0))
	l.data=length(zz)
	non.nulls=l.data-nulls
	yy=zz[which(zz!=0)]
	m.log=mean(log(yy))
	out <- list() ### prep an output file
	for(j in 1:n.rep){
		foo <- as.data.frame(matrix(sample(c(round(rlnorm(n=non.nulls, meanlog=m.log)),rep(0,nulls))),nrow(sp.plot),ncol(sp.plot)))
		rownames(foo) <- row.names(sp.plot)
		names(foo)=names(sp.plot)
		out[[j]]=foo
	}
	return(out)
}

