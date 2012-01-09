r.sp_plot <-
function(obj,n.rep=10){
	sp_plot=obj
	zz=unlist(unname(obj))
	nulls=length(which(zz==0))
	l.data=length(zz)
	non.nulls=l.data-nulls
	yy=zz[which(zz!=0)]
	m.log=mean(log(yy))
	out <- list() ### prep an output file
	for(j in 1:n.rep){
		foo <- as.data.frame(matrix(sample(c(round(rlnorm(n=non.nulls, meanlog=m.log)),rep(0,nulls))),nrow(sp_plot),ncol(sp_plot)))
		rownames(foo) <- row.names(sp_plot)
		names(foo)=names(sp_plot)
		out[[j]]=foo
	}
	return(out)
}

