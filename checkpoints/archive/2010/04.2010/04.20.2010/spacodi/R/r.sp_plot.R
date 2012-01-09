r.sp_plot <-
function(obj,n.rep=10){
	sp_plot=obj
	nulls=length(which(sp_plot==0))
	l.data=ncol(sp_plot)*nrow(sp_plot)
	non.nulls=l.data-nulls
	m.log=mean(log(sum(sp_plot)/non.nulls))
	out <- list() ### prep an output file
	for(j in 1:n.rep){
		out[[j]] <- as.data.frame(matrix(sample(c(round(rlnorm(n=non.nulls, meanlog=m.log)),rep(0,nulls))),nrow(sp_plot),ncol(sp_plot)))
		rownames(out[[j]]) <- row.names(sp_plot)
		names(out[[j]])=names(sp_plot)
	}
	return(out)
}

