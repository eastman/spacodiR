randomization.test=function(obs=obs, exp=exp, iter=10000, verbose=FALSE, two.tailed=TRUE){
	if(verbose && length(obs)<iter) warning("iterations exceeds observed data")
	if(verbose) warning(paste("p-value roughly approximated for ",length(exp), " values\n",sep=""))
	
	O=sample(as.numeric(obs[!is.na(obs)]), iter, replace=TRUE)
	E=sample(as.numeric(exp[!is.na(exp)]), iter, replace=TRUE)
	
	rr=cbind((O-E), NA)
	if(any(is.na(rr[,1])))stop("Cannot deal with NaN values")
	rr[,2]=ifelse(rr[,1]>0, 1, ifelse(rr[,1]<0, 0, sample(c(0,1),1)))
	p=rr[,2]
	
	p.g=1-(sum(p)/iter)
	p.l=1-p.g
	
	names(p.g)="greater"
	names(p.l)="lesser"
	
	if(!two.tailed) {
		res=list(p.g,p.l)
		names(res)=c("greater","lesser")
	} else {
		if(p.g!=p.l)res=min(c(p.g,p.l))*2 else res=1
		names(res)="p.value"
	}
	if(verbose) {
		res=list(c(O-E), res)
		names(res)=c("differences", "greater", "lesser")
	}
	return(res)
}
