randomization.test.sp=function(obs=obs, exp=exp, iter=10000, return.all=FALSE, two.tailed=TRUE){
	
	obs=obs[is.finite(as.numeric(obs))]
	exp=exp[is.finite(as.numeric(exp))]
	
	O=sample(as.numeric(obs), size=iter, replace=TRUE)
	E=sample(as.numeric(exp), size=iter, replace=TRUE)
	
	rr=cbind((O-E), NA)
	rr[,2]=ifelse(rr[,1]>0, 1, ifelse(rr[,1]<0, 0, sample(c(0,1),1)))
	p=rr[,2]
	
	p.g=1-(sum(p)/iter)
	p.l=1-p.g
	
	out.d=as.numeric(rr[,1])
	
	if(!two.tailed) {
		out.p=list(p.g,p.l)
	} else {
		if(p.g!=p.l)out.p=min(c(p.g,p.l))*2 else out.p=1
	}
	
	if(return.all) {
		r.out=list(out.d, out.p)
	} else {
		r.out=out.p
	}
	if(length(r.out)==1) {									#r.all=FALSE; two=TRUE
		names(r.out)="p.value"
	} else if(two.tailed==FALSE){
		if(return.all==TRUE) {
			names(r.out)=c("differences","p.value")         #r.all=TRUE; two=FALSE
			names(r.out[[2]])=c("greater","lesser")
		} else {
			names(r.out)=c("greater","lesser")				#r.all=FALSE; two=FALSE
		}	
	} else {												#r.all=TRUE; two=TRUE
		names(r.out)=c("differences","p.value")  
	}
	return(r.out)
}
