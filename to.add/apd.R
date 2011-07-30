apd=function(sp.plot, cophen, abundance=TRUE){
	PDij=cophen[upper.tri(cophen)]
	spl=sp.plot
	N=unique(dim(cophen))
	if(length(dim)!=1) stop("cophenetic matrix appears improperly formatted")
	if(!all(rownames(cophen)==rownames(spl))){
		spl=as.matrix(spl[match(rownames(spl),rownames(cophen)),],ncol=ncol(spl))
		spl=as.data.frame(spl)
		rownames(spl)=rownames(cophen)
	}
	DeltaP=sum(2*PDij)/(N*(N-1))
	if(abundance){
		relAbund=rowSums(spl)/sum(spl)
	} else {
		relAbund=apply(spl,1,function(x)sum(x!=0))/sum(spl!=0)
	}
	ffPDij<-rAij<-matrix(0,N,N)	
	for(r in 1:(N-1)){
		for(c in (r+1):N){
			rAij[r,c]=relAbund[r]*relAbund[c]
			ffPDij[r,c]=rAij[r,c]*cophen[r,c]
		}	
	}
	Db=sum(ffPDij)/sum(rAij)
	apd=(DeltaP-Db)/DeltaP
	return(round(apd,12))
}
