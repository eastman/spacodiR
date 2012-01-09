resamp.3t <-
function(obj, dmat=NULL) {
	if(is.null(dmat)) {
		dmat=as.data.frame(matrix(0,ncol(obj),ncol(obj)))
		row.names(dmat)=names(obj)
		names(dmat)=names(obj)
	}
	if(!is.null(dmat) && all(row.names(dmat)%in%names(dmat)) && row.names(dmat)%in%names(obj)) {
		dmat=dmat-min(dmat)
		if(all(dmat==0)) flag=TRUE else flag=FALSE
		dd=stack(dmat)
		dd=cbind(dd,row.names(dmat))
		names(dd)=c("dist","p.from","p.to")
		dd=dd[-which(dd[,2]==dd[,3]),]
		p.all=split(dd, dd$p.from)
		p.plots=array(dim=c(ncol(obj)))
		for(pp in 1:ncol(obj)) {
			p.plots[pp]=mean(p.all[[pp]][,1])
		}
		p.plots=as.data.frame(p.plots)
		p.plots=(p.plots-min(p.plots))/(max(p.plots)-min(p.plots))
		row.names(p.plots)=names(p.all)
	} else {stop("Distance matrix cannot be interpreted with supplied dataset")}
	
	for(ss in 1:nrow(obj)){
		if(!flag) {
			p.dists=unique(c(unname(p.plots))[[1]])
			r=runif(1)
			r.d=max(p.dists[which(r >= p.dists)])
			shuffle=which(p.plots==r.d)
			if(length(shuffle)>1) shuffle=sample(shuffle,1)
		} else {shuffle = sample(1:ncol(obj),1)}
		
		torus=rep(1:ncol(obj),2)

		t.array=array(dim=c(ncol(obj),2))
		t.array[,1]=1:ncol(obj)
		for(o in 1:ncol(obj)){
			tt=torus[shuffle+(o-1)]
			t.array[tt,2]=obj[ss,o]
		}
		obj[ss,]=t.array[,2]
	}
	if(flag || is.null(dmat)) warning("Plots were assumed to be equidistant from one another.")
	return(obj)
}

