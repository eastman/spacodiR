resamp.3t <-
function(obj, dmat=NULL) {
	if(is.null(dmat))flag=TRUE else flag=FALSE
	names.orig=names(obj)
	names(obj)=seq(1:ncol(obj))
	if(is.null(dmat)) {
		dmat=as.data.frame(matrix(0,ncol(obj),ncol(obj)))
		names(dmat)=names.orig
		row.names(dmat)=names.orig
	}
	dmat=as.data.frame(as.matrix(dmat))
	if(!all(names(dmat)%in%names.orig) || ncol(obj)!=ncol(dmat) || ncol(dmat)!=nrow(dmat))stop("Names in distance matrix do not correspond to plot names")
	row.names(dmat)=names(obj)
	names(dmat)=names(obj)
	
# find all distances from plot.tt to plot.tt + some shifter value (e.g., '3' would be plot1 to plot4, plot10 to plot3, ... plotN to plotN+3)
# tabulate these values and find the average distance from each plot to every plot+'shifter'
	torus=rep(1:ncol(obj),2)
	plus.array=array(dim=c(1,ncol(obj)))
	torus.array=array(dim=c(ncol(obj), ncol(obj)))
	for(plus in 1:ncol(obj)) {
		for(tt in 1:ncol(torus.array)) {
			from=tt
			to=torus[tt+plus]
			d.tt=dmat[from, to]
			torus.array[tt,plus]=d.tt
		}
		plus.array[1,plus]=mean(torus.array[,plus],na.rm=TRUE)
	}

	plus.array=as.data.frame(plus.array)
	names(plus.array)=names(obj)
	plus.array=plus.array[,order(plus.array)]

# randomly generate a value between 0 and maximum average distance from a plot to every other plot+'shifter'
# shift species abundances by the randomly chosen 'shifter' 
	for(ss in 1:nrow(obj)){
		if(!flag) {
			shifter=as.numeric(names(plus.array))[min(which(plus.array>=runif(1,min=min(plus.array), max=max(plus.array))))]
		} else {shifter = sample(as.numeric(names(plus.array)),1)}
				
		t.array=array(dim=c(ncol(obj),2))
		t.array[,1]=1:ncol(obj)
		for(o in 1:ncol(obj)){
			tt=torus[shifter+(o-1)]
			t.array[tt,2]=obj[ss,o]
		}
		obj[ss,]=t.array[,2]
	}
	if(flag) message("Plots were assumed to be equidistant from one another.")
	res=obj
	names(res)=names.orig
	return(res)
}

