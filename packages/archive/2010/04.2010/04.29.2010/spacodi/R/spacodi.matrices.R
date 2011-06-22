spacodi.matrices<-function(sp.plot, phy=NULL, sp.traits=NULL) {
	if((!missing(phy) && !missing(sp.traits)) || (missing(phy) && missing(sp.traits))) stop("Please supply either a phylogeny or a trait matrix, but not both.")

	orig=names(sp.plot)
	l.spp=nrow(sp.plot)
	drop.plots=vector()
	for(sp in 1:ncol(sp.plot)) {
		l.nulls=length(which(sp.plot[,sp]==0))
		if((l.spp-l.nulls)<2) {
			drop.plots=cbind(drop.plots, sp)
		}
	}
	
	if(length(drop.plots)!=0) {
		sp.plot=sp.plot[,-drop.plots[!is.na(drop.plots)]]
		warning("At least one plot was pruned from the matrix, there being fewer than two species sampled therein.")
	}
		
	plots=names(sp.plot)
	l.plots=length(plots)
	sp.array=array(dim=c(l.plots, l.plots, 4))
	out.array=array(dim=c(l.plots, 1, 4))
	r.out=list()
	
	ticker=c(1:40)*ceiling(l.plots/40)
	cat("\npairwise SPACoDi calculations in progress:\n")
	
	for(row in 1:l.plots) {
		rr=plots[row]
		r=which(plots==rr)
		for(cc in plots[which(plots!=rr)]) {
			c=which(plots==cc)
			keep=plots[!is.na(match(plots, c(cc,rr)))]
			sp.rc=sp.plot[,keep]
			if(!missing(phy)) {
				if(any(attr(phy, "names")=="node.label"))phy$node.label=NULL
				ss=spacodi.calc(sp.plot=sp.rc, phy=phy)
			} else if(!missing(sp.traits) && (length(intersect(row.names(sp.traits), r.spl<-row.names(sp.plot)))==length(r.spl)) ) {
				ss=spacodi.calc(sp.plot=sp.rc, sp.traits=sp.traits, all.together=TRUE)
			} else {stop("Poor structure in input object.")}
			class(ss)="vector"
			for(xx in 1:length(ss)) {
				sp.array[r,c,xx]=unlist(unname(ss[which(names(ss)==names(ss)[xx])]))
			}
		}
		for(yy in 1:4) {
			foo=sp.array[r,,yy]
			out.array[r,1,yy]=sum(foo[which(!is.na(foo))])/(l.plots-1)
		}
		
		if(row%in%ticker)cat(".")
	
	}
	cat("\n")
	for(zz in 1:4) {
		foo=as.data.frame(out.array[,,zz])
		row.names(foo)=plots
		d.matrix=dist(foo)
		r.out[[zz]]=d.matrix
	}
	if(!missing(phy)) {
		names(r.out)=c("Ist","Pst","Bst","PIst") 
	} else {
		names(r.out)=c("Ist","Tst","T*st","TAUst")
		if(ncol(sp.traits)>1) message("All traits were considered together for distance matrices")
	}
	return(r.out)
}