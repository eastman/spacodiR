spacodi.matrices<-function(sp.plot, phy=NULL, sp.traits=NULL) {
	if((!is.null(phy) && !is.null(sp.traits)) || (is.null(phy) && is.null(sp.traits))) stop("Please supply either a phylogeny or a trait matrix, but not both.")
	Bst.flag=FALSE
	if(!is.null(phy)) {
		sp.data=match.spacodi.data(sp.plot=sp.plot, phy=phy)
		sp.plot=sp.data$sp.plot
		phy=sp.data$sp.tree
		Bst.flag=TRUE
	}	
	if(!is.null(sp.traits)) {
		sp.data=match.spacodi.data(sp.plot=sp.plot, sp.traits=sp.traits)
		sp.plot=sp.data$sp.plot
		sp.traits=sp.data$sp.traits		
	}
	sp.plot=as.data.frame(sp.plot)
	plots=names(sp.plot)
	l.plots=length(plots)
	sp.array=array(dim=c(l.plots, l.plots, 4))
	out.array=array(dim=c(l.plots, 1, 4))
	r.out=list()
	
	ticker=c(1:40)*ceiling(l.plots/40)
	cat("\npairwise SPACoDi calculations between plots in progress:\n")
	
	for(row in 1:l.plots) {
		rr=plots[row]
		r=which(plots==rr)
		for(cc in plots[which(plots!=rr)]) {
			c=which(plots==cc)
			keep=plots[!is.na(match(plots, c(cc,rr)))]
			sp.rc=sp.plot[,keep]
			if(!is.null(phy)) {
				if(any(attr(phy, "names")=="node.label"))phy$node.label=NULL
				ss=suppressWarnings(spacodi.calc(sp.plot=sp.rc, phy=phy))
			} else if(!is.null(sp.traits) && (length(intersect(row.names(sp.traits), r.spl<-row.names(sp.plot)))==length(r.spl)) ) {
				ss=suppressWarnings(spacodi.calc(sp.plot=sp.rc, sp.traits=sp.traits, all.together=TRUE))
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
		r.out[[zz]]=dist(foo)
	}
	if(Bst.flag) {
		names(r.out)=c("Ist","Pst","Bst","PIst") 
	} else {
		names(r.out)=c("Ist","Tst","T*st","TAUst")
		if(ncol(sp.traits)>1) message("All traits were considered together for distance matrices")
	}
	return(r.out)
}