match.spacodi.data<-function(sp.plot, phy=NULL, sp.traits=NULL, prune=TRUE, verbose=FALSE) {
	if(is.null(row.names(sp.plot))) stop("Check that sp.plot has both row names as species.") 
	if(is.null(colnames(sp.plot))) {
		warning("sp.plot does not appear to have plots as column names.")
		names(sp.plot)=paste("plot",1:ncol(sp.plot),sep="")
	}
	sp.plot=as.matrix(sp.plot)
	if(any(!is.finite(sp.plot))) stop("Poor data values in the sp.plot: NA or NaN.")
	
# find missing species
	missing.spp=vector()
	if(!missing(phy)) {
		if(class(phy)!="phylo") stop("Tree does not appear to be of class 'phylo'")
		if(!is.null(phy$node.label)) phy$node.label=NULL
		missing.spp=c(row.names(sp.plot)[!row.names(sp.plot)%in%phy$tip.label],missing.spp)
	}
	if(!missing(sp.traits)) {
		missing.spp=c(row.names(sp.plot)[!row.names(sp.plot)%in%row.names(sp.traits)],missing.spp)
		drop.traits=row.names(sp.traits)%in%row.names(sp.plot)
	}

# find undersampled plots
	if(prune) {
		l.spp=nrow(sp.plot)
		drop.plots=vector()
		for(sp in 1:ncol(sp.plot)) {
			l.nulls=length(which(sp.plot[,sp]==0))
			if((l.spp-l.nulls)<2) {
				drop.plots=cbind(drop.plots, sp)
			}
		}
	
		plot.names.orig=colnames(sp.plot)
		dropped.plots=plot.names.orig[drop.plots]
		if(length(drop.plots)!=0) {
			if(verbose)message({cat("\nThe following plots were dropped from sp.plot:\n\t");cat(dropped.plots, sep=" "); cat("\n")})
			plot.names.orig=names(sp.plot)
			sp.plot=sp.plot[,-as.numeric(drop.plots[!is.na(drop.plots)])]
			if(is.null(ncol(sp.plot))) sp.plot=NULL
		}
	}
# prune all objects with set of matching species
	if(length(missing.spp)!=0) {
		species=intersect(row.names(sp.plot),row.names(sp.plot)[!row.names(sp.plot)%in%missing.spp])
		sp.plot=sp.plot[species,]
	}
	
#if(ncol(sp.plot)==0 || nrow(sp.plot)==0) stop("Datasets cannot be matched together.")

	if(!missing(phy)) {
		if(length(intersect(phy$tip.label, row.names(sp.plot)))!=length(phy$tip.label)) {
			phy=drop.tip(phy, phy$tip.label[!phy$tip.label%in%row.names(sp.plot)])
		}
	}
	
	if(!missing(sp.traits)) {
		if(length(intersect(row.names(sp.traits), row.names(sp.plot)))!=length(row.names(sp.traits))) {
			sp.traits=sp.traits[row.names(sp.plot),]
		}
	}
	
# prepare output
	sp.plot=as.data.frame(sp.plot)
	if(!missing(phy) && missing(sp.traits)) {
		r.out=list(sp.plot, phy)
		names(r.out)=c("sp.plot", "sp.tree")
	} else if(missing(phy) && !missing(sp.traits)) {
		r.out=list(sp.plot, sp.traits)
		names(r.out)=c("sp.plot", "sp.traits")
	} else if(!missing(phy) && !missing(sp.traits)) {
		r.out=list(sp.plot, phy, sp.traits)
		names(r.out)=c("sp.plot", "sp.tree", "sp.traits")
	} else {
		r.out=list(sp.plot)
		names(r.out)=c("sp.plot")
	}	
	return(r.out)
}	