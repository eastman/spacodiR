spacodi.calc <-function(sp.plot, phy = NULL, sp.traits = NA, all.together=TRUE, prune.tree=TRUE, prune.plot=TRUE){
	if(!missing(phy) && !missing(sp.traits)) stop("Please supply either a phylogeny or a trait matrix, but not both.")
	l.spp=nrow(sp.plot)
	drop.plots=vector()
	for(sp in 1:ncol(sp.plot)) {
		l.nulls=length(which(sp.plot[,sp]==0))
		if((l.spp-l.nulls)<2) {
			drop.plots=cbind(drop.plots, sp)
		}
	}
	plot.names.orig=names(sp.plot)
	if(prune.plot==TRUE && length(drop.plots)!=0) {
		plot.names.orig=names(sp.plot)
		if(ncol(sp.plot)<=2)stop("SPACoDi cannot handle this dataset: too few plots available for calculations.")
		sp.plot=sp.plot[,-as.numeric(drop.plots[!is.na(drop.plots)])]
		names(sp.plot)=plot.names.orig[-drop.plots]
	} else if(length(drop.plots)!=0 && prune.plot==FALSE) {
		stop("At least one plot has fewer than two species for which samples are present, causing spacodi.calc() to abort")
	}
	
	sp.plot <-	as.matrix(sp.plot)
	plots   <-	names(sp.plot)
	np      <-	ncol(sp.plot)
	ns      <-	nrow(sp.plot)

# determine type of abundance
	stripped=unname(unlist(c(sp.plot)))
	if(all(!is.na(match(stripped, c(0,1))))) {
		abundt = 0		# presence|absence
	} else if(all(!is.na(match(range(stripped), c(0,1)))) && length(unique(stripped))>2) {
		abundt = 1	# relative abundance
	} else {
		abundt = 2		# abundance (n.individuals)
	}
	
# INTERPRET whether to compute trait diversity or phylogenetic turnover
	if(missing(sp.traits)){ 
	# distmat is phylogenetic or is null: Bst or Ist
		if(!missing(phy) && class(phy)=="phylo") {
			if(any(attr(phy, "names")=="node.label"))phy$node.label=NULL
			l.phy	<-	length(phy$tip.label)
			u.match=union(phy$tip.label,row.names(sp.plot))
			i.match=intersect(phy$tip.label,row.names(sp.plot))
			if(length(u.match)!=length(i.match) && prune.tree==FALSE) stop("Not all species are present for both the tree and the sp.plot.")
			
			phy.match=row.names(sp.plot)%in%phy$tip.label
			# fail if any species in community sample are missing from tree
			if(!all(phy.match)) {
				cat("\n\nThe following are species unsampled in the tree, causing spacodi.calc() to abort:\n\n")
				cat(row.names(sp.plot)[!phy.match])
				stop(cat("\n\n"))
			}
			# if (at least) all community-sampled species are present in tree, prune tree to match
			if(prune.tree==TRUE) phy=prune.tree.sp(phy=phy, spacodi.object=sp.plot, resolve=FALSE)
			distmat <- cophenetic(phy)
		} else {
			distmat <- matrix(1, ncol = ns, nrow = ns)
		}
		diag(distmat) <- 0	
		dim <- nrow(distmat)
		out <- NA
		out <- .C("spacodi", 
				  np = as.integer(np),
				  ns = as.integer(ns),
				  sp.plot = as.double(as.vector(sp.plot)),
				  distmat = as.double(as.vector(as.matrix(distmat))),
				  abundtype = as.integer(abundt),
				  Ndclass = as.integer(0),
				  dclass = as.double(c(0,0)),
				  Ist = 0,
				  Pst = 0,
				  Bst = 0,
				  PIst = 0
				  )
	# compile results for phylogenetic turnover
		if(missing(phy)) {
			r.out=as.numeric(out[8])
			names(r.out)="Ist"
			r.out=as.data.frame(t(r.out))
			return(r.out)
		} else {
			r.out=as.numeric(c(out[8:11]))
			names(r.out)=c("Ist","Pst","Bst","PIst")
			r.out=as.data.frame(t(r.out))
			return(r.out)
		}
		
	} else if(all(row.names(sp.traits)%in%row.names(sp.plot))) { 
	# distmat is trait-based: Tst
		
		if(ncol(sp.traits)==1) all.together=TRUE
		if(all(is.null(names(sp.traits)))) names(sp.traits)=paste("trt",seq(1:ncol(sp.traits)),sep="")
		if(all.together==TRUE){
			distmat=as.matrix(dist(sp.traits))
			dim <- nrow(distmat)
			out <- NA
			out <- .C("spacodi", 
					  np = as.integer(np),
					  ns = as.integer(ns),
					  sp.plot = as.double(as.vector(sp.plot)),
					  distmat = as.double(as.vector(as.matrix(distmat))),
					  abundtype = as.integer(abundt),
					  Ndclass = as.integer(0),
					  dclass = as.double(c(0,0)),
					  Ist = 0,
					  Pst = 0,
					  Bst = 0,
					  PIst = 0
					  )
		# compile results for multivariate trait diversity
			r.out=as.numeric(c(out[8:11]))
			names(r.out)=c("Ist","Tst","T*st","TAUst")
			r.out=as.data.frame(t(r.out))
			return(r.out)
		} else {
			out.array=array(dim=c(ncol(sp.traits), 4))
			for(tt in 1:ncol(sp.traits)) {
				trait=data.frame(sp.traits[,tt])
				row.names(trait)=row.names(sp.traits)
				distmat <- as.matrix(dist(trait))
				dim <- nrow(distmat)
				out <- NA
				out <- .C("spacodi", 
					  np = as.integer(np),
					  ns = as.integer(ns),
					  sp.plot = as.double(as.vector(sp.plot)),
					  distmat = as.double(as.vector(as.matrix(distmat))),
					  abundtype = as.integer(abundt),
					  Ndclass = as.integer(0),
					  dclass = as.double(c(0,0)),
					  Ist = 0,
					  Pst = 0,
					  Bst = 0,
					  PIst = 0
					  )
		# compile results for single sp.traits trait diversity
				out.array[tt,]=as.numeric(c(out[8:11]))
			} 
			r.out=as.data.frame(out.array)
			names(r.out)=c("Ist","Tst","T*st","TAUst")
			row.names(r.out)=names(sp.traits)
			return(r.out)
		}
	} else stop("Supplied datasets cannot be matched to each other.") 	
}

