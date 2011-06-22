spacodi.calc <-function(sp.plot, phy = NULL, sp.traits = NA, all.together=TRUE, prune=TRUE, ...){
	if(!missing(phy) && !missing(sp.traits)) stop("Please supply either a phylogeny or a trait matrix, but not both.")
		
# INTERPRET whether to compute trait diversity or phylogenetic turnover or species turnover
	
	if(missing(sp.traits)){ 
		
		if(!missing(phy) && class(phy)=="phylo") {
		  # distmat is phylogenetic: Bst
			sp.data=match.spacodi.data(sp.plot=sp.plot, phy=phy, prune=prune, ...)
			sp.plot=sp.data$sp.plot
			phy=sp.data$sp.tree
			distmat <- cophenetic(phy)
			
		} else {
		  # distmat is null: Ist
			sp.data=match.spacodi.data(sp.plot=sp.plot, prune=prune, ...)
			sp.plot=sp.data$sp.plot
			distmat <- matrix(1, ncol = nrow(sp.plot), nrow = nrow(sp.plot))
		}
		
	} else if(!missing(sp.traits)) { 
		if(ncol(sp.traits)==1) all.together=TRUE
		if(all(is.null(names(sp.traits)))) names(sp.traits)=paste("trt",seq(1:ncol(sp.traits)),sep="")	
		
		if(all.together==TRUE) {
		  # distmat is trait-based: Tst
			sp.data=match.spacodi.data(sp.plot=sp.plot, sp.traits=sp.traits, prune=prune, ...)
			sp.plot=sp.data$sp.plot
			sp.traits=sp.data$sp.traits
			distmat=as.matrix(dist(sp.traits))
			dim <- nrow(distmat)
		
		} else if(all.together==FALSE) {
		  # distmat is trait-based: Tst for separate traits
			out.array=array(dim=c(ncol(sp.traits), 4))
			for(tt in 1:ncol(sp.traits)) {
				trait.tt=data.frame(sp.traits[,tt])
				row.names(trait.tt)=row.names(sp.traits)
				sp.data=match.spacodi.data(sp.plot=sp.plot, sp.traits=trait.tt, prune=prune, ...)
				distmat <- as.matrix(dist(trait.tt))
				dim <- nrow(distmat)
				  
			  # determine type of abundance
				stripped=unname(unlist(c(sp.plot)))
				if(all(!is.na(match(stripped, c(0,1))))) {
					abundt = 0		# presence|absence
				} else if(all(!is.na(match(range(stripped), c(0,1)))) && length(unique(stripped))>2) {
					abundt = 1	# relative abundance
				} else {
					abundt = 2		# abundance (n.individuals)
				}
				  
				sp.plot <-	as.matrix(sp.plot)
				plots   <-	names(sp.plot)
				ns      <-	nrow(sp.plot)
				np      <-	ncol(sp.plot)
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
				
				out.array[tt,]=as.numeric(c(out[8:11]))
			} 
			r.out=out.array
		} 
	} else {
		stop("Supplied dataset(s) do not appear appropriate for spacodi.calc().") 	
	}		  
			
  # determine type of abundance
	stripped=unname(unlist(c(sp.plot)))
	if(all(!is.na(match(stripped, c(0,1))))) {
		abundt = 0		# presence|absence
	} else if(all(!is.na(match(range(stripped), c(0,1)))) && length(unique(stripped))>2) {
		abundt = 1	# relative abundance
	} else {
		abundt = 2		# abundance (n.individuals)
	}
	
	sp.plot <-	as.matrix(sp.plot)
	plots   <-	names(sp.plot)
	ns      <-	nrow(sp.plot)
	np      <-	ncol(sp.plot)
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
	if(missing(phy) && missing(sp.traits)) {
		r.out=as.numeric(out[8])
		names(r.out)="Ist"
		r.out=as.data.frame(t(r.out))

	} else if(missing(sp.traits)){
		r.out=as.numeric(c(out[8:11]))
		names(r.out)=c("Ist","Pst","Bst","PIst")
		r.out=as.data.frame(t(r.out))

	} else if(!missing(sp.traits) && all.together==TRUE) {
		r.out=as.numeric(c(out[8:11]))
		names(r.out)=c("Ist","Tst","T*st","TAUst")
		r.out=as.data.frame(t(r.out))
	} else {
		r.out=as.data.frame(r.out)
		names(r.out)=c("Ist","Tst","T*st","TAUst")
		row.names(r.out)=names(sp.traits)
	}
	return(r.out)
}
		
