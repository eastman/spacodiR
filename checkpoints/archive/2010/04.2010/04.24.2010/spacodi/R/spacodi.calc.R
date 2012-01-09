spacodi.calc <-
function(sp_plot, phy = NULL, traits = NA, all.together=TRUE, prune.tree=FALSE, prune.plot=TRUE){
	if(!missing(phy) && !missing(traits)) stop("Please supply either a phylogeny or a trait matrix, but not both.")
	
	l.spp=nrow(sp_plot)
	drop.plots=vector()
	for(sp in 1:ncol(sp_plot)) {
		l.nulls=length(which(sp_plot[,sp]==0))
		if((l.spp-l.nulls)<2) {
			drop.plots[sp]=sp
		}
	}
	
	if(prune.plot==TRUE && length(drop.plots)!=0) {
		sp_plot=sp_plot[,-drop.plots[!is.na(drop.plots)]]
		warning("At least one plot was pruned from the matrix, there being fewer than two species sampled therein.")
	} else if(length(drop.plots)!=0 && prune.plot==FALSE) {
		stop("At least one plot has fewer than two species for which samples are present, causing spacodi.calc() to abort")
	}
	
	sp_plot <-	as.matrix(sp_plot)
	plots   <-	names(sp_plot)
	np      <-	ncol(sp_plot)
	ns      <-	nrow(sp_plot)

# determine type of abundance
	stripped=unname(unlist(c(sp_plot)))
	if(all(!is.na(match(stripped, c(0,1))))) {
		abundt = 0		# presence|absence
	} else if(all(!is.na(match(range(stripped), c(0,1)))) && length(unique(stripped))>2) {
		abundt = 1	# relative abundance
	} else {
		abundt = 2		# abundance (n.individuals)
	}
	
# INTERPRET whether to compute trait diversity or phylogenetic turnover
	if(missing(traits)){ 
	# distmat is phylogenetic or is null: Bst or Ist
		if(!missing(phy) && class(phy)=="phylo") {
			if(!is.binary.tree(phy)) {
				phy=multi2di(phy)
				warning("Supplied tree was not dichotomous and was randomly resolved with ape:::multi2di(). ") 
			}
			l.phy	<-	length(phy$tip.label)
			u.match=union(phy$tip.label,row.names(sp_plot))
			i.match=intersect(phy$tip.label,row.names(sp_plot))
			if(length(u.match)!=length(i.match) && prune.tree==FALSE) stop("Not all species are present for both the tree and the sp_plot.")
			
			phy.match=row.names(sp_plot)%in%phy$tip.label
			# fail if any species in community sample are missing from tree
			if(!all(phy.match)) {
				cat("\n\nThe following are species unsampled in the tree, causing spacodi.calc() to abort:\n\n")
				cat(row.names(sp_plot)[!phy.match])
				stop(cat("\n\n"))
			}
			# if (at least) all community-sampled species are present in tree, prune tree to match
			if(prune.tree==TRUE) {
				tree.drop=phy$tip.label%in%i.match
				if(any(!tree.drop)) phy=drop.tip(phy, phy$tip.label[!tree.drop])
				sp_plot=sp_plot[which(i.match==row.names(sp_plot)),]
				ns=nrow(sp_plot)
				if(length(phy$tip.label)!=l.phy) warning("At least one species was absent from the community samples: tree was pruned to match dataset.")
			}
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
				  sp_plot = as.double(as.vector(sp_plot)),
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
		
	} else if(all(row.names(traits)%in%row.names(sp_plot))) { 
	# distmat is trait-based: Tst
		
		if(ncol(traits)==1) all.together=TRUE
		if(all(is.null(names(traits)))) names(traits)=paste("trt",seq(1:ncol(traits)),sep="")

		if(all.together==TRUE){
			distmat=as.matrix(dist(traits))
			dim <- nrow(distmat)
			out <- NA
			out <- .C("spacodi", 
					  np = as.integer(np),
					  ns = as.integer(ns),
					  sp_plot = as.double(as.vector(sp_plot)),
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
			out.array=array(dim=c(ncol(traits), 4))
			for(tt in 1:ncol(traits)) {
				trait=data.frame(traits[,tt])
				row.names(trait)=row.names(traits)
				distmat <- as.matrix(dist(trait))
				dim <- nrow(distmat)
				out <- NA
				out <- .C("spacodi", 
					  np = as.integer(np),
					  ns = as.integer(ns),
					  sp_plot = as.double(as.vector(sp_plot)),
					  distmat = as.double(as.vector(as.matrix(distmat))),
					  abundtype = as.integer(abundt),
					  Ndclass = as.integer(0),
					  dclass = as.double(c(0,0)),
					  Ist = 0,
					  Pst = 0,
					  Bst = 0,
					  PIst = 0
					  )
		# compile results for single traits trait diversity
				out.array[tt,]=as.numeric(c(out[8:11]))
			} 
			r.out=as.data.frame(out.array)
			names(r.out)=c("Ist","Tst","T*st","TAUst")
			row.names(r.out)=names(traits)
			return(r.out)
		}
	} else stop("Supplied datasets cannot be matched to each other.") 	
}

