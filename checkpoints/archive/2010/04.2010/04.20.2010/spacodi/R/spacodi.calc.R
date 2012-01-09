spacodi.calc <-
function(sp_plot, phy = NULL, sptraits = NA, which.plots = NA){
# dyn.load("./spacodi.Wrk/src/spacodi.so")
# preceding line should not be necessary with zzz.R; uncomment and re-install if issues
	if(!is.na(which.plots[1])) sp_plot <- sp_plot[,which(names(sp_plot) %in% which.plots)] # eliminate cetain plots, if desired
	sp_plot <- as.matrix(sp_plot)
	plots   <- names(sp_plot)
	np      <- ncol(sp_plot)
	ns      <- nrow(sp_plot)  # how many spp
	if(!missing(phy) & !missing(sptraits)) stop("Please supply either a phylogeny or a trait matrix, but not both.")
	
	# generate the distance matrix to be used
	if(!missing(phy)){ # distmat is phylogenetic: Bst
		distmat <- cophenetic(phy)
	} else {
		if(!all(is.na(sptraits))) { #distmat is trait-based: Tst
			distmat <- as.matrix(dist(sptraits))
		} else { #distmat is species ID: Ist
			distmat <- matrix(1, ncol = ns, nrow = ns)
			diag(distmat) <- 0	
		}
	}
	dim <- nrow(distmat)
	out <- NA
	out <- .C("spacodi", 
		np = as.integer(np),
		ns = as.integer(ns),
		sp_plot = as.double(as.vector(sp_plot)),
		distmat = as.double(as.vector(as.matrix(distmat))),
		abundtype = as.integer(2),
		Ndclass = as.integer(0),
		dclass = as.double(c(0,0)),
		Ist = 0,
		Pst = 0,
		Bst = 0,
		PIst = 0
		)
	if(!all(is.na(sptraits))) names(out)[9:11] <- c("Tst", "T*st", "TAUst")
	if(missing(phy) & missing(sptraits)){
		return(out[8])
		} else {
		return(out[8:11])
		} 	
}

