Bst.all.nodes <-
function(sp_plot, phy, return.all=TRUE){
	require(ape)
	if (!is.binary.tree(phy)) {
		phy=multi2di(phy)
		warning("Supplied tree was not dichotomous and was randomly resolved with ape:::multi2di(). ") 
	}
	sub=subtrees(phy)
	sp_plot <- as.matrix(sp_plot)
	n.plots <- ncol(sp_plot)
	if(is.null(rownames(sp_plot))) stop("Species-plot matrix must have row names corresponding to tip labels of the phylogeny.")
	n.nodes <- phy$Nnode
	out <- array(dim=c(length(sub), 4)) ### prep an output file
	for (i in 1:n.nodes) {
		sp_plot.i <- sp_plot[rownames(sp_plot) %in% sub[[i]]$tip.label,]
		if(nrow(sp_plot.i) > 0){
			bb <- try(spacodi.calc(sp_plot.i, phy = sub[[i]])$Bst, TRUE)
			Bst.i <- ifelse(class(bb)!="try-error", bb, NA)
		} else {
			Bst.i <- NA	
		}
		out[i, 1] <- Bst.i
		out[i, 2] <- Ntip(sub[[i]])
		out[i, 3] <- max(branching.times(sub[[i]]))
		out[i, 4] <- min(sub[[i]]$node)
	}
	array.cols<-c("Bst","tips","node.time","node.ID")
	dimnames(out)=list(NULL, array.cols)
	if(!return.all) {
		out=out[which(!is.na(out[,1])),]
	}
	return(out)
}

