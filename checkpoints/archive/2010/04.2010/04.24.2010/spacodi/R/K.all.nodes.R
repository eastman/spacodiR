K.all.nodes <-
function(phy, traits, return.all=TRUE){
	require(picante)
	if (!is.binary.tree(phy)) {
		phy=multi2di(phy)
		warning("Supplied tree was not dichotomous and was randomly resolved with ape:::multi2di(). ") 
	}
	if (length(phy$tip.label)!=nrow(traits)) stop("Tree cannot be matched to data")
	sub=subtrees(phy)
	n.traits <- ncol(traits)
	if(length(unique(colnames(traits))) < n.traits) {trait.names=paste("trt",seq(1:n.traits),sep=".")} else {trait.names = colnames(traits)}
	if(is.null(rownames(traits))||!all(row.names(traits)%in%phy$tip.label)) stop("Trait matrix must have row names corresponding to tip labels of the phylogeny.")
	n.nodes <- phy$Nnode
	out <- array(dim=c(length(sub), 4, n.traits)) ### prep an output file
	array.cols<-c("blomberg.K","tips","node.time","node.ID")
	for (i in 1:n.nodes) {### cycle over subtrees
		trt.i <- data.frame(traits[match(sub[[i]]$tip.label, rownames(traits)),]) #extract the relevant individual from trait matrix, put them into tiplabels order. 
		for (j in 1:n.traits){ 
			trt.ij <- trt.i[,j]
			nn<-sub[[i]]$Nnode
			if(length(trt.ij) > 0 & var(trt.ij) != 0 & nn > 1){ # eliminate basal polytomies and cases with no trait variance
				K.ij <- Kcalc(trt.ij, sub[[i]], FALSE)
				} else {
				K.ij <- NA	
			}
			out[i, 1, j] <- K.ij
			out[i, 2, j] <- nn+1
			out[i, 3, j] <- max(branching.times(sub[[i]]))
			out[i, 4, j] <- min(sub[[i]]$node)
		}
		
	}
	dimnames(out)=list(NULL, array.cols, trait.names)
	if(!return.all) {
		out=out[which(!is.na(out[,1,1])),,]
		cat("Clades with two taxa were excluded as Blomberg's K is incalculable for small clades.\n\tIf this is unintended behavior, switch the 'return.all' operator to 'true'.\n\n")
	}
	return(out)
}

