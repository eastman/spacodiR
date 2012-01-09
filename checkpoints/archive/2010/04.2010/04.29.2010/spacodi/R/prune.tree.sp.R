prune.tree.sp<-function(phy, spacodi.object, resolve=FALSE) {

	if(is.null(row.names(spacodi.object)))stop("SPACoDi.prune.tree cannot proceed: supplied object must have species as row.names")
	if(class(phy)=="phylo") {
		if(any(attr(phy, "names")=="node.label"))phy$node.label=NULL
	} else {stop("Tree does not appear to be a 'phylo' object")}
	
	l.phy	<-	length(phy$tip.label)
	u.match=union(phy$tip.label,row.names(spacodi.object))
	i.match=intersect(phy$tip.label,row.names(spacodi.object))
	tree.drop=phy$tip.label%in%i.match
	if(any(!tree.drop)) phy=drop.tip(phy, phy$tip.label[!tree.drop])
	
	n.phy <- length(phy$tip.label)
	
	if(n.phy!=l.phy) message("At least one species was pruned to match dataset.")
	if(n.phy<nrow(spacodi.object)) warning("At least one species from the spacodi.object cannot be found in the tree.")
	if(!is.binary.tree(phy) && resolve==TRUE) {
		phy=multi2di(phy)
		warning("Supplied tree was not dichotomous and was randomly resolved with ape:::multi2di(). ") 
	}
	return(phy)
}