plot.K.through.time <-
function(phy, traits, outfile=NULL){
	if(!is.null(outfile)){
		pdf(outfile)
	}
	if(ncol(traits)>1)stop("Cannot plot more than a single trait at a time")
	as.data.frame(K.all.nodes(phy,traits,return.all=FALSE))->raw
	raw$bt=-raw$node.time
	plot(raw$blomberg.K~raw$bt,ylab="Blomberg's K",xlab="branching.time", main="Blomberg's K through time")
	if(!is.null(outfile)){
		dev.off()
	}
}

