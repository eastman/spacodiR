plot.Bst.through.time <-
function(phy, sp_plot, outfile=NULL){
	if(!is.null(outfile)){
		pdf(outfile)
	}
	as.data.frame(Bst.all.nodes(phy,sp_plot))->raw
	raw$bt=-raw$node.time
	plot(raw$Bst~raw$bt,ylab="Bst",xlab="branching.time", main="Bst through time")
	if(!is.null(outfile)){
		dev.off()
	}
}

