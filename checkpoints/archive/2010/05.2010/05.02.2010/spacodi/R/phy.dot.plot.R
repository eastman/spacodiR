phy.dot.plot<-function(spacodi.object, phy, edge.width=0.2, lab.adj=c(0, 0), tips.adj=c(0.5, 3), tips.cex=1.0, pch=21, print.labs=FALSE, outfile=NULL,...) {
	op=par(no.readonly=TRUE)
	on.exit(par(op))
	if(!is.null(outfile))pdf(outfile)
	
	obj=spacodi.object
	tips=phy$tip.label
	
	obj=as.phylocom(obj, picante=FALSE)
	obj=obj[which(obj$samples!=0),]
	plots=split(obj, obj$plot)
	n=ceiling(sqrt(length(plots)->l.plot))
	layout(matrix(1:n^2,n,n))
	if(n^2>16)warning("Tip labels are unlikely to plot properly with this may trees.")
	for(i in 1:l.plot) {
		tip.labs=vector()
		for(j in 1:length(tips)){
			if(!is.na(match(as.vector(tips[j]), as.vector(plots[[i]]$species)))) tip.labs[j]=1 else tip.labs[j]=0
		}
		plot(phy, no.margin = TRUE, show.tip.label = FALSE, direction="u", edge.width=edge.width, ...)
		tiplabels(pch=pch, bg=ifelse(tip.labs == 1, "black", NA), col=NA, adj=tips.adj, cex=tips.cex)
		legend("bottomleft", legend = names(plots)[i], fill = NULL, bg = NULL, box.lty="blank", border="white", adj=lab.adj)
	}
	if(!is.null(outfile))dev.off()
	if(print.labs) {res=list(plots, phy$tip.label); names(res)=c("groups", "tip.labels"); print(res); dev.new(); plot(phy,direction="u")}
	invisible()
}

