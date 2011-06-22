Bst.permutation.plot<-function(Bst.permutations, cex=list(1.6, 0.1, 0.6), col=list("black", "lightgray"), bg=list("white", "lightgray", "black"), all.points=TRUE, add.id=TRUE, sig.plot=TRUE, cut.off=0.05, envelope=FALSE, outfile=NULL, ...) {
	op=par(no.readonly=TRUE)
	on.exit(par(op))
		
	if(is.atomic(cex))cex=list(cex, cex/10, cex/2)
	if(length(cex)==2 && (all.points==FALSE && add.id==TRUE))cex=list(cex[[1]],cex[[1]],cex[[2]])
	
# data gathering
	o.orig=Bst.permutations$observed.Bst[which(!is.na(Bst.permutations$observed.Bst$Bst)),]
	o=as.vector(o.orig$Bst)
	e=Bst.permutations$expected.Bst

# find and define plot limits
	lim.set<-function(obj, scl){
		out=lapply(obj, function(x){foo=x+(x*scl); return(foo)})
		return(unlist(out))
	}
	x.lim=lim.set(c(min(-o.orig$node.time), max(-o.orig$node.time)), 0.05)
	if(all.points==TRUE) {
		range.array=array(dim=c(nrow(e),2))
		for(rr in 1:nrow(e)) {
			r.foo=e[rr, which(!is.na(e[rr,]))]
			if(length(r.foo)>0) {
				range.array[rr,1]=min(r.foo)
				range.array[rr,2]=max(r.foo)
			} else {
				lapply(range.array[rr,], function(obj)obj=NA)
			}
		}
		
		Bst.lim.y=c(l.b<-min(as.vector(o),as.vector(unlist(e)),na.rm=TRUE), u.b<-max(as.vector(o),as.vector(unlist(e)),na.rm=TRUE))
		y.scl=0.35
		y.lim.exp=lim.set(Bst.lim.y, y.scl)
		y.lim=lim.set(c(min(o, na.rm=TRUE), max(o, na.rm=TRUE)), y.scl)
		if(diff(y.lim.exp)>diff(y.lim))warning("Some expected Bsts were excluded from the plot.")

	} else {
		y.lim=lim.set(c(min(o, na.rm=TRUE), max(o, na.rm=TRUE)), y.scl)
	}
# begin plotting
	if(!is.null(outfile))pdf(file=outfile)
	plot(-o.orig$node.time, o.orig$Bst, main="Bst permutations through time", xlab="branching times", ylab="Bst", ylim=y.lim, xlim=c(x.lim[1],0), type="n")
	if(all.points==TRUE) {
		e=as.data.frame(e)
		for(pp in 1:nrow(e)) {
			b.time=o.orig$node.time[which(o.orig$node.ID==row.names(e)[pp])]
			points(rep(-b.time, length(e[pp,])), e[pp,], cex=cex[[2]], col=col[[2]], bg=bg[[2]], pch=21)
		}
	}
	if(envelope==TRUE) {
		poor=which(is.na(range.array[,1]))
		if(length(poor)>0) min.array=range.array[-poor,1] else min.array=range.array[,1]
		if(length(poor)>0) max.array=range.array[-poor,2] else max.array=range.array[,2]
		if(length(poor)>0) o.nodes=-o.orig$node.time[-poor] else o.nodes=-o.orig$node.time
		
		lines(smooth.spline(o.nodes, min.array, ...),lty=2)
		lines(smooth.spline(o.nodes, max.array, ...),lty=2)
	}
	if(sig.plot==TRUE) {
		if(any(names(Bst.permutations)=="randomization.test")) {
			rr=Bst.permutations$randomization.test
		} else {stop("Must have randomization.test results supplied with Bst.permutations.")}
		bg.sig=vector()
		if(length(bg)==3) bg.cols=bg else bg.cols=list("white", "black", "gray")
		for(node in 1:nrow(rr)) {
			if(rr[node, "p.value"] > cut.off || is.na(rr[node, "p.value"])) {
				bg.sig[node]=bg.cols[[1]]
			} else if(rr[node, "p.value"] <= cut.off && ((rr[node, "obs.Bst"] - rr[node, "m.exp.Bst"]) > 0)) {
				bg.sig[node]=bg.cols[[2]]
			} else {
				bg.sig[node]=bg.cols[[3]]
			} 
		}
		points(-o.orig$node.time, o.orig$Bst, cex=cex[[1]], col=col[[1]], pch=21, bg=bg.sig, ...)
		legend("topleft", legend=c("Bst: as expected", "Bst: larger than expected", "Bst: smaller than expected"), fill=c(bg.cols[[1]], bg.cols[[2]], bg.cols[[3]]), bg = NULL, box.lty="blank", border="gray", cex=0.5)
	} else if(sig.plot==FALSE) {
		points(-o.orig$node.time, o.orig$Bst, cex=cex[[1]], col=col[[1]], pch=21, bg=bg[[1]], ...)
	}
	if(add.id==TRUE) {
		require(calibrate)
		textxy(-o.orig$node.time, o.orig$Bst, o.orig$node.ID, m=c(mean(-o.orig$node.time), mean(o.orig$Bst)), cx=cex[[3]])
	}	
	if(!is.null(outfile))dev.off()
	invisible()
}
