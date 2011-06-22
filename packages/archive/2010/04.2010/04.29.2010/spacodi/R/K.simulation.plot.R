K.simulation.plot<-function(K.simulations, cex=list(1.6, 0.1, 0.6), col=list("black", "lightgray"), bg=list("white", "lightgray", "black"), all.points=TRUE, add.id=TRUE, sig.plot=TRUE, cut.off=0.05, envelope=FALSE, outfile=NULL, ...) {
	op=par(no.readonly=TRUE)
	on.exit(par(op))
		
	if(is.atomic(cex))cex=list(cex, cex/10, cex/2)
	if(length(cex)==2 && (all.points==FALSE && add.id==TRUE))cex=list(cex[[1]],cex[[1]],cex[[2]])
	
# data gathering
	o.orig=K.simulations$observed.K
	o.orig=K.simulations$observed.K[which(!is.na(K.simulations$observed.K$Blomberg.K)),]

	o=as.vector(o.orig$Blomberg.K)
	e=K.simulations$expected.K

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
		
		K.lim.y=c(l.b<-min(as.vector(o),as.vector(unlist(e)),na.rm=TRUE), u.b<-max(as.vector(o),as.vector(unlist(e)),na.rm=TRUE))
		y.scl=0.35
		y.lim.exp=lim.set(K.lim.y, y.scl)
		y.lim=lim.set(c(min(o, na.rm=TRUE), max(o, na.rm=TRUE)), y.scl)
#	if(diff(y.lim.exp)>diff(y.lim))message("Some expected Ks were excluded from the plot (but were used in the credible-envelope calculation).")

	} else {
		y.lim=lim.set(c(min(o, na.rm=TRUE), max(o, na.rm=TRUE)), y.scl)
	}
# begin plotting
	if(!is.null(outfile))pdf(file=outfile)
	plot(-o.orig$node.time, o.orig$Blomberg.K, main="K simulations through time", xlab="branching times", ylab="Blomberg's K", ylim=K.lim.y, xlim=c(x.lim[1],0), type="n")
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
		
		lines(smooth.spline(o.nodes, min.array),lty=2)
		lines(smooth.spline(o.nodes, max.array),lty=2)
	}
	if(sig.plot==TRUE) {
		if(any(names(K.simulations)=="randomization.test")) {
			rr=K.simulations$randomization.test
		} else {stop("Must have randomization.test results supplied with K.simulations.")}
		bg.sig=vector()
		if(length(bg)==3) bg.cols=bg else bg.cols=list("white", "black", "gray")
		for(node in 1:nrow(rr)) {
			if(rr[node, "p.value"] > cut.off || is.na(rr[node, "p.value"])) {
				bg.sig[node]=bg.cols[[1]]
			} else if(rr[node, "p.value"] <= cut.off && ((rr[node, "obs.K"] - rr[node, "m.exp.K"]) > 0)) {
				bg.sig[node]=bg.cols[[2]]
			} else {
				bg.sig[node]=bg.cols[[3]]
			} 
		}
		points(-o.orig$node.time, o.orig$Blomberg.K, cex=cex[[1]], col=col[[1]], pch=21, bg=bg.sig, ...)
		legend("topleft", legend=c("K: as expected", "K: larger than expected", "K: smaller than expected"), fill=c(bg.cols[[1]], bg.cols[[2]], bg.cols[[3]]), bg = NULL, box.lty="blank", border="gray", cex=0.5)
	} else if(sig.plot==FALSE) {
		points(-o.orig$node.time, o.orig$Blomberg.K, cex=cex[[1]], col=col[[1]], pch=21, bg=bg[[1]], ...)
	}
	if(add.id==TRUE) {
		require(calibrate)
		textxy(-o.orig$node.time, o.orig$Blomberg.K, o.orig$node.ID, m=c(mean(-o.orig$node.time), mean(o.orig$Blomberg.K)), cx=cex[[3]])
	}
	if(!is.null(outfile))dev.off()
	invisible()
}
