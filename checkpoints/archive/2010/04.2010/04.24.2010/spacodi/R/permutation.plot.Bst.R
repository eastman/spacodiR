permutation.plot.Bst<-function(phy, sp_plot, n.rep=10, method="1s", all.points=TRUE, envelope=TRUE, add.id=FALSE, cex=list(2, 0.2, 1), outfile=NULL, ...) {
	op=par(no.readonly=TRUE)
	on.exit(par(op))
	
	if(length(cex)==1)cex=list(cex, cex/10, cex/2)
	
	if(is.null(sp_plot) && is.null(phy)) {
		stop("Must supply both an sp_plot and tree.")
	}
	
	# get Bsts for nodes						# for empirical data plot
	o.foo=Bst.all.nodes(sp_plot=sp_plot, phy=phy, return.all=FALSE)
	o.orig=as.data.frame(o.foo)
	o=o.foo[,"Bst"]
	names(o)=o.foo[,"node.ID"]
	o=o[order(o)]
	n.o=names(o)
		
	# BST PLOT: randomization
	r.out=array(dim=c(length(o), n.rep))		# for all randomization
	for(bb in 1:n.rep) {
		b.foo=as.data.frame(Bst.permutation(sp_plot=sp_plot, phy=phy, n.rep=1, method=method)[,,1])
		b.match=match(as.numeric(n.o), b.foo$node.ID)
		r.out[,bb]=as.numeric(b.foo[b.match,"Bst"])
	}
	
	range.array=array(dim=c(length(o), 2))		# for spline
	for(rr in 1:length(o)) {
		r.foo=r.out[rr, which(!is.na(r.out[rr,]))]
		if(length(r.foo)>0) {
			range.array[rr,1]=min(r.foo)
			range.array[rr,2]=max(r.foo)
		} else {
			lapply(range.array[bb,], function(obj)obj=NA)
		}
	}
	Bst.lim.y=c(l.b<-min(as.vector(o),as.vector(r.out),na.rm=TRUE), u.b<-max(as.vector(o),as.vector(r.out),na.rm=TRUE))
	Bst.lim.x=c(min(-o.orig$node.time), max(-o.orig$node.time))
	lim.set<-function(obj, scl){obj[1]=obj[1]-(scl*obj[1]); obj[2]=obj[2]+(scl*obj[2]); return(obj)}
	
	y.lim=lim.set(Bst.lim.y, 0.25)
	x.lim=lim.set(Bst.lim.x, 0.20)
	
	if(!is.null(outfile))pdf(file=outfile)
	plot(-o.orig$node.time, o.orig$Bst, main="Bst permutations through time", xlab="branching times", ylab="Bst", ylim=y.lim, xlim=c(x.lim[1],0), cex=cex[[1]],...)
	if(all.points==TRUE) {
		r.out=as.data.frame(r.out)
		row.names(r.out)=as.numeric(n.o)
		for(pp in 1:nrow(r.out)) {
			b.time=o.orig$node.time[which(o.orig$node.ID==row.names(r.out)[pp])]
			points(rep(-b.time, length(r.out[pp,])), r.out[pp,], cex=cex[[2]],...)
		}
	}
	if(envelope==TRUE) {
		lines(smooth.spline(-o.orig$node.time, range.array[,1], df=3),lty=2)
		lines(smooth.spline(-o.orig$node.time, range.array[,2], df=3),lty=2)
	}
	if(add.id==TRUE) {
		require(calibrate)
		textxy(-o.orig$node.time, o.orig$Bst, o.orig$node.ID, m=c(mean(-o.orig$node.time), mean(o.orig$Bst)), cx=cex[[3]])
	}
	if(!is.null(outfile))dev.off()
	invisible()
}
