spacodi.matrices<-function(sp_plot, phy) {
	orig=names(sp_plot)
	l.spp=nrow(sp_plot)
	drop.plots=vector()
	for(sp in 1:ncol(sp_plot)) {
		l.nulls=length(which(sp_plot[,sp]==0))
		if((l.spp-l.nulls)<2) {
			drop.plots[sp]=sp
		}
	}
	
	if(length(drop.plots)!=0) {
		sp_plot=sp_plot[,-drop.plots[!is.na(drop.plots)]]
		warning("At least one plot was pruned from the matrix, there being fewer than two species sampled therein.")
	} 
	
	plots=names(sp_plot)
	l.plots=length(plots)
	sp.array=array(dim=c(l.plots, l.plots, 4))
	out.array=array(dim=c(l.plots, 1, 4))
	r.out=list()
	for(row in 1:l.plots) {
		rr=plots[row]
		r=which(plots==rr)
		for(cc in plots[which(plots!=rr)]) {
			c=which(plots==cc)
			keep=plots[!is.na(match(plots, c(cc,rr)))]
			sp.rc=sp_plot[,keep]
			ss=try(spacodi.calc(sp_plot=sp.rc, phy=phy, prune.tree=TRUE, prune.plot=FALSE),silent=TRUE)
			if(class(ss)!="try-error") {
				class(ss)="vector"
				s.names=names(ss)
				for(xx in 1:length(ss)) {
					sp.array[r,c,xx]=unlist(unname(ss[which(names(ss)==names(ss)[xx])]))
				}
			} #else {lapply(sp.array[r,c,], function(obj)obj=NA)}
		}
		for(yy in 1:4) {
			foo=sp.array[r,,yy]
			out.array[r,1,yy]=sum(foo[which(!is.na(foo))])/(l.plots-1)
		}
	}
	for(zz in 1:4) {
		foo=as.data.frame(out.array[,,zz])
		row.names(foo)=plots
		d.matrix=dist(foo)
		r.out[[zz]]=d.matrix
	}
	names(r.out)=c("Ist","Pst","Bst","PIst")
	return(r.out)
}