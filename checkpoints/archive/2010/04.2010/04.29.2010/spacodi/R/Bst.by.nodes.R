Bst.by.nodes<-function(sp.plot, phy, obs.only=FALSE, return.all=TRUE, n.rep=10, method="1s", parm=NULL, dmat=NULL, rand.test=TRUE, r.rep=10000) {

# parameter settings
	if(is.null(sp.plot) && is.null(phy)) {
		stop("Must supply both an sp.plot and tree.")
	}
	if(is.null(r.rep) && rand.test)stop("Must specify 'r.rep' if using the randomization test.")

	if(obs.only==TRUE){
		rand.test=FALSE 
		exp.Bst=FALSE
	} else {exp.Bst=TRUE}
	
# internal function calling spacodi.calc()	
	Bst.all.nodes <- function(sp.plot, phy, return.all=TRUE){
		if(any(attr(phy, "names")=="node.label"))phy$node.label=NULL

		if(length(intersect(row.names(sp.plot), phy$tip.label)->i.match)!=nrow(sp.plot) || length(i.match)!=length(phy$tip.label))stop("SPACoDi cannot proceed: dataset cannot be matched.")
		sub=subtrees(phy)
		n.nodes <- phy$Nnode
		out <- array(dim=c(length(sub), 4)) ### prep an output file
		for (i in 1:n.nodes) {
			sp.plot.i <- sp.plot[rownames(sp.plot) %in% sub[[i]]$tip.label,]
			if(nrow(sp.plot.i) > 0){
				bb <- try(spacodi.calc(sp.plot = sp.plot.i, phy = sub[[i]], prune.tree=FALSE, prune.plot=TRUE), silent=TRUE)
				if(class(bb)!="try-error") Bst.i=bb$Bst else Bst.i=NaN
			} else {
				Bst.i <- NA	# should never see 'NA' in the output array
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

# internal function preparing permuted datasets for Bst.all.nodes()
	Bst.exp.nodes <- function(sp.plot, phy, n.rep = 10, method="1s", parm=NULL, dmat=NULL,...) {
		if(any(attr(phy, "names")=="node.label"))phy$node.label=NULL

		sub=subtrees(phy)
		n.plot=ncol(sp.plot)
		n.taxa=Ntip(phy)
		out <- array(dim=c(length(sub), 4, n.rep)) ### prep an output file
		if(method=="r.sp.plot"){
			for(j in 1:n.rep){
				foo <- Bst.all.nodes(phy = phy, sp.plot = r.sp.plot(sp.plot, n.rep=1)[[1]])
				out[,,j] <- foo
			}
		} else if(method=="1s"){
			for(j in 1:n.rep){
				out[,,j] <- foo <- Bst.all.nodes(phy = phy, sp.plot = resamp.1s(sp.plot))
			}
		} else if(method=="1a" && !is.null(parm$abund.class.ratio)){
			for(j in 1:n.rep){
				out[,,j] <- foo <- Bst.all.nodes(phy = phy, sp.plot = resamp.1a(sp.plot, parm$abund.class.ratio))
			}
		} else if(method=="2s"){
			for(j in 1:n.rep){
				out[,,j] <- foo <- Bst.all.nodes(phy = phy, sp.plot = resamp.2s(sp.plot))
			}
		} else if(method=="2x" && !is.null(parm$level)){
			for(j in 1:n.rep){
				out[,,j] <- foo <- Bst.all.nodes(phy = phy, sp.plot = resamp.2x(sp.plot, parm$level))
			}
		} else if(method=="3i"){
			for(j in 1:n.rep){
				out[,,j] <- foo <- Bst.all.nodes(phy = phy, sp.plot = resamp.3i(sp.plot))
			}
		} else if(method=="3t"){
			if(!is.null(dmat)) dmat=dmat else dmat=NULL
			for(j in 1:n.rep){
				out[,,j] <- foo <- Bst.all.nodes(phy = phy, sp.plot = resamp.3t(sp.plot, dmat))
			}
		} else if(method=="3x" && !is.null(parm$level)){
			for(j in 1:n.rep){
				out[,,j] <- foo <- Bst.all.nodes(phy = phy, sp.plot = resamp.3x(sp.plot, parm$level))
			}
		} else stop(cat("Unrecognized method declaration or insufficient information supplied:\n\tIs 'parm' set to a non-null value?\n\n"))
		dimnames(out)=list(NULL, names(as.data.frame(foo)), paste("iter",seq(1:n.rep),sep="."))
		return(out)
	}
	
# get Bsts for nodes							
	o.foo=Bst.all.nodes(sp.plot=sp.plot, phy=phy, return.all=ifelse(obs.only==TRUE, TRUE, FALSE))
	o=o.foo[,"Bst"]
	names(o)=o.foo[,"node.ID"]
	o=o[order(o)]
	n.o=names(o)
	
	o.orig=as.data.frame(o.foo)
	row.names(o.orig)=o.orig$node.ID

	
	if(exp.Bst==TRUE) {	
# BST permutation
		ticker=c(1:35)*ceiling(n.rep/35)
		cat("\nSPACoDi calculations in progress:\n")
			                                           							
		e.out=array(dim=c(length(o), n.rep))		
		for(bb in 1:n.rep) {
			b.foo=as.data.frame(Bst.exp.nodes(sp.plot=sp.plot, phy=phy, n.rep=1, method=method, parm=parm)[,,1])
			b.match=match(as.numeric(n.o), b.foo$node.ID)
			e.out[,bb]=as.numeric(b.foo[b.match,"Bst"])
			
			if(bb%in%ticker)cat(".")
		}
		cat("\n")
# dataset matching
		e.out=as.data.frame(e.out)
		row.names(e.out)=as.numeric(n.o)
		names(e.out)=paste("iter",seq(1:n.rep),sep="")
		match(as.numeric(o.orig$node.ID),as.numeric(row.names(e.out)))->orig.match
		e.out=e.out[orig.match,]
	
# randomization test
		if(rand.test==TRUE) {
			rand.array=array(dim=c(nrow(o.orig),5))
			for(node in 1:nrow(o.orig)) {
				obs=o.orig$Bst[node]
				nn=o.orig$node.ID[node]
				exp=e.out[node,which(!is.na(e.out[node,]))]
				rand.array[node,1]=ifelse(length(exp)!=0, randomization.test.sp(obs=obs,exp=exp, iter=r.rep, return.all=FALSE, two.tailed=TRUE), NA)
				rand.array[node,2]=nn
				rand.array[node,3]=obs
				rand.array[node,4]=ifelse(length(exp)!=0, mean(exp), NA)
				rand.array[node,5]=length(exp)
			}
			r.test=as.data.frame(rand.array)
			names(r.test)=c("p.value","node.ID","obs.Bst","m.exp.Bst","valid.comparisons")
			row.names(r.test)=r.test$node.ID
			if(return.all==TRUE) {
				o.orig=as.data.frame(Bst.all.nodes(sp.plot=sp.plot, phy=phy, return.all=TRUE))
				row.names(o.orig)=o.orig$node.ID
			} else {o.orig=o.orig}
			Bst.out=list(o.orig, e.out, r.test)
			names(Bst.out)=c("observed.Bst","expected.Bst","randomization.test")
			return(Bst.out)
		} else {
			if(return.all==TRUE) {
				o.orig=as.data.frame(Bst.all.nodes(sp.plot=sp.plot, phy=phy, return.all=TRUE))
				row.names(o.orig)=o.orig$node.ID
			} else {o.orig=o.orig}
			Bst.out=list(o.orig, e.out)
			names(Bst.out)=c("observed.Bst","expected.Bst")
			return(Bst.out)
		}
	} else {
		Bst.out=list(o.orig)
		names(Bst.out)="observed.Bst"
		return(Bst.out)
	}
}