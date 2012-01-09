K.by.nodes<-function(sp.traits, phy, obs.only=FALSE, return.all=TRUE, n.rep=10, rand.test=TRUE, r.rep=10000) {
	
# parameter settings
	if(is.null(sp.traits) && is.null(phy)) {
		stop("Must supply both sp.traits and tree.")
	}
	if(is.null(r.rep) && rand.test)stop("Must specify 'r.rep' if using the randomization test.")
	
	if(obs.only==TRUE){
		rand.test=FALSE 
		exp.K=FALSE
	} else {exp.K=TRUE}
	
# internal function calling Kcalc()	
	K.all.nodes <-function(sp.traits, phy, return.all=TRUE){
		if(ncol(sp.traits)!=1)stop("Cannot handle more than a single trait at a time.")
		if (length(phy$tip.label)!=nrow(sp.traits)) stop("Tree cannot be matched to data.")
		tt=as.data.frame(as.numeric(as.vector(sp.traits[,1])))
		row.names(tt)=row.names(sp.traits)
		names(tt)="trait"
		sp.traits=tt
		sub=subtrees(phy)
		if(is.null(rownames(sp.traits))||!all(row.names(sp.traits)%in%phy$tip.label)) stop("Trait matrix must have row names corresponding to tip labels of the phylogeny.")
		n.nodes <- phy$Nnode
		out <- array(dim=c(length(sub), 4)) ### prep an output file
		for (i in 1:n.nodes) {### cycle over subtrees
			trt.i <- data.frame(sp.traits[match(sub[[i]]$tip.label, rownames(sp.traits)),]) #extract the relevant individual from trait matrix, put them into tiplabels order. 
			nn<-sub[[i]]$Nnode
			if(length(trt.i) > 0 & var(trt.i) != 0 & nn > 1){ # eliminate basal polytomies and cases with no trait variance
				K.i <- Kcalc(trt.i, sub[[i]], FALSE)
			} else {
				K.i <- NA	
			}
			out[i, 1] <- K.i
			out[i, 2] <- nn+1
			out[i, 3] <- max(branching.times(sub[[i]]))
			out[i, 4] <- min(sub[[i]]$node)
		}
		array.cols=c("Blomberg.K","tips","node.time","node.ID")
		dimnames(out)=list(NULL, array.cols)
		if(!return.all) {
			out=out[which(!is.na(out[,1])),]
		}
		return(out)
	}
	
# internal function preparing permuted datasets for K.all.nodes()
	K.exp.nodes <- function(phy, n.rep = 10, Brownian=TRUE,...) {
		out <- array(dim=c(phy$Nnode, 4, n.rep)) ### prep an output file
	
		if(Brownian==TRUE){
			trt.sims <- as.data.frame(sim.char(phy=phy, model.matrix=as.matrix(1), nsims=n.rep, model="brownian"))
			names(trt.sims)=paste("trt",seq(1:n.rep),sep=".")
		} else {
			trt.sims <- as.data.frame(replicate(n.rep, runif(length(phy$tip.label))))
			row.names(trt.sims)=phy$tip.label
			names(trt.sims)=paste("trt",seq(1:n.rep),sep=".")
		}
		for(k in 1:n.rep) {
			k.sim=as.data.frame(trt.sims[,k])
			row.names(k.sim)=row.names(trt.sims)
			out[,,k] <- foo <- K.all.nodes(phy=phy, sp.traits=k.sim, return.all=TRUE)
		}	
		dimnames(out)=list(NULL, names(as.data.frame(foo)), paste("iter",seq(1:n.rep),sep="."))
		return(out)
	}
	
# get K for nodes							
	o.foo=K.all.nodes(sp.traits=sp.traits, phy=phy, return.all=ifelse(obs.only==TRUE, TRUE, FALSE))
	o=o.foo[,"Blomberg.K"]
	names(o)=o.foo[,"node.ID"]
	o=o[order(o)]
	n.o=names(o)
	
	o.orig=as.data.frame(o.foo)
	row.names(o.orig)=o.orig$node.ID
	
	
	if(exp.K==TRUE) {	
		
# K permutation								
		e.out=array(dim=c(length(o), n.rep))		
		for(bb in 1:n.rep) {
			b.foo=as.data.frame(K.exp.nodes(sp.traits=sp.traits, phy=phy, n.rep=1)[,,1])
			b.match=match(as.numeric(n.o), b.foo$node.ID)
			e.out[,bb]=as.numeric(b.foo[b.match,"Blomberg.K"])
		}
		
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
				obs=o.orig$Blomberg.K[node]
				nn=o.orig$node.ID[node]
				exp=e.out[node,which(!is.na(e.out[node,]))]
				rand.array[node,1]=ifelse(length(exp)>0, randomization.test.sp(obs=obs,exp=exp, iter=r.rep, return.all=FALSE, two.tailed=TRUE), NA)
				rand.array[node,2]=nn
				rand.array[node,3]=obs
				rand.array[node,4]=ifelse(length(exp)>0, mean(exp), NA)
				rand.array[node,5]=length(exp)
			}
			r.test=as.data.frame(rand.array)
			names(r.test)=c("p.value","node.ID","obs.K","m.exp.K","valid.comparisons")
			row.names(r.test)=r.test$node.ID
			if(return.all==TRUE) {
				o.orig=as.data.frame(K.all.nodes(sp.traits=sp.traits, phy=phy, return.all=TRUE))
				row.names(o.orig)=o.orig$node.ID
			} else {o.orig=o.orig}
			K.out=list(o.orig, e.out, r.test)
			names(K.out)=c("observed.K","expected.K","randomization.test")
			return(K.out)
		} else {
			if(return.all==TRUE) {
				o.orig=as.data.frame(K.all.nodes(sp.traits=sp.traits, phy=phy, return.all=TRUE))
				row.names(o.orig)=o.orig$node.ID
			} else {o.orig=o.orig}
			K.out=list(o.orig, e.out)
			names(K.out)=c("observed.K","expected.K")
			return(K.out)
		}
	} else {
		K.out=list(o.orig)
		names(K.out)="observed.K"
		return(K.out)
	}
}


