# in a shell, run  
# R CMD SHLIB spacodi.c
# to make the shared library that gets called from within R

###########################################################################################################
# Ist: see Hardy and Senterre (2007) - a measure of local species-identity excess between individuals expressing species turnover. 
#		It is a form of spatial partition of Gini-Simpson diversity (equivalent to Fst in population genetics). If large, many species are unshared between sites.
#	I-statistics (given when abundances rather than just presence/absence data are provided in the species-plot matrix) are based on probabilities of species identity for pairs of individuals (they partition Gini-Simpson diversity into alpha and beta components; Hardy & Senterre 2007).
#		Diw is the probability that two individuals from a plot belong to distinct species 
#		Dia is the probability that two individuals from different plots belong to distinct species
#		Ist = (Dia-Diw)/Dia and thus is an analogue of Fst in population genetics
#
# Pst: see Hardy and Senterre (2007)
#	P-statistics (given when a phylogeny and abundances data are provided in the species-plot matrix) are based on mean phylogenetic distances between individuals (Hardy & Senterre 2007).
#		Dpw is the mean phylogenetic distance between individuals from the same plot
#		Dpa is the mean phylogenetic distance between individuals from different plots
#		Pst = (Dpa-Dpw)/Dpa  is an analogue of Nst in population genetics
#
# Bst: ~P*st: see Hardy and Jost (2008) and Hardy and Senterre (2007) - proportional to Pst - Ist
#		relative increase of mean phylogenetic distance between INDIVIDUALS from different species, sampled among versus within sites
#	P*-statistics (given when a phylogeny and abundances data are provided in the species-plot matrix) are based on mean phylogenetic distances between individuals of distinct species (Hardy & Jost 2008).
#		Dp*w is the mean phylogenetic distance between individuals of distinct species from the same plot
#		Dp*a is the mean phylogenetic distance between individuals of distinct species from different plots
#		P*st = (Dp*a - Dp*w) / Dp*a  
#
# PIst: see Hardy and Senterre (2007) - relative increase of mean phylogenetic distance between SPECIES, sampled among versus within sites.
#		A Pst analogue for presence/absence data (written “PIst” in the result file), expressing phylogenetic turnover (independently of species turnover).
#	PI-statistics (given when a phylogeny and a species-plot matrix are provided) are based on mean phylogenetic distances between species, not accounting for species abundances (Hardy & Senterre 2007).
#		DELTAw is the mean phylogenetic distance between species from the same plot
#		DELTAa is the mean phylogenetic distance between species from different plots
#		PIst = (DELTAa-DELTAw)/DELTAa  
#
# INTERPRETATION of RESULTS
#	Spatial phylogenetic CLUSTERING 
#		species within plots are MORE related on average than species from distinct plots 
#		supported where Pst > Ist;  P*st > 0; or  PIst > 0
#	Spatial phylogenetic overdispersion 
#		species within plots are LESS related on average than species from distinct plots 
#		occurs when Pst < Ist;  P*st < 0; or  PIst < 0
#
# 
###########################################################################################################

### MAIN FUNCTIONS ###

# calculates Hardy Bst, Ist, Pst, and analogues
spacodi.calc <- function(sp_plot, phy = NULL, sptraits = NA, which.plots = NA){
	dyn.load("./spacodi.Wrk/src/spacodi.so")
	if(!is.na(which.plots[1])) sp_plot <- sp_plot[,which(names(sp_plot) %in% which.plots)] # eliminate cetain plots, if desired
	sp_plot <- as.matrix(sp_plot)
	plots   <- names(sp_plot)
	np      <- ncol(sp_plot)
	ns      <- nrow(sp_plot)  # how many spp
	if(!missing(phy) & !missing(sptraits)) stop("Please supply either a phylogeny or a trait matrix, but not both.")
	
	# generate the distance matrix to be used
	if(!missing(phy)){ # distmat is phylogenetic: Bst
		distmat <- cophenetic(phy)
	} else {
		if(!all(is.na(sptraits))) { #distmat is trait-based: Tst
			distmat <- as.matrix(dist(sptraits))
		} else { #distmat is species ID: Ist
			distmat <- matrix(1, ncol = ns, nrow = ns)
			diag(distmat) <- 0	
		}
	}
	dim <- nrow(distmat)
	out <- NA
	out <- .C("spacodi", 
		np = as.integer(np),
		ns = as.integer(ns),
		sp_plot = as.double(as.vector(sp_plot)),
		distmat = as.double(as.vector(as.matrix(distmat))),
		abundtype = as.integer(2),
		Ndclass = as.integer(0),
		dclass = as.double(c(0,0)),
		Ist = 0,
		Pst = 0,
		Bst = 0,
		PIst = 0
		)
	if(!all(is.na(sptraits))) names(out)[9:11] <- c("Tst", "T*st", "TAUst")
	if(missing(phy) & missing(sptraits)){
		return(out[8])
		} else {
		return(out[8:11])
		} 	
}


# Calculate Hardy Bst for all subnodes of a phylogeny
Bst.obs.nodes <- Bst.all.nodes <- function(phy, sp_plot, return.all=TRUE){
	require(ape)
	if (!is.binary.tree(phy)) {
		phy=multi2di(phy)
		warning("Supplied tree was not dichotomous and was randomly resolved with ape:::multi2di(). ") 
	}
	sub=subtrees(phy)
	sp_plot <- as.matrix(sp_plot)
	n.plots <- ncol(sp_plot)
	if(is.null(rownames(sp_plot))) stop("Species-plot matrix must have row names corresponding to tip labels of the phylogeny.")
	n.nodes <- phy$Nnode
	out <- array(dim=c(length(sub), 4)) ### prep an output file
	for (i in 1:n.nodes) {
		sp_plot.i <- sp_plot[rownames(sp_plot) %in% sub[[i]]$tip.label,]
		if(nrow(sp_plot.i) > 0){
			Bst.i <- spacodi.calc(sp_plot.i, phy = sub[[i]])$Bst
		} else {
			Bst.i <- NA	
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

# generates a distribution of expected Bst values at each node based on randomization of (observed or simulated) individuals across plots and species 
# 'method' should be one of c("r.sp_plot","1s","1a","2s","2x","3i","3t","3x")
# parm=list(abund.class.ratio=value) for method="1a"
# parm=list(level=proportion) for method="2x" 
Bst.permutation <- function(phy, sp_plot, n.rep = 10, method="1s", parm=NULL) {
	require(geiger)
	sub=subtrees(phy)
	n.plot=ncol(sp_plot)
	n.taxa=Ntip(phy)
	out <- array(dim=c(length(sub), 4, n.rep)) ### prep an output file
	if(method=="r.sp_plot"){
		for(j in 1:n.rep){
			out[,,j] <- foo <- Bst.all.nodes(phy = phy, sp_plot = as.data.frame(r.sp_plot(sp_plot,n.rep=1)))
		}
	} else if(method=="1s"){
		for(j in 1:n.rep){
			out[,,j] <- foo <- Bst.all.nodes(phy = phy, sp_plot = resamp.1s(sp_plot))
		}
	} else if(method=="1a" && !is.null(parm$abund.class.ratio)){
		for(j in 1:n.rep){
			out[,,j] <- foo <- Bst.all.nodes(phy = phy, sp_plot = resamp.1a(sp_plot, parm$abund.class.ratio))
		}
	} else if(method=="2s"){
		for(j in 1:n.rep){
			out[,,j] <- foo <- Bst.all.nodes(phy = phy, sp_plot = resamp.2s(sp_plot))
		}
	} else if(method=="2x" && !is.null(parm$level)){
		for(j in 1:n.rep){
			out[,,j] <- foo <- Bst.all.nodes(phy = phy, sp_plot = resamp.2x(sp_plot, parm$level))
		}
	} else if(method=="3i"){
		for(j in 1:n.rep){
			out[,,j] <- foo <- Bst.all.nodes(phy = phy, sp_plot = resamp.3i(sp_plot))
		}
	} else if(method=="3t"){
		warning("Plots were assumed to be arranged linearly in space")
		for(j in 1:n.rep){
			out[,,j] <- foo <- Bst.all.nodes(phy = phy, sp_plot = resamp.3t(sp_plot))
		}
	} else if(method=="3x" && !is.null(parm$level)){
		for(j in 1:n.rep){
			out[,,j] <- foo <- Bst.all.nodes(phy = phy, sp_plot = resamp.3x(sp_plot, parm$level))
		}
	} else stop(cat("Unrecognized method declaration or insufficient information supplied:\n\tIs 'parm' set to a non-null value?\n\n"))
	dimnames(out)=list(NULL, names(as.data.frame(foo)), paste("iter",seq(1:n.rep),sep="."))
	return(out)
}

# randomly generates abundance data assuming log-normal abundance-distribution with the mean of the distribution determined by empirical data and number of empty cells maintained
# 'FIXME': add additional expected-abundance distributions
r.sp_plot<-function(sp_plot,n.rep=10){
	nulls=length(which(sp_plot==0))
	l.data=ncol(sp_plot)*nrow(sp_plot)
	non.nulls=l.data-nulls
	m.log=mean(log(sum(sp_plot)/non.nulls))
	out <- list() ### prep an output file
	for(j in 1:n.rep){
		out[[j]] <- as.data.frame(matrix(sample(c(round(rlnorm(n=non.nulls, meanlog=m.log)),rep(0,nulls))),nrow(sp_plot),ncol(sp_plot)))
		rownames(out[[j]]) <- row.names(sp_plot)
		names(out[[j]])=names(sp_plot)
	}
	return(out)
}

# implementation of Hardy '1a' resampling for a species-plot matrix (Hardy 2008 JoE). Shuffles species within abundance classes
resamp.1a <- function(input, abund.class.ratio = 4) {
	abund <- rowSums(input)
	n.ind   <- length(abund) 
	classes <- runif(1)* abund.class.ratio^(0:7)
	class <- rep(NA, n.ind)
	for(i in 1:6){
		class[abund > classes[i]] <- i
	}
	new_name <- rep(NA, n.ind) 
	for(i in unique(class)){
		new_name[class == i] <- sample(rownames(input[class == i,]))
	}
	row.names(input) <- new_name
	return(input[order(rownames(input)),])
}

# implementation of Hardy '1s' resampling for a species-plot matrix (Hardy 2008 JoE). Shuffles species all species, regardless of abundance
resamp.1s <- function(obj) {
	orig=obj
	row.names(obj) <- sample(row.names(obj))
	return(obj[order(match(row.names(obj),row.names(orig))),])
}

# implementation of Hardy '2s' resampling for a species-plot matrix (Hardy 2008 JoE). Shuffles all species, regardless of abundance, within plots
resamp.2s <- function(obj) {
	for(nn in 1:length(names(obj))){
		obj[,nn]=sample(obj[,nn])
	}
	return(obj)
}

# implementation of Hardy '2x' resampling for a species-plot matrix (Hardy 2008 JoE). Shuffles abundances within each of a pair of plots for a pair of species (i.e., Gotelli swapping)
resamp.2x <- function(obj, level=0.1) {
	swaps=round(level*ncol(obj)*nrow(obj))
	orig=obj
	for(swap in 1:swaps){
		rcol=sample(1:ncol(obj),2)
		rspp=sample(1:nrow(obj),2)
		obj[rspp[1],rcol[1]]=orig[rspp[2],rcol[1]]
		obj[rspp[2],rcol[1]]=orig[rspp[1],rcol[1]]
		obj[rspp[1],rcol[2]]=orig[rspp[2],rcol[2]]
		obj[rspp[2],rcol[2]]=orig[rspp[1],rcol[2]]
		orig=obj
	}
	if(sum(obj)!=sum(orig))warning("A poor result is likely.")
	return(obj)
}

# implementation of Hardy '3i' resampling for a species-plot matrix (Hardy 2008 JoE). Shuffles abundances within species among plots
resamp.3i <- function(obj) {
	for(ss in 1:nrow(obj)){
		obj[ss,]=sample(sp_plot[ss,])
	}
	return(obj)
}

# implementation of Hardy '3t' resampling for a species-plot matrix (Hardy 2008 JoE). Shuffles abundances to neighboring plots (with even spacing) for each species
resamp.3t <- function(obj) {
	for(ss in 1:nrow(obj)){
		torus=rep(1:ncol(obj),2)
		r=sample(1:(ncol(obj)-1),1)
		t.array=array(dim=c(ncol(obj),2))
		t.array[,1]=1:ncol(obj)
		for(o in 1:ncol(obj)){
			tt=torus[r+o]
			t.array[tt,2]=obj[ss,o]
		}
		obj[ss,]=t.array[,2]
	}
	return(obj)
}

# implementation of Hardy '3x' resampling for a species-plot matrix (Hardy 2008 JoE). Shuffles abundances within each of a pair of species for a pair of plots (i.e., Gotelli swapping)
resamp.3x <- function(obj, level=0.1) {
	swaps=round(level*ncol(obj)*nrow(obj))
	orig=obj
	for(swap in 1:swaps){
		rcol=sample(1:ncol(obj),2)
		rspp=sample(1:nrow(obj),2)
		obj[rspp[1],rcol[1]]=orig[rspp[1],rcol[2]]
		obj[rspp[1],rcol[2]]=orig[rspp[1],rcol[1]]
		obj[rspp[2],rcol[1]]=orig[rspp[2],rcol[2]]
		obj[rspp[2],rcol[2]]=orig[rspp[2],rcol[1]]
		orig=obj
	}
	if(sum(obj)!=sum(orig))warning("A poor result is likely.")
	return(obj)
}	

# calculate Blomberg K statistic for all subtrees of a phylogenetic tree; returns array with branching.times, descendants, Blomberg K, and node names
K.obs.nodes <- K.all.nodes <- function(phy, traits, return.all=TRUE){
	require(picante)
	if (!is.binary.tree(phy)) {
		phy=multi2di(phy)
		warning("Supplied tree was not dichotomous and was randomly resolved with ape:::multi2di(). ") 
	}
	if (length(phy$tip.label)!=nrow(traits)) stop("Tree cannot be matched to data")
	sub=subtrees(phy)
	n.traits <- ncol(traits)
	if(length(unique(colnames(traits))) < n.traits) {trait.names=paste("trt",seq(1:n.traits),sep=".")} else {trait.names = colnames(traits)}
	if(is.null(rownames(traits))||!all(row.names(traits)%in%phy$tip.label)) stop("Trait matrix must have row names corresponding to tip labels of the phylogeny.")
	n.nodes <- phy$Nnode
	out <- array(dim=c(length(sub), 4, n.traits)) ### prep an output file
	array.cols<-c("blomberg.K","tips","node.time","node.ID")
	for (i in 1:n.nodes) {### cycle over subtrees
		trt.i <- data.frame(traits[match(sub[[i]]$tip.label, rownames(traits)),]) #extract the relevant individual from trait matrix, put them into tiplabels order. 
		for (j in 1:n.traits){ 
			trt.ij <- trt.i[,j]
			nn<-sub[[i]]$Nnode
			if(length(trt.ij) > 0 & var(trt.ij) != 0 & nn > 1){ # eliminate basal polytomies and cases with no trait variance
				K.ij <- Kcalc(trt.ij, sub[[i]], F)
				} else {
				K.ij <- NA	
			}
			out[i, 1, j] <- K.ij
			out[i, 2, j] <- nn+1
			out[i, 3, j] <- max(branching.times(sub[[i]]))
			out[i, 4, j] <- min(sub[[i]]$node)
		}
		
	}
	dimnames(out)=list(NULL, array.cols, trait.names)
	if(!return.all) {
		out=out[which(!is.na(out[,1,1])),,]
		cat("Clades with two taxa were excluded as Blomberg's K is incalculable for small clades.\n\tIf this is unintended behavior, switch the 'return.all' operator to 'true'.\n\n")
	}
	return(out)
}


# generates a distribution of expected K values for each node, using a given phylogeny, and either brownian-motion trait values (Kexp = 1) or random trait values (Kexp = 0)
K.exp.nodes <- function(phy, n.rep = 10, brownian = T,...) {
	require(geiger)
	out <- matrix(ncol = n.rep, nrow = Nnode(phy))
	trt <- matrix(ncol = n.rep, nrow = length(phy$tip.label))
	if(brownian == T){
			trt.sims <- as.data.frame(sim.char(phy = phy, as.matrix(1), n.rep))
			names(trt.sims)=paste("trt",seq(1:n.rep),sep=".")
		} else {
			trt.sims <- as.data.frame(replicate(n.rep, runif(Ntip(phy))))
			row.names(trt.sims)=phy$tip.label
			names(trt.sims)=paste("trt",seq(1:n.rep),sep=".")
		}
	out <- K.obs.nodes(phy = phy, traits = trt.sims,...)
	return(out)
}

### UTILITIES ###
# returns nodes that are within some fraction of time from the root of the tree
node.time.extractor<-function(phy,time.fraction=0.5,return.times=F) {
	require(ape)
	xx=branching.times(phy)
	xx=abs(xx-max(xx))
	
	if(any(attr(phy,"names")=="node.label")) {nn=names(xx)[keep<-which(xx<=time.fraction*max(xx))] 
	} else {nn=as.numeric(names(xx)[keep<-which(xx<=time.fraction*max(xx))])}
	
	if(return.times==TRUE) {return(list(nn,unname(xx[keep])))
	} else {return(nn, as.numeric(nn))}
}

# returns nodes sorted by time such that the first element of the list is the root, last entry is most recent split in tree
node.time.sort<-function(phy,return.times=F) {
	require(ape)
	xx=branching.times(phy)
	xx=abs(xx-max(xx))
	xx=xx[order(xx)]
	if(return.times==FALSE) {
		return(list(as.numeric(names(xx))))
	} else {return(xx)}
}

# gets rid of last split (of zero length) when using bd.tree() or geiger:::birthdeath.tree() with 'taxa.stop' criterion; chooses one of the two lineages to prune from tree
prune.last.split<-function(phy) {
	require(ape)
	nn=as.numeric(names(xx<-branching.times(phy))[which(xx==min(xx))])
	for(node in 1:length(nn)) {	
		sp=sample(phy$edge[which(phy$edge[,1]==nn[node]),2],1)
		obj=drop.tip(phy,phy$tip.label[sp])
	}
	return(obj)
}

# simulates species-by-plot abundance data
r.plot<-function(species=100, plots=30, missing.prop=0.6, mean.abund=15, sim.tree=FALSE){
	require(geiger)
	l.data=species*plots
	nulls=round(l.data*missing.prop)+1
	non.nulls=l.data-nulls
	m.log=log(mean.abund)
	out <- as.data.frame(matrix(sample(c(round(rlnorm(n=non.nulls, meanlog=m.log)),rep(0,nulls))),species,plots))
	rownames(out) <- paste("sp",seq(1:species),sep="")
	names(out) <- paste("plot",seq(1:plots),sep="")
	if(sim.tree) {
		phy=prune.last.split(bd.tree(b=0.02,d=0.0,taxa.stop=(1+species),return.all.extinct=F))
		phy$tip.label=row.names(out)
		return(list(out,phy))
	} else {return(out)}
}


# convert to phylocom format
as.phylocom<-function(data, outfile=NULL){
	if(class(data)!="data.frame"){
		if(class(data)=="matrix") data=as.data.frame(data) else if(class(data)=="vector") data=as.data.frame(data) else stop("Function cannot handle data that are not in data.frame format")
	}
	if(length(unique(row.names(data)))!=nrow(data))warning("Data do not appear to be in proper format.  Results may be nonsensical.")
	plots=ncol(data)
	spp=nrow(data)
	out=array(dim=c(spp*plots,3))
	for(plot in 1:plots){
		start=((plot-1)*spp)+1
		end=plot*spp
		out[start:end,]=cbind(names(data)[plot], data[,plot], row.names(data))
	}
	out=as.data.frame(out)
	names(out)=c("plot","samples","species")
	if(!is.null(outfile)){write.table(out,outfile,quote=F,row=F,col=F,sep="\t")}
	return(out)
}

# convert to SPaCODi format
as.spacodi<-function(data, outfile=NULL){
	if(ncol(data)!=3)stop(cat("Format does not appear to be in 'triplet'.\n\tColumns should be plot, sample(s), species.\n\n"))
	if(class(data)!="data.frame"){
		if(class(data)=="matrix") data=as.data.frame(data) else if(class(data)=="vector") data=as.data.frame(data) else stop("Function cannot handle data that are not in data.frame format")
	}
	species=as.character(unique(data[,3]))
	dd=split(data,data[,1])
	out=array(dim=c(length(species)->spp,length(dd)->plots))
	for(plot in 1:plots){
		cur.array=array(dim=c(spp, 1))
		ind=as.numeric(as.vector(dd[[plot]][,2]))
		dd[[plot]]->cur.plot
		names(ind)=as.character(cur.plot[,3])
		for(r in 1:nrow(cur.plot)){
			cur.array[which(names(ind[r])==species)]=ind[r]
		}
		out[,plot]=as.numeric(cur.array)
	}
	out=as.data.frame(out)
	names(out)=names(dd)
	row.names(out)=species
	if(!is.null(outfile)){write.spacodi.data(out,outfile)}
	return(out)
}

# exports data.frame of species-by-plot abundance data in format readable by 'spacodi.exe'								
write.spacodi.data<-function(data, outfile){
	if(file.exists(outfile)){
		warning("Overwrote existing outfile.")
		unlink(outfile)
	}
	names=names(data)
	for(n in 1:length(names)){cat(c("\t",names[n]),file=outfile,append=T,sep="")}
	cat("\n",file=outfile,append=T,sep="")
	write.table(data,outfile,quote=F,col=F,append=T,sep="\t")
}

# simulates a phylogeny under time-homogenous birth-death process.  Modified from geiger:::birthdeath.tree()
bd.tree <- function (b=1, d=0, time.stop = 0, taxa.stop = 15, return.all.extinct = TRUE) {
	if (time.stop == 0 & taxa.stop == 0) 
	stop("Must have stopping criterion\n")
    while (1) {
        edge <- rbind(c(1, 2), c(1, 3))
        edge.length <- rep(NA, 2)
        stem.depth <- numeric(2)
        alive <- rep(TRUE, 2)
        t <- 0
        next.node <- 4
        repeat {
            if (taxa.stop) 
			if (sum(alive) >= taxa.stop) 
			break
            if (sum(alive) == 0) 
			break
            dt <- rexp(1, sum(alive) * (b + d))
            t <- t + dt
            if (time.stop) 
			if (t >= time.stop) {
				t <- time.stop
				break
			}
            r <- runif(1)
            if (r <= b/(b + d)) {
                random_lineage <- round(runif(1, min = 1, max = sum(alive)))
                e <- matrix(edge[alive, ], ncol = 2)
                parent <- e[random_lineage, 2]
                alive[alive][random_lineage] <- FALSE
                edge <- rbind(edge, c(parent, next.node), c(parent, 
															next.node + 1))
                next.node <- next.node + 2
                alive <- c(alive, TRUE, TRUE)
                stem.depth <- c(stem.depth, t, t)
                x <- which(edge[, 2] == parent)
                edge.length[x] <- t - stem.depth[x]
                edge.length <- c(edge.length, NA, NA)
            }
            else {
                random_lineage <- round(runif(1, min = 1, max = sum(alive)))
                edge.length[alive][random_lineage] <- t - stem.depth[alive][random_lineage]
                alive[alive][random_lineage] <- FALSE
            }
        }
        if (return.all.extinct == T | sum(alive) > 1) 
		break
    }
    edge.length[alive] <- t - stem.depth[alive]
    n <- -1
    for (i in 1:max(edge)) {
        if (any(edge[, 1] == i)) {
            edge[which(edge[, 1] == i), 1] <- n
            edge[which(edge[, 2] == i), 2] <- n
            n <- n - 1
        }
    }
    edge[edge > 0] <- 1:sum(edge > 0)
    tip.label <- 1:sum(edge > 0)
    mode(edge) <- "character"
    mode(tip.label) <- "character"
    obj <- list(edge = edge, edge.length = edge.length, tip.label = tip.label)
    class(obj) <- "phylo"
    obj <- old2new.phylo(obj)
    obj <- read.tree(text = write.tree(obj))
    obj$tip.label=as.character(paste("sp",seq(1,length(obj$tip.label)),sep=""))
	obj
}

### PLOTTING FUNCTIONS ###

#plot.quant plots points at the tips sized by the value of the datum; data should be nameless array
plot.quant<-function(phy, data, data.names, shrinker=1, cex=0.5, adj=c(15,5), lwd=0.5, font=1){
	
	f=match(phy$tip.label, data.names)
	if(is.na(sum(f))) {
		stop("tips cannot be matched to tree")
	}
	else {
		data=data[f]
		phy$trait=data
		plot.phylo(phy, show.tip.label=T, label.offset=adj[1], cex=cex, font=font)
		tiplabels(text=NULL, pch=21, adj=adj[2], cex=c(phy$trait/shrinker), bg="white",lwd=lwd)
		} 
}

# Calculate Hardy Bst for all subnodes of a phylogeny.. plot Bst through time
plot.Bst.through.time <- function(phy, sp_plot, outfile=NULL){
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

# plots color-coded Bst on tree for every node
plot.Bst <- function(phy, sp_plot, outfile = NULL){
	B=data.frame(Bst.all.nodes(phy,sp_plot,return.all=T)[,1])
	if(!is.null(outfile)){
		pdf(outfile)
	}
	par(mfrow = c(1, ncol(B)))
	for(i in 1:ncol(B)){
		bst.i <- B[,i]
		col.i <- 1+(bst.i-min(bst.i, na.rm = T))*99/(max(bst.i, na.rm = T)-min(bst.i, na.rm = T))
		col.i <- round(col.i)
		cols  <- topo.colors(100)
		col.i <- cols[c(col.i, recursive = T)]
		col.i[is.na(col.i)] <- "gray"
		if(i != ncol(B)){
			plot(phy, show.tip.label = F, no.margin = T)
		} else {
			plot(phy, cex = 0.2, no.margin = T)
		}
		nodelabels.tim(bg = col.i, border = col.i, text = "")
		legend("bottomleft", legend = sprintf("%1.2f", seq(min(bst.i, na.rm = T), max(bst.i, na.rm = T), length = 10)), fill = topo.colors(10), bg = "white")
	}
	title(main = "Values of Bst", outer = T, line = -1)
	if(!is.null(outfile)){
		dev.off()
	}
}

# Calculate Blomberg K for all subnodes of a phylogeny.. plot K through time
plot.K.through.time <- function(phy, traits, outfile=NULL){
	if(!is.null(outfile)){
		pdf(outfile)
	}
	if(ncol(traits)>1)stop("Cannot plot more than a single trait at a time")
	as.data.frame(K.obs.nodes(phy,traits,return.all=F))->raw
	raw$bt=-raw$node.time
	plot(raw$blomberg.K~raw$bt,ylab="Blomberg's K",xlab="branching.time", main="Blomberg's K through time")
	if(!is.null(outfile)){
		dev.off()
	}
}

# plots color-coded K on tree for every node
plot.K <- function(phy, traits, outfile = NULL){
	K=data.frame(K.all.nodes(phy,traits,return.all=T)[,1,])
	names(K)=names(traits)
	if(!is.null(outfile)){
		pdf(outfile)
	}
	par(mfrow = c(1, ncol(K)))
	for(i in 1:ncol(K)){
		trt.i <- K[,i]
		col.i <- 1+(trt.i-min(trt.i, na.rm = T))*99/(max(trt.i, na.rm = T)-min(trt.i, na.rm = T))
		col.i <- round(col.i)
		cols  <- topo.colors(100)
		col.i <- cols[c(col.i, recursive = T)]
		col.i[is.na(col.i)] <- "gray"
		if(i != ncol(K)){
			plot(phy, show.tip.label = F, no.margin = T)
		} else {
			plot(phy, cex = 0.2, no.margin = T)
		}
		nodelabels.tim(bg = col.i, border = col.i, text = "")
		legend("bottomleft", legend = sprintf("%1.2f", seq(min(trt.i, na.rm = T), max(trt.i, na.rm = T), length = 10)), fill = topo.colors(10), bg = "white")
		title(main = names(K)[i], line = -2)
	}
	title(main = "Values of Blomberg's K", outer = T, line = -1)
	if(!is.null(outfile)){
		dev.off()
	}
}

# for plot.K()
nodelabels.tim <- function (text, node, adj = c(0.5, 0.5), frame = "rect", pch = NULL, thermo = NULL, pie = NULL, piecol = NULL, col = "black", bg = "lightblue", border = "black", ...) {
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    if (missing(node)) 
	node <- (lastPP$Ntip + 1):length(lastPP$xx)
    XX <- lastPP$xx[node]
    YY <- lastPP$yy[node]
    BOTHlabels.tim(text, node, XX, YY, adj, frame, pch, thermo, pie, piecol, col, bg, border, ...)
}

# for plot.K()
BOTHlabels.tim <- function (text, sel, XX, YY, adj, frame, pch, thermo, pie, piecol, col, bg, border, ...) {
    if (missing(text)) 
	text <- NULL
    if (length(adj) == 1) 
	adj <- c(adj, 0.5)
    if (is.null(text) && is.null(pch) && is.null(thermo) && is.null(pie)) 
	text <- as.character(sel)
    frame <- match.arg(frame, c("rect", "circle", "none"))
    args <- list(...)
    CEX <- if ("cex" %in% names(args)) 
	args$cex
    else par("cex")
    if (frame != "none" && !is.null(text)) {
        if (frame == "rect") {
            width <- strwidth(text, units = "inches", cex = CEX)
            height <- strheight(text, units = "inches", cex = CEX)
            if ("srt" %in% names(args)) {
                args$srt <- args$srt%%360
                if (args$srt == 90 || args$srt == 270) {
					tmp <- width
					width <- height
					height <- tmp
                }
                else if (args$srt != 0) 
				warning("only right angle rotation of frame is supported;\n         try  `frame = \"n\"' instead.\n")
            }
            width <- xinch(width)
            height <- yinch(height)
            xl <- XX - xinch(0.10)
            xr <- XX + xinch(-0.01)
            yb <- YY - yinch(0.02)
            yt <- YY + yinch(0.02)
            rect(xl, yb, xr, yt, col = bg, border = border)
        }
        if (frame == "circle") {
            radii <- 0.8 * apply(cbind(strheight(text, units = "inches", 
												 cex = CEX), strwidth(text, units = "inches", 
																	  cex = CEX)), 1, max)
            symbols(XX, YY, circles = radii, inches = max(radii), add = TRUE, bg = bg)
        }
    }
    if (!is.null(thermo)) {
        parusr <- par("usr")
        width <- CEX * (parusr[2] - parusr[1])/40
        height <- CEX * (parusr[4] - parusr[3])/15
        if (is.vector(thermo)) 
		thermo <- cbind(thermo, 1 - thermo)
        thermo <- height * thermo
        xl <- XX - width/2
        xr <- xl + width
        yb <- YY - height/2
        if (is.null(piecol)) 
		piecol <- rainbow(ncol(thermo))
        rect(xl, yb, xr, yb + thermo[, 1], border = NA, col = piecol[1])
        for (i in 2:ncol(thermo)) rect(xl, yb + rowSums(thermo[, 
														1:(i - 1), drop = FALSE]), xr, yb + rowSums(thermo[, 
																									1:i]), border = NA, col = piecol[i])
        rect(xl, yb, xr, yb + height, border = "black")
        segments(xl, YY, xl - width/5, YY)
        segments(xr, YY, xr + width/5, YY)
    }
    if (!is.null(pie)) {
        if (is.vector(pie)) 
		pie <- cbind(pie, 1 - pie)
        xrad <- CEX * diff(par("usr")[1:2])/50
        xrad <- rep(xrad, length(sel))
        for (i in 1:length(sel)) floating.pie.asp(XX[i], YY[i], 
												  pie[i, ], radius = xrad[i], col = piecol)
    }
    if (!is.null(text)) 
	text(XX, YY, text, adj = adj, col = col, ...)
    if (!is.null(pch)) 
	points(XX + adj[1] - 0.5, YY + adj[2] - 0.5, pch = pch, 
		   col = col, bg = bg, ...)
}

