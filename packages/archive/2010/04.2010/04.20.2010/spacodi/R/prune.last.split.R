prune.last.split <-
function(phy) {
	require(ape)
	nn=as.numeric(names(xx<-branching.times(phy))[which(xx==min(xx))])
	for(node in 1:length(nn)) {	
		sp=sample(phy$edge[which(phy$edge[,1]==nn[node]),2],1)
		obj=drop.tip(phy,phy$tip.label[sp])
	}
	return(obj)
}

