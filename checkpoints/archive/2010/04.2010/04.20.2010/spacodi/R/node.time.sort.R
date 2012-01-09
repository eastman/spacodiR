node.time.sort <-
function(phy,return.times=FALSE) {
	require(ape)
	xx=branching.times(phy)
	xx=abs(xx-max(xx))
	xx=xx[order(xx)]
	if(return.times==FALSE) {
		return(list(as.numeric(names(xx))))
	} else {return(xx)}
}

