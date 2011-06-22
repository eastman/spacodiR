node.time.extractor <-
function(phy,time.fraction=0.5,return.times=FALSE) {
	require(ape)
	xx=branching.times(phy)
	xx=abs(xx-max(xx))
	
	if(any(attr(phy,"names")=="node.label")) {nn=names(xx)[keep<-which(xx<=time.fraction*max(xx))] 
	} else {nn=as.numeric(names(xx)[keep<-which(xx<=time.fraction*max(xx))])}
	
	if(return.times==TRUE) {return(list(nn,unname(xx[keep])))
	} else {return(nn, as.numeric(nn))}
}

