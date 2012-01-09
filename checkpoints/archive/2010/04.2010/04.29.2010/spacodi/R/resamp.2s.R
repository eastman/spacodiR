resamp.2s <-
function(obj) {
	for(nn in 1:ncol(obj)){
		obj[,nn]=sample(obj[,nn])
	}
	return(obj)
}

