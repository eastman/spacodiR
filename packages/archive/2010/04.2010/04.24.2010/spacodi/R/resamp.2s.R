resamp.2s <-
function(obj) {
	for(nn in 1:length(names(obj))){
		obj[,nn]=sample(obj[,nn])
	}
	return(obj)
}

