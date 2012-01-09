resamp.3i <-
function(obj) {
	for(ss in 1:nrow(obj)){
		obj[ss,]=sample(obj[ss,])
	}
	return(obj)
}

