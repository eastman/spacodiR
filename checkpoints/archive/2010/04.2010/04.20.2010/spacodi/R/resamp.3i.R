resamp.3i <-
function(obj) {
	for(ss in 1:nrow(obj)){
		obj[ss,]=sample(sp_plot[ss,])
	}
	return(obj)
}

