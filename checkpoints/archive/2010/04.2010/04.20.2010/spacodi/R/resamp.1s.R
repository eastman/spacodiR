resamp.1s <-
function(obj) {
	orig=obj
	row.names(obj) <- sample(row.names(obj))
	return(obj[order(match(row.names(obj),row.names(orig))),])
}

