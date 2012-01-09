resamp.1a <-
function(obj, abund.class.ratio = 4) {
	abund <- rowSums(obj)
	n.ind   <- length(abund) 
	classes <- runif(1)* abund.class.ratio^(0:7)
	class <- rep(NA, n.ind)
	for(i in 1:6){
		class[abund > classes[i]] <- i
	}
	new_name <- rep(NA, n.ind) 
	for(i in unique(class)){
		new_name[class == i] <- sample(rownames(obj[class == i,]))
	}
	row.names(obj) <- new_name
	return(obj[order(rownames(obj)),])
}

