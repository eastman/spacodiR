resamp.1a <-
function(obj, abund.class.ratio = 4) {
	if(abund.class.ratio<=1)stop("Supplied abundance-class ratio did not appear sensible: choose a value greater than 1.")
	orig=obj
	abund <- rowSums(obj)
	n.spp   <- length(abund)
	
	aa <- abund.class.ratio
	
	while(1) {
		classes <- aa^(0:ifelse(aa<2, n.spp*(1/(aa-1)), n.spp))
		if(length(classes)>=n.spp) {
			classes <- unique(round(runif(1)*classes))
			break()
		}
	}
	
	class <- rep(NA, n.spp)
	for(i in 1:length(classes)){
		class[abund > classes[i]] <- i
	}
	if(any(is.na(class))) class[which(is.na(class))]=1
	new_name <- rep(NA, n.spp) 
	for(i in unique(class)){
		new_name[class == i] <- sample(rownames(obj[class == i,]))
	}
	row.names(obj) <- new_name
	return(obj[order(match(row.names(obj),row.names(orig))),])
}