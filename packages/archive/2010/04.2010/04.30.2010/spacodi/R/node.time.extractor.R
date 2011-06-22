node.time.extractor <-
function(phy, start.time=0, stop.time=1, return.times=TRUE, proportion=TRUE, ape.ID=TRUE) {
	flag=FALSE
	if(!is.null(phy$node.label)) {
		nodes.orig=phy$node.label
		phy$node.label=NULL
		flag=TRUE
	} 
	xx = branching.times(phy)
	if (proportion && any(c(start.time, stop.time) > 1)) stop("Times do not appear to be proportions; at least one element exceeds 1.")
	if (proportion) times = c(start.time * max(xx), stop.time * max(xx)) else times = c(start.time, stop.time)
	stt = times[m <- which(c(start.time, stop.time) == min(c(start.time, stop.time)))]
	stp = times[-m]
	if (stp > max(xx)) {
		stp = max(xx)
		warning("Supplied time was not sensible: root age was substituted as a bound on the time-slice.")
	}
	new <- xx[which(xx <= stp)]
	old <- xx[which(xx >= stt)]
	nn <- as.numeric(intersect(names(new), names(old)))
	keep <- intersect(new, old)
	names(keep) = nn
	srt = sort(keep, decreasing = TRUE)
	if (return.times == TRUE) {
		if(flag && !ape.ID) {
			foo = list(nodes.orig[match(names(srt), names(xx))], unname(keep))
		} else {
			foo = list(as.numeric(names(srt)), unname(keep))
		}
		names(foo) = c("nodes", "times")
		return(foo)
	}
	else {
		if(flag && !ape.ID) {
			foo = nodes.orig[match(names(srt), names(xx))]
		} else {
			foo = as.numeric(names(srt))
		}
		return(foo)
	}
}


