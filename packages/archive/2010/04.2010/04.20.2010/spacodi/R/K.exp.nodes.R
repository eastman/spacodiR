K.exp.nodes <-
function(phy, n.rep = 10, brownian = TRUE,...) {
	require(geiger)
	out <- matrix(ncol = n.rep, nrow = Nnode(phy))
	trt <- matrix(ncol = n.rep, nrow = length(phy$tip.label))
	if(brownian == TRUE){
			trt.sims <- as.data.frame(sim.char(phy = phy, as.matrix(1), n.rep))
			names(trt.sims)=paste("trt",seq(1:n.rep),sep=".")
		} else {
			trt.sims <- as.data.frame(replicate(n.rep, runif(Ntip(phy))))
			row.names(trt.sims)=phy$tip.label
			names(trt.sims)=paste("trt",seq(1:n.rep),sep=".")
		}
	out <- K.all.nodes(phy = phy, traits = trt.sims,...)
	return(out)
}

