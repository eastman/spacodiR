Bst.permutation <-
function(phy, sp_plot, n.rep = 10, method="1s", parm=NULL) {
	require(geiger)
	sub=subtrees(phy)
	n.plot=ncol(sp_plot)
	n.taxa=Ntip(phy)
	out <- array(dim=c(length(sub), 4, n.rep)) ### prep an output file
	if(method=="r.sp_plot"){
		for(j in 1:n.rep){
			out[,,j] <- foo <- Bst.all.nodes(phy = phy, sp_plot = as.data.frame(r.sp_plot(sp_plot,n.rep=1)))
		}
	} else if(method=="1s"){
		for(j in 1:n.rep){
			out[,,j] <- foo <- Bst.all.nodes(phy = phy, sp_plot = resamp.1s(sp_plot))
		}
	} else if(method=="1a" && !is.null(parm$abund.class.ratio)){
		for(j in 1:n.rep){
			out[,,j] <- foo <- Bst.all.nodes(phy = phy, sp_plot = resamp.1a(sp_plot, parm$abund.class.ratio))
		}
	} else if(method=="2s"){
		for(j in 1:n.rep){
			out[,,j] <- foo <- Bst.all.nodes(phy = phy, sp_plot = resamp.2s(sp_plot))
		}
	} else if(method=="2x" && !is.null(parm$level)){
		for(j in 1:n.rep){
			out[,,j] <- foo <- Bst.all.nodes(phy = phy, sp_plot = resamp.2x(sp_plot, parm$level))
		}
	} else if(method=="3i"){
		for(j in 1:n.rep){
			out[,,j] <- foo <- Bst.all.nodes(phy = phy, sp_plot = resamp.3i(sp_plot))
		}
	} else if(method=="3t"){
		warning("Plots were assumed to be arranged linearly in space")
		for(j in 1:n.rep){
			out[,,j] <- foo <- Bst.all.nodes(phy = phy, sp_plot = resamp.3t(sp_plot))
		}
	} else if(method=="3x" && !is.null(parm$level)){
		for(j in 1:n.rep){
			out[,,j] <- foo <- Bst.all.nodes(phy = phy, sp_plot = resamp.3x(sp_plot, parm$level))
		}
	} else stop(cat("Unrecognized method declaration or insufficient information supplied:\n\tIs 'parm' set to a non-null value?\n\n"))
	dimnames(out)=list(NULL, names(as.data.frame(foo)), paste("iter",seq(1:n.rep),sep="."))
	return(out)
}

