spacodi.tree.plot <-function(spacodi.object, phy, sig.plot=TRUE, cut.off=0.05, cols=list("white", "gray", "black"), main=TRUE, outfile=NULL, ...){
	op=par(no.readonly=TRUE)
	on.exit(par(op))
	
	if(class(spacodi.object)=="list") base=as.data.frame(spacodi.object[[1]]) else base=as.data.frame(spacodi.object)
	if(nrow(base)!=phy$Nnode)stop("SPACoDi.plot requires data (whether NA not) for every node in the spacodi.object.") 
	
	base=base[order(base$node.ID),]
	
	if(sig.plot==TRUE) {
		rand=as.data.frame(spacodi.object$randomization.test)
		if(length(rand)!=0) {
			rand=rand 
		} else {
			stop("Cannot perform significance plot without there being a randomization.test within the spacodi.object")
		}
	}
	
	if(any(names(base)=="Bst")) {
	   # BST PLOT
	   B=as.matrix(base$Bst)
	   if(!is.null(outfile))pdf(outfile)
	   par(mfrow = c(1, ncol(B)))
	   for(i in 1:ncol(B)){
			bst.i <- B[,i]
			if(sig.plot==TRUE) {
				col.i=array(dim=c(nrow(base), 2))
				col.i[,1]=as.numeric(base$node.ID)
				for(x in 1:nrow(base)){
					if(!is.na(match(base$node.ID[x],rand$node.ID))) {
						r.row=rand[which(rand$node.ID==base$node.ID[x]),]
						if((r.row$obs.Bst-rand$m.exp.Bst) > 0 && r.row$p.value <= cut.off) {
							col.i[x,2]=cols[[2]]
						} else if((r.row$obs.Bst-r.row$m.exp.Bst) < 0 && r.row$p.value <= cut.off) {
							col.i[x,2]=cols[[3]]
						} else {
							col.i[x,2]=NA
						}
					}
					else {col.i[x,2]=NA}
				}
				col.i=as.vector(col.i[,2])
			} else {
				col.i <- 1+(bst.i-min(bst.i, na.rm = TRUE))*99/(max(bst.i, na.rm = TRUE)-min(bst.i, na.rm = TRUE))
				col.i <- round(col.i)
				cols  <- topo.colors(100)
				col.i <- cols[c(col.i, recursive = TRUE)]
				col.i[is.na(col.i)] <- "gray"
			}
			if(i != ncol(B)){
				plot(phy, show.tip.label = FALSE, no.margin = TRUE, ...)
			} else {
				plot(phy, cex = 0.2, no.margin = TRUE, ...)
			}
		    if(sig.plot) {
				nodelabels(bg=col.i, col="gray", frame="circle", pch=21, cex=1.5)
				legend("bottomleft", legend=c("Bst: as expected", "Bst: larger than expected", "Bst: smaller than expected"), fill=c(NA, cols[[2]],cols[[3]]), bg = NULL, box.lty="blank", border="gray", cex=0.5)
		    } else {
				nodelabels.sp(bg = col.i, border = col.i, text = "")
				legend("bottomleft", legend = sprintf("%1.2f", seq(min(bst.i, na.rm = TRUE), max(bst.i, na.rm = TRUE), length = 10)), fill = topo.colors(10), bg = "white")
		    }
	   }
		if(main)title(main = "Values of Bst", outer = TRUE, line = -1)
		if(!is.null(outfile))dev.off()
	} else if(any(names(base)=="Blomberg.K")) {
		# K PLOT
		K=as.matrix(base$Blomberg.K)
		
		if(!is.null(outfile)){
			pdf(outfile)
		}
		par(mfrow = c(1, ncol(K)))
		for(i in 1:ncol(K)){
			trt.i <- K[,i]
			if(sig.plot==TRUE) {
				col.i=array(dim=c(nrow(base), 2))
				col.i[,1]=as.numeric(base$node.ID)
				for(x in 1:nrow(base)){
					if(!is.na(match(base$node.ID[x],rand$node.ID))) {
						r.row=rand[which(rand$node.ID==base$node.ID[x]),]
						if((r.row$obs.K-rand$m.exp.K) > 0 && r.row$p.value <= cut.off) {
							col.i[x,2]=cols[[2]]
						} else if((r.row$obs.K-r.row$m.exp.K) < 0 && r.row$p.value <= cut.off) {
							col.i[x,2]=cols[[3]]
						} else {
							col.i[x,2]=NA
						}
					}
					else {
						col.i[x,2]=NA
					}
				}
				col.i=as.vector(col.i[,2])
			} else {
				col.i <- 1+(trt.i-min(trt.i, na.rm = TRUE))*99/(max(trt.i, na.rm = TRUE)-min(trt.i, na.rm = TRUE))
				col.i <- round(col.i)
				cols  <- topo.colors(100)
				col.i <- cols[c(col.i, recursive = TRUE)]
				col.i[is.na(col.i)] <- "gray"				
			}
			if(i != ncol(K)){
				plot(phy, show.tip.label = FALSE, no.margin = TRUE, ...)
			} else {
				plot(phy, cex = 0.2, no.margin = TRUE, ...)
			}
			if(sig.plot) {
				nodelabels(bg=col.i, col="gray", frame="circle", pch=21, cex=1.5)
				legend("bottomleft", legend=c("K: as expected", "K: larger than expected", "K: smaller than expected"), fill=c(NA, cols[[2]],cols[[3]]), bg = NULL, box.lty="blank", border="gray", cex=0.5)
			} else {
				nodelabels.sp(bg = col.i, border = col.i, text = "")
				legend("bottomleft", legend = sprintf("%1.2f", seq(min(trt.i, na.rm = TRUE), max(trt.i, na.rm = TRUE), length = 10)), fill = topo.colors(10), bg = "white")
			}
			if(main)title(main = names(K)[i], line = -2)
		}
		title(main = "Values of Blomberg's K", outer = TRUE, line = -1)
		if(!is.null(outfile))dev.off()
	} else {stop("Cannot interpret spacodi.object for plotting.")}
	invisible()
}	

