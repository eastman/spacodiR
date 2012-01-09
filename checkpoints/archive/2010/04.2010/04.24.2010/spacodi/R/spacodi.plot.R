spacodi.plot <-
function(phy=NULL, sp_plot=NULL, trait=NULL, color.plot=FALSE, time.plot=FALSE, outfile=NULL, ...){
	op=par(no.readonly=TRUE)
	on.exit(par(op))
	if((color.plot && time.plot) || (!is.null(sp_plot) && !is.null(trait)) || (!time.plot && !color.plot)) {
	   stop("Choose either color.plot or time.plot; must also supply either an sp_plot or trait along with the phylogeny.")
	}
	if(!is.null(phy) && !is.null(sp_plot) && color.plot) {
	   # BST PLOT
	   B=data.frame(Bst.all.nodes(phy=phy,sp_plot=sp_plot,return.all=TRUE)[,1])
	   if(!is.null(outfile)){pdf(outfile)}
	   par(mfrow = c(1, ncol(B)))
	   for(i in 1:ncol(B)){
			bst.i <- B[,i]
			col.i <- 1+(bst.i-min(bst.i, na.rm = TRUE))*99/(max(bst.i, na.rm = TRUE)-min(bst.i, na.rm = TRUE))
			col.i <- round(col.i)
			cols  <- topo.colors(100)
			col.i <- cols[c(col.i, recursive = TRUE)]
			col.i[is.na(col.i)] <- "gray"
			if(i != ncol(B)){
				plot(phy, show.tip.label = FALSE, no.margin = TRUE, ...)
			} else {
				plot(phy, cex = 0.2, no.margin = TRUE, ...)
			}
			nodelabels.sp(bg = col.i, border = col.i, text = "")
			legend("bottomleft", legend = sprintf("%1.2f", seq(min(bst.i, na.rm = TRUE), max(bst.i, na.rm = TRUE), length = 10)), fill = topo.colors(10), bg = "white")
		}
		title(main = "Values of Bst", outer = TRUE, line = -1)
		if(!is.null(outfile)){dev.off()}
	} else if(!is.null(phy) && !is.null(sp_plot) && time.plot) {
		# BST PLOT
		if(!is.null(outfile)){
			pdf(outfile)
		}
		as.data.frame(Bst.all.nodes(phy=phy,sp_plot=sp_plot))->raw
		raw$bt=-raw$node.time
		plot(raw$Bst~raw$bt,ylab="Bst",xlab="branching.time", main="Bst through time", ...)
		if(!is.null(outfile)){
			dev.off()
		}
	} else if(!missing(trait) && !missing(phy) && color.plot) {
		# K PLOT
		print(trait)
		K=data.frame(K.all.nodes(phy=phy,traits=trait,return.all=TRUE)[,1,])
		names(K)=names(trait)
		print(K)
		if(!is.null(outfile)){
			pdf(outfile)
		}
		par(mfrow = c(1, ncol(K)))
		for(i in 1:ncol(K)){
			trt.i <- K[,i]
			col.i <- 1+(trt.i-min(trt.i, na.rm = TRUE))*99/(max(trt.i, na.rm = TRUE)-min(trt.i, na.rm = TRUE))
			col.i <- round(col.i)
			cols  <- topo.colors(100)
			col.i <- cols[c(col.i, recursive = TRUE)]
			col.i[is.na(col.i)] <- "gray"
			if(i != ncol(K)){
				plot(phy, show.tip.label = FALSE, no.margin = TRUE, ...)
			} else {
				plot(phy, cex = 0.2, no.margin = TRUE, ...)
			}
			nodelabels.sp(bg = col.i, border = col.i, text = "")
			legend("bottomleft", legend = sprintf("%1.2f", seq(min(trt.i, na.rm = TRUE), max(trt.i, na.rm = TRUE), length = 10)), fill = topo.colors(10), bg = "white")
			title(main = names(K)[i], line = -2)
		}
		title(main = "Values of Blomberg's K", outer = TRUE, line = -1)
		if(!is.null(outfile)){
			dev.off()
		}
	} else if (!missing(trait) && !missing(phy) && time.plot) { 
		# K PLOT
		if(!is.null(outfile)){
			pdf(outfile)
		}
		if(is.atomic(trait))stop("Cannot match traits with phylogeny: supply trait as a dataframe with row.names as species")
		if(ncol(trait)>1)stop("Cannot plot more than a single trait at a time.")
		as.data.frame(K.all.nodes(phy=phy,traits=trait,return.all=FALSE))->raw
		raw$bt=-raw$node.time
		plot(raw$blomberg.K~raw$bt,ylab="Blomberg's K",xlab="branching.time", main="Blomberg's K through time", ...)
		if(!is.null(outfile)){
			dev.off()
		}
	}
	invisible()
}	

