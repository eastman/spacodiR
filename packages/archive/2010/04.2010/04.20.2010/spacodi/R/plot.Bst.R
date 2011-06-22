plot.Bst <-
function(phy, sp_plot, outfile = NULL){
	B=data.frame(Bst.all.nodes(phy,sp_plot,return.all=TRUE)[,1])
	if(!is.null(outfile)){
		pdf(outfile)
	}
	par(mfrow = c(1, ncol(B)))
	for(i in 1:ncol(B)){
		bst.i <- B[,i]
		col.i <- 1+(bst.i-min(bst.i, na.rm = TRUE))*99/(max(bst.i, na.rm = TRUE)-min(bst.i, na.rm = TRUE))
		col.i <- round(col.i)
		cols  <- topo.colors(100)
		col.i <- cols[c(col.i, recursive = TRUE)]
		col.i[is.na(col.i)] <- "gray"
		if(i != ncol(B)){
			plot(phy, show.tip.label = FALSE, no.margin = TRUE)
		} else {
			plot(phy, cex = 0.2, no.margin = TRUE)
		}
		nodelabels.tim(bg = col.i, border = col.i, text = "")
		legend("bottomleft", legend = sprintf("%1.2f", seq(min(bst.i, na.rm = TRUE), max(bst.i, na.rm = TRUE), length = 10)), fill = topo.colors(10), bg = "white")
	}
	title(main = "Values of Bst", outer = TRUE, line = -1)
	if(!is.null(outfile)){
		dev.off()
	}
}

