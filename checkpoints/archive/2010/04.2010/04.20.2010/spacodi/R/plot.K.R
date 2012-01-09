plot.K <-
function(phy, traits, outfile = NULL){
	K=data.frame(K.all.nodes(phy,traits,return.all=TRUE)[,1,])
	names(K)=names(traits)
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
			plot(phy, show.tip.label = FALSE, no.margin = TRUE)
		} else {
			plot(phy, cex = 0.2, no.margin = TRUE)
		}
		nodelabels.tim(bg = col.i, border = col.i, text = "")
		legend("bottomleft", legend = sprintf("%1.2f", seq(min(trt.i, na.rm = TRUE), max(trt.i, na.rm = TRUE), length = 10)), fill = topo.colors(10), bg = "white")
		title(main = names(K)[i], line = -2)
	}
	title(main = "Values of Blomberg's K", outer = TRUE, line = -1)
	if(!is.null(outfile)){
		dev.off()
	}
}

