nodelabels.tim <-
function (text, node, adj = c(0.5, 0.5), frame = "rect", pch = NULL, thermo = NULL, pie = NULL, piecol = NULL, col = "black", bg = "lightblue", border = "black", ...) {
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    if (missing(node)) 
	node <- (lastPP$Ntip + 1):length(lastPP$xx)
    XX <- lastPP$xx[node]
    YY <- lastPP$yy[node]
    BOTHlabels.tim(text, node, XX, YY, adj, frame, pch, thermo, pie, piecol, col, bg, border, ...)
}

