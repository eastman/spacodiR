write.spacodi.data <-
function(data, outfile){
	if(file.exists(outfile)){
		warning("Overwrote existing outfile.")
		unlink(outfile)
	}
	names=names(data)
	for(n in 1:length(names)){cat(c("\t",names[n]),file=outfile,append=TRUE,sep="")}
	cat("\n",file=outfile,append=TRUE,sep="")
	write.table(data,outfile,quote=FALSE,col=FALSE,append=TRUE,sep="\t")
}

