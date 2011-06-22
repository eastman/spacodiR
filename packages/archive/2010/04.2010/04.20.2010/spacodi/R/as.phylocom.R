as.phylocom <-
function(data, outfile=NULL){
	if(class(data)!="data.frame"){
		if(class(data)=="matrix") data=as.data.frame(data) else if(class(data)=="vector") data=as.data.frame(data) else stop("Function cannot handle data that are not in data.frame format")
	}
	if(length(unique(row.names(data)))!=nrow(data))warning("Data do not appear to be in proper format.  Results may be nonsensical.")
	plots=ncol(data)
	spp=nrow(data)
	out=array(dim=c(spp*plots,3))
	for(plot in 1:plots){
		start=((plot-1)*spp)+1
		end=plot*spp
		out[start:end,]=cbind(names(data)[plot], data[,plot], row.names(data))
	}
	out=as.data.frame(out)
	names(out)=c("plot","samples","species")
	if(!is.null(outfile)){write.table(out,outfile,quote=FALSE,row=FALSE,col=FALSE,sep="\t")}
	return(out)
}

