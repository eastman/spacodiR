as.spacodi <-
function(data, outfile=NULL){
	if(ncol(data)!=3)stop(cat("Format does not appear to be in 'triplet'.\n\tColumns should be plot, sample(s), species.\n\n"))
	if(class(data)!="data.frame"){
		if(class(data)=="matrix") data=as.data.frame(data) else if(class(data)=="vector") data=as.data.frame(data) else stop("Function cannot handle data that are not in data.frame format")
	}
	species=as.character(unique(data[,3]))
	dd=split(data,data[,1])
	out=array(dim=c(length(species)->spp,length(dd)->plots))
	for(plot in 1:plots){
		cur.array=array(dim=c(spp, 1))
		ind=as.numeric(as.vector(dd[[plot]][,2]))
		dd[[plot]]->cur.plot
		names(ind)=as.character(cur.plot[,3])
		for(r in 1:nrow(cur.plot)){
			cur.array[which(names(ind[r])==species)]=ind[r]
		}
		out[,plot]=as.numeric(cur.array)
	}
	out=as.data.frame(out)
	names(out)=names(dd)
	row.names(out)=species
	if(!is.null(outfile)){write.spacodi.data(out,outfile)}
	return(out)
}

