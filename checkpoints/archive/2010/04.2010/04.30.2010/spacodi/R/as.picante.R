as.picante <-
function(data, outfile=NULL){
	if(ncol(data)==3) {
		message("Formatting assumed to be that used with phylocom")
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
			cur.array[which(is.na(cur.array))]=0
			out[,plot]=as.numeric(cur.array)
		}
		out=as.data.frame(t(out))
		names(out)=species
		row.names(out)=names(dd)
		if(!is.null(outfile)){write.spacodi.data(out,outfile)}
		return(out)
	} else {
		message("Formatting assumed to be that used with spacodi")
		if(class(data)!="data.frame"){
			if(class(data)=="matrix") data=as.data.frame(data) else if(class(data)=="vector") data=as.data.frame(data) else stop("Function cannot handle data that are not in data.frame format")
		}
		if(is.null(row.names(data)))row.names(data)=paste("plot",seq(1:nrow(data)))
		if(is.null(names(data)))names(data)=paste("species",seq(1:nrow(data)))
	
		out=as.matrix(t(data))
		if(!is.null(outfile)){write.table(out,outfile,quote=FALSE)}
		return(out)
	}
}
