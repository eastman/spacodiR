resamp.3x <-
function(obj, level=0.1) {
	swaps=round(level*ncol(obj)*nrow(obj))
	orig=obj
	for(swap in 1:swaps){
		rcol=sample(1:ncol(obj),2)
		rspp=sample(1:nrow(obj),2)
		obj[rspp[1],rcol[1]]=orig[rspp[1],rcol[2]]
		obj[rspp[1],rcol[2]]=orig[rspp[1],rcol[1]]
		obj[rspp[2],rcol[1]]=orig[rspp[2],rcol[2]]
		obj[rspp[2],rcol[2]]=orig[rspp[2],rcol[1]]
		orig=obj
	}
	if(sum(obj)!=sum(orig))warning("A poor result is likely.")
	return(obj)
}

