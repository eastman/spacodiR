resamp.3t <-
function(obj) {
	for(ss in 1:nrow(obj)){
		torus=rep(1:ncol(obj),2)
		r=sample(1:(ncol(obj)-1),1)
		t.array=array(dim=c(ncol(obj),2))
		t.array[,1]=1:ncol(obj)
		for(o in 1:ncol(obj)){
			tt=torus[r+o]
			t.array[tt,2]=obj[ss,o]
		}
		obj[ss,]=t.array[,2]
	}
	return(obj)
}

