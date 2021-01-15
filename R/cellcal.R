cellcal<-function(posg,negg=c(),expma,MM){
	cellpos=cellneg=rep(0,MM)
	if(length(posg)>0){
		cellgeneid = match(posg,rownames(expma))
		cellpos = expma[cellgeneid[1],]
		if(length(posg)>1){
			cellpos = colSums(expma[cellgeneid,])
		}
	}
	if(length(negg)>0){
		cellgeneid = match(negg,rownames(expma))
		cellneg = expma[cellgeneid[1],]
		if(length(negg)>1){
			cellneg = colSums(expma[cellgeneid,])
		}
	}
	cellss=(cellpos-cellneg)/(length(posg)+length(negg))
	cellss
}
