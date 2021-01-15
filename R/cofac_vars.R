cofac_vars<-function(facts,clinical,vars){
	if(length(facts)==0){
		B1 = vars
	} else {
		facts=strsplit(facts,",")[[1]]
		factid=match(facts,colnames(clinical))
		B = clinical[,factid, drop=F]
		B1=cbind(B,vars)
	}
	B1=as.data.frame(B1)
	B1
}
