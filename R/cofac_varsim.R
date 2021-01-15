cofac_varsim<-function(facts,B,vars){
	if(length(facts)==0){
		B1 = as.data.frame(vars)
	} else {
		B1=as.data.frame(cbind(B,vars))
	}
	B1
}
