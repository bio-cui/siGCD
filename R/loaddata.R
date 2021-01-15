loaddata<-function(expfile,clifile){
	datalist=list()
	load(paste('data/',expfile,sep=''))
	datalist[[1]]=expma
	load(paste('data/',clifile,sep=''))
	datalist[[2]]=clinical
	datalist
}
