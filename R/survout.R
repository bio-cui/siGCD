survout<-function(surv,exprs_data_cutoff,path_output,tulis){
	kaplan_meier <- survfit(surv ~ exprs_data_cutoff,conf.type = "log")
	xxx=summary(kaplan_meier)$table
	outinfo=data.frame(xxx[,c(1,4,7)])
    colnames(outinfo)=c("Number of Cases, Total","Number of Cases, Deceased","Median Survival Time")
	#rownames(outinfo)=strsplit(tulis,"\\|")[[1]]
	outinfo$Name = strsplit(tulis,"\\|")[[1]]
    outinfo<-select(outinfo,4,1,2,3)
	write.table(outinfo,file=paste(path_output,"outinfo.txt",sep=''),sep='\t',quote=F,row.names=F)
	NN=length(kaplan_meier$strata)
	timepoints=as.numeric(c(1,kaplan_meier$strata))
	kmall=cbind(kaplan_meier$surv,kaplan_meier$time,kaplan_meier$n.event,kaplan_meier$n.censor)
	for(i in 1:NN){
		low_matrix <- matrix(nrow=kaplan_meier$strata[i],ncol=4)
		if(kaplan_meier$strata[i]>1){
			idl=sum(timepoints[1:i])
			idh=sum(timepoints[1:(i+1)])-1
			low_matrix <- kmall[idl:idh,]
			colnames(low_matrix) <- c("surv", "time", "event", "censor")
			write.table(low_matrix,file=paste(path_output,"survival_",as.character(i),".txt",sep = ""),sep='\t',quote=F,row.names = FALSE)
		} else {
			nopat=c(0,0,0,0)
			names(nopat)=c("surv", "time", "event", "censor")
			write.table(t(nopat),file=paste(path_output,"survival_",as.character(i),".txt",sep = ""),sep='\t',quote=F,row.names=F)
		}
	}
}
