############load DATA
loaddata<-function(expfile,clifile){
	datalist=list()
	load(paste('data/',expfile,sep=''))
	datalist[[1]]=expma
	load(paste('data/',clifile,sep=''))
	datalist[[2]]=clinical
	datalist
}

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
cofac_varsim<-function(facts,B,vars){
	if(length(facts)==0){
		B1 = as.data.frame(vars)
	} else {
		B1=as.data.frame(cbind(B,vars))
	}
	B1
}

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

coxcal<-function(surv,datas){
    errflag = F
	coxph.fit = tryCatch(coxph(surv~., data=datas),
		                 error = function(e) errflag <<- T)
	if(!errflag){
		reg.summary = summary(coxph.fit)$coef
	}
	reg.summary
}
coxcalsim<-function(surv,datas){
    errflag = F
	coxph.fit = tryCatch(coxph(surv~., data=datas),
		                 error = function(e) errflag <<- T)
	if(!errflag){
		reg.summary = summary(coxph.fit)$coef
	}
	reg.summary["Interaction", c("z", "Pr(>|z|)")]
}

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

