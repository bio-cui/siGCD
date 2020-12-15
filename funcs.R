############load DATA
loaddata<-function(expfile,clifile){
	datalist=list()
	load(expfile)
	datalist[[1]]=expma
	load(clifile)
	datalist[[2]]=clinical
	datalist
}

cofac_vars<-function(facts,clinical,vars){
	if(facts=="none"){
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
	if(facts=="none"){
		B1 = as.data.frame(vars)
	} else {
		B1=as.data.frame)cbind(B,vars))
	}
	B1
}

survout<-function(surv,exprs_data_cutoff,path_output,tulis){
	kaplan_meier <- survfit(surv ~ exprs_data_cutoff,conf.type = "log") #conf.type="log-log"?
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
    # print(errflag)
	if(!errflag){
		reg.summary = summary(coxph.fit)$coef
	}
	reg.summary
}
coxcalsim<-function(surv,datas){
    errflag = F
	coxph.fit = tryCatch(coxph(surv~., data=datas),
		                 error = function(e) errflag <<- T)
    # print(errflag)
	if(!errflag){
		reg.summary = summary(coxph.fit)$coef
	}
	reg.summary["Interaction", c("z", "Pr(>|z|)")]
}

cellcal<-function(posg,negg,expma,MM){
	gene_pos=gene_neg=c()
	cellpos=cellneg=rep(0,MM)
	if(posg!="none"){
		gene_pos=strsplit(posg,",")[[1]]
		cellgeneid = match(gene_pos,rownames(expma)) # 参数1：细胞阳性基因列表，gene symbol
		cellpos = expma[cellgeneid[1],]
		if(length(gene_pos)>1){
			cellpos = colSums(expma[cellgeneid,])   #应该是MM长度的向量，每个病人该类细胞的含量
		}
	}
	if(negg!="none"){
		gene_neg=strsplit(negg,",")[[1]]
		cellgeneid = match(gene_neg,rownames(expma)) # 参数2：细胞阴性基因列表，gene symbol
		cellneg = expma[cellgeneid[1],]
		if(length(gene_neg)>1){
			cellneg = colSums(expma[cellgeneid,])
		}
	}
	cellss=(cellpos-cellneg)/(length(gene_pos)+length(gene_neg))   #MM长度的向量，每个病人该类细胞的含量
	cellss
}

