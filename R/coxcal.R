coxcal<-function(surv,datas){
    errflag = F
	coxph.fit = tryCatch(coxph(surv~., data=datas),
		                 error = function(e) errflag <<- T)
	if(!errflag){
		reg.summary = summary(coxph.fit)$coef
	}
	reg.summary
}
