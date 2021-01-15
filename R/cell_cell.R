#!/usr/bin/env Rscript
cell_cell <- function(posg1,negg1=c(),posg2,negg2=c(),expfile,clifile,confoundf=c()) {
#################################################################################################################
#    input:  posg1: positive genes of pivot cell
#            negg1: negative genes of pivot cell  ("none" if there are no negative genes)
#            posg2: positive genes of partner cell
#            negg2: negative genes of partner cell  ("none" if there are no negative genes)
#            expfile: gene expression data file (eg. "data/RNAseqLIHC.Rdata")
#            clifile: clinical data file (eg. "data/ClinicalLIHC.Rdata")
#            confoundf: cofounding factors (eg. "Age" or c("Age","TNM staging"))
#    output: a list which contains "survival interaction", "piv_survival total","piv_survival high sub","piv_survival low sub"
#            "part_survival total","part_survival high sub","part_survival high sub","survival 4 group"
################################################################################################################

	output=list()
	suppressPackageStartupMessages(library(survival))
	suppressPackageStartupMessages(library(dplyr))
	source("funcs.R")
	datalist=loaddata(expfile,clifile)
	expma=datalist[[1]]
	clinical=datalist[[2]]
	NN = dim(expma)[1]
	MM = dim(expma)[2]
	cellss1=cellcal(posg1,negg1,expma,MM)
	cellss2=cellcal(posg2,negg2,expma,MM)

	##################### expression level
	colnames(clinical)=gsub("\\."," ",colnames(clinical))
	sur2=c("Survival Time","Death")
	surv = Surv(as.numeric(clinical[,match(sur2[1],colnames(clinical))]),as.numeric(clinical[,match(sur2[2],colnames(clinical))]))

	cellsh1=rep(1,length(cellss1))
	cellsh1[which(cellss1>median(cellss1))]=2
	cellsh2=rep(1,length(cellss2))
	cellsh2[which(cellss2>median(cellss2))]=2

	Interaction=cellsh1 * cellsh2
	B3=cofac_vars(confoundf,clinical,cbind(cellsh1,cellsh2,Interaction))
	coxout=coxcal(surv,B3)
	result_interaction = coxout["Interaction", c("z", "Pr(>|z|)")]
	names(result_interaction)=c("z.Score","p.value")
	output[["surv_interaction"]]=result_interaction

	#####################################
	B1 = cofac_vars(confoundf,clinical,cellsh1)
	coxout=coxcal(surv,B1)
	output[["piv_surv_total"]] = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")] #pivot cell

	############################################
	idgh=which(cellsh2==2)
	coxout=coxcal(surv[idgh,],as.data.frame(B1[idgh,]))
	output[["piv_survival high sub"]] = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]

	###########################################
	idgl=which(cellsh2==1)
	coxout=coxcal(surv[idgl,],as.data.frame(B1[idgl,]))
	output[["piv_surv_low_sub"]] = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]

	############################################
	B2 = cofac_vars(confoundf,clinical,cellsh2)
	coxout=coxcal(surv,B2)
	output[["part_surv_total"]] = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]

	############################################
	idgh=which(cellsh1==2)
	coxout=coxcal(surv[idgh,],as.data.frame(B2[idgh,]))
	output[["part_surv_high_sub"]] = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]

	###########################################
	idgl=which(cellsh1==1)
	coxout=coxcal(surv[idgl,],as.data.frame(B2[idgl,]))
	output[["part_surv_low_sub"]] = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]

	##############################################
	class4=rep(1,MM)
	class4[intersect(which(cellsh2==2),which(cellsh1==1))]=2
	class4[intersect(which(cellsh2==1),which(cellsh1==2))]=3
	class4[intersect(which(cellsh2==2),which(cellsh1==2))]=4
	B2 = cofac_vars(confoundf,clinical,class4)
	coxout=coxcal(surv,B2)
	output[["surv_4_group"]] = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]

	output
}
