#!/usr/bin/env Rscript
drug_drug <- function(drug1=c(),tar1=c(),drug2=c(),tar2=c(),expfile,clifile,confoundf) {
######################################################################################################################
#   *** make sure that the drugbank.txt file is in the work directory
#   input:   drug1: pivot drug name
#            tar1: targets of pivot drug ("none" if use targets of drug in DrugBank database)
#            drug2: partner drug ("none" if use all drugs in drugbank)
#            tar2: targets of partner drug ("none" if use targets of drug in DrugBank database)
#            expfile: gene expression data file (eg. "data/RNAseqLIHC.Rdata")
#            clifile: clinical data file (eg. "data/ClinicalLIHC.Rdata")
#            confoundf: cofounding factors (eg. "Age" or c("Age","TNM staging"))
#    output: a list which contains "survival interaction", "piv_survival total","piv_survival high sub","piv_survival low sub"
#            "part_survival total","part_survival high sub","part_survival high sub","survival 4 group"
######################################################################################################################
	output=list()
	suppressPackageStartupMessages(library(survival))
	suppressPackageStartupMessages(library(dplyr))
	output=list()
	#source("funcs.R")

	datalist=loaddata(expfile,clifile)
	expma=datalist[[1]]
	clinical=datalist[[2]]

	NN = dim(expma)[1]  #number of genes
	MM = dim(expma)[2]  #number of patients
	drugbank=load("data/drugbank.RData")
	idx=which(as.character(drugbank[,1])==drug1)
	if(length(idx)>0){
		targets=as.character(drugbank[idx,2])
	}
	if(length(tar1)>0){
		targets=drugtar
	}
	cellss=cellcal(targets,c(),expma,MM)
	cellsh=rep(1,length(cellss))
	cellsh[which(cellss>median(cellss))]=2

	###############################################################
	sur2=c("Survival Time","Death")
	colnames(clinical)=gsub("\\."," ",colnames(clinical))
	surv = Surv(as.numeric(clinical[,match(sur2[1],colnames(clinical))]),as.numeric(clinical[,match(sur2[2],colnames(clinical))]))

	if(length(confoundf)>0){
		factid=match(confoundf,colnames(clinical))
		B = clinical[,factid, drop=F]
	}

	################################################################
	tar2list=list()
	if(length(drug2)==0){
		partnergs=unique(as.character(drugbank[,1]))
		for(n in 1:length(partnergs)){
			tar2list[[n]] = as.character(drugbank[which(as.character(drugbank[,1])==partnergs[n]),2])
		}
	}else{
		idx=which(as.character(drugbank[,1])==drug2)
		if(length(idx)>0){
			tar2list[[1]]=as.character(drugbank[idx,2])
		}
		if(length(tar2)>0){
			tar2list[[1]]=tar2
		}
	}

	respz1=respz2=respz3=respz4=respz5=respz6=respz7=matrix(0,length(tar2list),3)
	result_interaction=matrix(0,length(tar2list),2)
	for(m in 1:length(tar2list)){
		drugv=cellcal(tar2list[[m]],c(),expma,MM)
		partner=rep(1,length(cellss))
		partner[which(drugv>median(drugv))]=2
		Interaction=cellsh * partner
		B3=cofac_vars(confoundf,clinical,cbind(cellsh,partner,Interaction))
		result_interaction[m,]=coxcalsim(surv,B3)

		B1 = cofac_vars(confoundf,clinical,cellsh)
		coxout=coxcal(surv,B1)
		respz1[m,] = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]

		idgh=which(partner==2)
		coxout=coxcal(surv[idgh,],as.data.frame(B1[idgh,]))
		respz2[m,] = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]

		idgl=which(partner==1)
		coxout=coxcal(surv[idgl,],as.data.frame(B1[idgl,]))
		respz3[m,] = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]

		B2 = cofac_vars(confoundf,clinical,partner)
		coxout=coxcal(surv,B2)
		respz4[m,] = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]

		idgh=which(cellsh==2)
		coxout=coxcal(surv[idgh,],as.data.frame(B2[idgh,]))
		respz5[m,] = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]

		idgl=which(cellsh==1)
		coxout=coxcal(surv[idgl,],as.data.frame(B2[idgl,]))
		respz6[m,] = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]

		class4=rep(1,MM)
		class4[intersect(which(partner==2),which(cellsh==1))]=2
		class4[intersect(which(partner==1),which(cellsh==2))]=3
		class4[intersect(which(partner==2),which(cellsh==2))]=4
		B2 = cofac_vars(confoundf,clinical,class4)
		coxout=coxcal(surv,B2)
		respz7[m,] = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]
		colnames(respz1)=colnames(respz2)=colnames(respz3)=colnames(respz4)=colnames(respz5)=colnames(respz6)=colnames(respz7)=c("exp(coef)","z", "Pr(>|z|)")
	}
	output[["surv_interaction"]]=result_interaction
	output[["piv_surv_total"]]=respz1
	output[["piv_surv_high_sub"]]=respz2
	output[["piv_surv_low_sub"]]=respz3
	output[["part_surv_total"]]=respz4
	output[["part_surv_high_sub"]]=respz5
	output[["part_surv_high_sub"]]=respz6
	output[["surv_4_group"]]=respz7
	output
}
