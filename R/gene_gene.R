#!/usr/bin/env Rscript
gene_gene <- function(pivg,partg=c(),expfile,clifile,confoundf=c()) {
#################################################################################################################
#    input:  pivg: pivot gene
#            partg: partner genes ("none" if use all genes as partner genes)
#            expfile: gene expression data file (eg. "data/RNAseqLIHC.Rdata")
#            clifile: clinical data file (eg. "data/ClinicalLIHC.Rdata")
#            confoundf: cofounding factors (eg. "Age" or c("Age","TNM staging"))
#    output: a list which contains "survival interaction", "piv_survival total","piv_survival high sub","piv_survival low sub"
#            "part_survival total","part_survival high sub","part_survival high sub","survival 4 group"
################################################################################################################
	output=list()
	suppressPackageStartupMessages(library(survival))
	suppressPackageStartupMessages(library(dplyr))
	#source("funcs.R")

	datalist=loaddata(expfile,clifile) #4 expression data, 5 clinical
	expma=datalist[[1]]
	clinical=datalist[[2]]

	NN = dim(expma)[1]  #number of genes
	MM = dim(expma)[2]  #number of patients

	cellss=expma[match(pivg,rownames(expma)),]

	##################### expression correlation
	if(length(partg)==0){
		partg=rownames(expma)
	}
	geneintid = match(partg,rownames(expma))

	##################### survival interaction
	sur2=c("Survival Time","Death")
	colnames(clinical)=gsub("\\."," ",colnames(clinical))
	surv = Surv(as.numeric(clinical[,match(sur2[1],colnames(clinical))]),as.numeric(clinical[,match(sur2[2],colnames(clinical))]))

	cellsh=rep(1,length(cellss))
	cellsh[which(cellss>median(cellss))]=2

	if(length(confoundf)>0){
		factid=match(facts,colnames(clinical))
		B = clinical[,factid, drop=F]
	}

	if(length(partg)==0){
		result_interaction=matrix(0,NN,2)
		rankss=apply(t(expma),2,order)
		halfs=round(MM/2)
		for(i in 1:NN){
			partner=rep(2,MM)
			partner[rankss[1:halfs,i]]=1
			Interaction=cellsh * partner
			B3=cofac_varsim(confoundf,B,cbind(cellsh,partner,Interaction))
			result_interaction[i,]=coxcalsim(surv,B3)
		}
	}else{
		result_interaction=matrix(0,length(geneintid),2)
		for(n in 1:length(geneintid)){
			partner=rep(1,MM)
			partner[which(expma[geneintid[n],]>median(expma[geneintid[n],]))]=2
			Interaction=cellsh * partner
			B3=cofac_vars(confoundf,clinical,cbind(cellsh,partner,Interaction))
			result_interaction[n,]=coxcalsim(surv,B3)
		}
	}
	output[["surv_interaction"]]=result_interaction

	respz1=respz2=respz3=respz4=respz5=respz6=respz7=matrix(0,length(geneintid),3)
	for(m in 1:length(geneintid)){
		cellss2=expma[geneintid[m],]
		cellsh2=rep(1,MM)
		cellsh2[which(cellss2>median(cellss2))]=2

		B1 = cofac_vars(confoundf,clinical,cellsh)
		coxout=coxcal(surv,B1)
		respz1[m,] = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]
		############################################
		idgh=which(cellsh2==2)
		coxout=coxcal(surv[idgh,],as.data.frame(B1[idgh,]))
		respz2[m,] = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]

		###########################################
		idgl=which(cellsh2==1)
		coxout=coxcal(surv[idgl,],as.data.frame(B1[idgl,]))
		respz3[m,] = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]

		############################
		B2 = cofac_vars(confoundf,clinical,cellsh2)
		coxout=coxcal(surv,B2)
		respz4[m,] = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]

		############################################
		idgh=which(cellsh==2)
		coxout=coxcal(surv[idgh,],as.data.frame(B2[idgh,]))
		respz5[m,] = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]

		###########################################
		idgl=which(cellsh==1)
		coxout=coxcal(surv[idgl,],as.data.frame(B2[idgl,]))
		respz6[m,] = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]

		##############################################
		class4=rep(1,MM)
		class4[intersect(which(cellsh2==2),which(cellsh==1))]=2
		class4[intersect(which(cellsh2==1),which(cellsh==2))]=3
		class4[intersect(which(cellsh2==2),which(cellsh==2))]=4
		B2 = cofac_vars(confoundf,clinical,class4)
		coxout=coxcal(surv,B2)
		respz7[m,] = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]
		colnames(respz1)=colnames(respz2)=colnames(respz3)=colnames(respz4)=colnames(respz5)=colnames(respz6)=colnames(respz7)=c("exp(coef)","z", "Pr(>|z|)")
	}
	output[["piv_surv_total"]]=respz1
	output[["piv_surv_high_sub"]]=respz2
	output[["piv_surv_low_sub"]]=respz3
	output[["part_surv_total"]]=respz4
	output[["part_surv_high_sub"]]=respz5
	output[["part_surv_low_sub"]]=respz6
	output[["surv_4_group"]]=respz7
	output
}
