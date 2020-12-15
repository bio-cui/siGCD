#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(dplyr))

ptm <- proc.time()
######################################################################################################################
#   7 parameters: 1.positive genes of pivot cell (separated by a comma)
#                 2.negative genes of pivot cell ("none" if there are no negative genes, separated by a comma)
#                 3.partner genes ("none" if use all genes as partner genes)
#                 4.gene expression data file (eg. "data/RNAseqLIHC.Rdata")
#                 5.clinical data file (eg. "data/ClinicalLIHC.Rdata")
#                 6.cofounding factors (eg. "Age" or c("Age","TNM staging:))
#                 7.project name
#       output:   1.cell_exp.csv
#                 2.cor_rprr.csv
#                 4.lable2/pic1/pic.txt
#                 5.lable2/pic2/download.txt
######################################################################################################################

##################### expression level
commands=commandArgs(trailingOnly=T)

source("funcs.R")

datalist=loaddata(commands[4],commands[5]) #4 expression data, 5 clinical
expma=datalist[[1]]
clinical=datalist[[2]]

NN = dim(expma)[1]  #number of genes
MM = dim(expma)[2]  #number of patients

cellss=cellcal(commands[1],commands[2],expma,MM)
file_name = commands[7]
gene_or_cell = "gene"

write.table (cellss,file=paste0(file_name,"cell_exp.csv"),row.names=F,col.names=F,quote=F)
##################### expression correlation
if(gene_or_cell == "gene" | gene_or_cell == "cell" | gene_or_cell =="DrugGene"){
	if(commands[3]=="none"){
		partnergs=rownames(expma)
		partnergsout=partnergs
		cor_rprr=matrix(0,NN,3)
		for(m in 1:NN){
			corres=cor.test(cellss,expma[m,],method="pearson")
			cor_rprr[m,1]=corres$estimate
			cor_rprr[m,2]=corres$p.value
			lm.reg=lm(cellss~expma[m,])
			cor_rprr[m,3]=summary(lm.reg)$r.squared
		}
	}else{
		partnergs=strsplit(commands[3],",")[[1]]
		partnergsout=partnergs
		geneintid = match(partnergs,rownames(expma))
		cor_rprr=matrix(0,length(geneintid),3)
		for(n in 1:length(geneintid)){
			corres=cor.test(cellss,expma[geneintid[n],],method="pearson")
			cor_rprr[n,1]=corres$estimate
			cor_rprr[n,2]=corres$p.value
			lm.reg=lm(cellss~expma[geneintid[n],])
			cor_rprr[n,3]=summary(lm.reg)$r.squared
		}
	}
}

ordss=order(abs(cor_rprr[,1]),decreasing = T)
newcor_rprr=cbind(partnergsout,cor_rprr)
colnames(newcor_rprr)=c("Genes","pearson's correlation","p.value","r.squared")

if(length(ordss)==1){
	cor_rprr = t(as.matrix(newcor_rprr[ordss,]))
}else{
	cor_rprr = newcor_rprr[ordss,]
}

write.csv(cor_rprr,row.names=F,file=paste0(file_name,"cor_rprr.csv"),quote=F)


##################### survival interaction
sur2=c("Survival Time","Death")
colnames(clinical)=gsub("\\."," ",colnames(clinical))
surv = Surv(as.numeric(clinical[,match(sur2[1],colnames(clinical))]),as.numeric(clinical[,match(sur2[2],colnames(clinical))]))

cellsh=rep(1,length(cellss))
cellsh[which(cellss>median(cellss))]=2

if(commands[6]!="none"){
	facts=strsplit(commands[6],",")[[1]]
	factid=match(facts,colnames(clinical))
	B = clinical[,factid, drop=F]
}

if(gene_or_cell == "gene" | gene_or_cell == "cell" | gene_or_cell =="DrugGene"){
	if(commands[3]=="none"){
		result_interaction=matrix(0,NN,2)
		rankss=apply(t(expma),2,order)
		halfs=round(MM/2)
		for(i in 1:NN){
			partner=rep(2,MM)
			partner[rankss[1:halfs,i]]=1
			Interaction=cellsh * partner
			B3=cofac_varsim(commands[6],B,cbind(cellsh,partner,Interaction))
			result_interaction[i,]=coxcalsim(surv,B3)
		}
	}else{
		result_interaction=matrix(0,length(geneintid),2)   #交互作用的z和p值
		for(n in 1:length(geneintid)){
			partner=rep(1,MM)
			partner[which(expma[geneintid[n],]>median(expma[geneintid[n],]))]=2
			Interaction=cellsh * partner
			B3=cofac_vars(commands[6],clinical,cbind(cellsh,partner,Interaction))
			result_interaction[n,]=coxcalsim(surv,B3)
		}
	}
}

ordss=order(abs(result_interaction[,1]),decreasing = T)
newresult_int=cbind(partnergsout,result_interaction)
colnames(newresult_int)=c("Genes","z.Score","p.value")

if(length(ordss)==1){
	result_interaction =t(as.matrix(newresult_int[ordss,]))
}else{
	result_interaction = newresult_int[ordss,]
}
write.csv(result_interaction,row.names=F,file=paste0(file_name,"result_interaction.csv"),quote=F)
first_gene2 <- partnergs[ordss[1]]
if(length(partnergs)==1){ #单药的情况
	xx=strsplit(commands[3],"-")[[1]]
	if(length(xx)>1){
		first_gene2 = xx[2]
	}
}

###########################################
cellss=as.numeric(read.table(paste0(file_name,"cell_exp.csv"))[[1]])

#################### 计算基因表达与细胞的含量对预后的交互作用
colnames(clinical)=gsub("\\."," ",colnames(clinical))
sur2=c("Survival Time","Death")
surv = Surv(as.numeric(clinical[,match(sur2[1],colnames(clinical))]),as.numeric(clinical[,match(sur2[2],colnames(clinical))]))

# 2 groups
if(gene_or_cell == "gene" | gene_or_cell == "cell" | gene_or_cell =="DrugGene"){
	cellss2=expma[match(first_gene2,rownames(expma)),]
}
######################################
diff=sum(abs(cellss2-cellss))
if(diff<0.1){
	cellss2=cellss2[sample(1:MM,MM)]
}
######################################
cellsh2=rep(1,MM)
cellsh2[which(cellss2>median(cellss2))]=2
#
B1 = cofac_vars(commands[6],clinical,cellsh)
coxout=coxcal(surv,B1)
resultcell_pz = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]
if(gene_or_cell=="cell" | gene_or_cell=="CellDrug"){
	tuli="Pivot Cell Low|Pivot Cell High"
}else if(gene_or_cell=="gene"){
	tuli="Pivot Gene Low|Pivot Gene High"
}else{
	tuli="Pivot Drug Low|Pivot Drug High"
}
if(!file.exists(paste0(file_name,"lable2")))	dir.create(path=paste0(file_name,"lable2"))
if(!file.exists(paste0(file_name,"lable2/pic1")))	dir.create(path=paste0(file_name,"lable2/pic1"))
survout(surv,cellsh,paste0(file_name,"lable2/pic1/"),tuli)
write.table(c("Total Cohort",tuli,resultcell_pz),file=paste0(file_name,"lable2/pic1/pic.txt"),row.names=F,quote=F)
res=cbind(cellss,clinical)
colnames(res)=c("Variable",colnames(clinical))
write.table(res,file=paste0(file_name,"lable2/pic1/download.txt"),sep='\t',quote=F)
############################################
idgh=which(cellsh2==2)
coxout=coxcal(surv[idgh,],as.data.frame(B1[idgh,]))
resultcell_pz = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]

if(!file.exists(paste0(file_name,"lable2")))	dir.create(path=paste0(file_name,"lable2"))
if(!file.exists(paste0(file_name,"lable2/pic2")))	dir.create(path=paste0(file_name,"lable2/pic2"))
survout(surv[idgh,],cellsh[idgh],paste0(file_name,"lable2/pic2/"),tuli)
write.table(c("Partner Gene High Subcohort",tuli,resultcell_pz),file=paste0(file_name,"lable2/pic2/pic.txt"),row.names=F,quote=F)
res=cbind(cellss[idgh],clinical[idgh,])
colnames(res)=c("Variable",colnames(clinical))
write.table(res,file=paste0(file_name,"lable2/pic2/download.txt"),,sep='\t',quote=F)
###########################################
idgl=which(cellsh2==1)
coxout=coxcal(surv[idgl,],as.data.frame(B1[idgl,]))
resultcell_pz = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]

if(!file.exists(paste0(file_name,"lable2")))	dir.create(path=paste0(file_name,"lable2"))
if(!file.exists(paste0(file_name,"lable2/pic3")))	dir.create(path=paste0(file_name,"lable2/pic3"))
survout(surv[idgl,],cellsh[idgl],paste0(file_name,"lable2/pic3/"),tuli)
write.table(c("Partner Gene Low Subcohort",tuli,resultcell_pz),file=paste0(file_name,"lable2/pic3/pic.txt"),row.names=F,quote=F)
res=cbind(cellss[idgl],clinical[idgl,])
colnames(res)=c("Variable",colnames(clinical))
write.table(res,file=paste0(file_name,"lable2/pic3/download.txt"),sep='\t',quote=F)

############################
B2 = cofac_vars(commands[6],clinical,cellsh2)
coxout=coxcal(surv,B2)
result_pg = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]

if(!file.exists(paste0(file_name,"lable3")))	dir.create(path=paste0(file_name,"lable3"))
if(!file.exists(paste0(file_name,"lable3/pic1")))	dir.create(path=paste0(file_name,"lable3/pic1"))
survout(surv,cellsh2,paste0(file_name,"lable3/pic1/"),"Partner Gene Low|Partner Gene High")
write.table(c("Total Cohort","Partner Gene Low|Partner Gene High",result_pg),file=paste0(file_name,"lable3/pic1/pic.txt"),row.names=F,quote=F)
res=cbind(cellss2,clinical)
colnames(res)=c("Variable",colnames(clinical))
write.table(res,file=paste0(file_name,"lable3/pic1/download.txt"),sep='\t',quote=F)
############################################
idgh=which(cellsh==2)
coxout=coxcal(surv[idgh,],as.data.frame(B2[idgh,]))
resultcell_pz = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]
if(gene_or_cell=="cell" | gene_or_cell=="CellDrug"){
	labless="Pivot Cell High Subcohort"
}else if(gene_or_cell=="gene"){
	labless="Pivot Gene High Subcohort"
}else{
	labless="Pivot Drug High Subcohort"
}
if(!file.exists(paste0(file_name,"lable3/pic2")))	dir.create(path=paste0(file_name,"lable3/pic2"))
survout(surv[idgh,],cellsh2[idgh],paste0(file_name,"lable3/pic2/"),"Partner Gene Low|Partner Gene High")
write.table(c(labless,"Partner Gene Low|Partner Gene High",resultcell_pz),file=paste0(file_name,"lable3/pic2/pic.txt"),row.names=F,quote=F)
res=cbind(cellss2[idgh],clinical[idgh,])
colnames(res)=c("Variable",colnames(clinical))
write.table(res,file=paste0(file_name,"lable3/pic2/download.txt"),,sep='\t',quote=F)

###########################################
idgl=which(cellsh==1)
coxout=coxcal(surv[idgl,],as.data.frame(B2[idgl,]))
resultcell_pz = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]

if(gene_or_cell=="cell" | gene_or_cell=="CellDrug"){
	labless="Pivot Cell Low Subcohort"
}else if(gene_or_cell=="gene"){
	labless="Pivot Gene Low Subcohort"
}else{
	labless="Pivot Drug Low Subcohort"
}
if(!file.exists(paste0(file_name,"lable3/pic3")))	dir.create(path=paste0(file_name,"lable3/pic3"))
survout(surv[idgl,],cellsh2[idgl],paste0(file_name,"lable3/pic3/"),"Partner Gene Low|Partner Gene High")
write.table(c(labless,"Partner Gene Low|Partner Gene High",resultcell_pz),file=paste0(file_name,"lable3/pic3/pic.txt"),row.names=F,quote=F)
res=cbind(cellss2[idgl],clinical[idgl,])
colnames(res)=c("Variable",colnames(clinical))
write.table(res,file=paste0(file_name,"lable3/pic3/download.txt"),sep='\t',quote=F)

##############################################
class4=rep(1,MM)
class4[intersect(which(cellsh2==2),which(cellsh==1))]=2
class4[intersect(which(cellsh2==1),which(cellsh==2))]=3
class4[intersect(which(cellsh2==2),which(cellsh==2))]=4
B2 = cofac_vars(commands[6],clinical,class4)
coxout=coxcal(surv,B2)
result_pg = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]
labels4="Pivot Cell Low & Partner Gene Low|Pivot Cell Low & Partner Gene High|Pivot Cell High & Partner Gene Low|Pivot Cell High & Partner Gene High"
if(gene_or_cell=="gene"){
	labels4="Pivot Gene Low & Partner Gene Low|Pivot Gene Low & Partner Gene High|Pivot Gene High & Partner Gene Low|Pivot Gene High & Partner Gene High"
}
if(gene_or_cell=="CellDrug"){
	labels4="Pivot Cell Low & Partner Drug Low|Pivot Cell Low & Partner Drug High|Pivot Cell High & Partner Drug Low|Pivot Cell High & Partner Drug High"
}
if(gene_or_cell=="DrugGene"){
	labels4="Pivot Drug Low & Partner Gene Low|Pivot Drug Low & Partner Gene High|Pivot Drug High & Partner Gene Low|Pivot Drug High & Partner Gene High"
}
if(gene_or_cell=="DrugDrug"){
	labels4="Pivot Drug Low & Partner Drug Low|Pivot Drug Low & Partner Drug High|Pivot Drug High & Partner Drug Low|Pivot Drug High & Partner Drug High"
}
if(!file.exists(paste0(file_name,"lable1")))	dir.create(path=paste0(file_name,"lable1"))
survout(surv,class4,paste0(file_name,"lable1/"),labels4)
write.table(c("Overall Survival of 4 Groups",labels4,result_pg),file=paste0(file_name,"lable1/pic.txt"),row.names=F,quote=F)
res=cbind(class4,clinical)
colnames(res)=c("Variable",colnames(clinical))
write.table(res,file=paste0(file_name,"lable1/download.txt"),sep='\t',quote=F)
####################################################################
if(gene_or_cell=="gene"){
	load("./data/ceres.Rdata")
	load("./data/demeter.Rdata")
	load("./data/slpair.Rdata")

	pivotg=commands[1]
	
	depceres_piv=ceres[match(pivotg,rownames(ceres)),]
	if(is.na(depceres_piv[1])==1){
		depceres_piv=rep(0,563)
		names(depceres_piv)=colnames(ceres)
	}
	depdemeter_piv=demeter[match(pivotg,rownames(demeter)),]
	if(is.na(depdemeter_piv[1])==1){
		depdemeter_piv=rep(0,501)
		names(depdemeter_piv)=colnames(demeter)
	}
	write.table(depceres_piv,file=paste0(file_name,"/depceres_piv.txt"),quote=F,col.names=F,sep='\t')  #改目录
	write.table(depdemeter_piv,file=paste0(file_name,"/depdemeter_piv.txt"),quote=F,col.names=F,sep='\t') #改目录
	depceres_par=ceres[match(partnergs,rownames(ceres)),]
	depdemeter_par=demeter[match(partnergs,rownames(demeter)),]
	cor_rprrd=matrix(0,length(partnergs),3)
	cor_rprrc=matrix(0,length(partnergs),3)
	for(m in 1:length(partnergs)){
		
		sec=depdemeter_par
		if(length(partnergs)>1){
			sec=depdemeter_par[m,]
		}
		if(is.na(sec)==0){
			corres=cor.test(depdemeter_piv,sec,method="pearson")
			cor_rprrd[m,1]=corres$estimate  #相关系数
			cor_rprrd[m,2]=corres$p.value   #p值
			lm.reg=lm(depdemeter_piv~sec)
			cor_rprrd[m,3]=summary(lm.reg)$r.squared  #回归R方值
		}
		
		sec=depceres_par
		if(length(partnergs)>1){
			sec=depceres_par[m,]
		}
		if(is.na(sec)==0){
			corres=cor.test(depceres_piv,sec,method="pearson")
			cor_rprrc[m,1]=corres$estimate  #相关系数
			cor_rprrc[m,2]=corres$p.value   #p值
			lm.reg=lm(depceres_piv~sec)
			cor_rprrc[m,3]=summary(lm.reg)$r.squared  #回归R方值
		}
	}

	ordss=order(abs(cor_rprrc[,1]),decreasing = T)
	newcor_rprrc=cbind(partnergs,cor_rprrc[,1],cor_rprrd[,1])
	#newcor_rprrd=cbind(partnergs,cor_rprrd)
	colnames(newcor_rprrc)=c("Genes","CRISPER co-dependence","RNAi co-dependence")
	#colnames(newcor_rprrd)=c("Genes","pearson's correlation","p.value","r.squared")
	
	if(length(ordss)==1){
		cor_rprrc = t(as.matrix(newcor_rprrc[ordss,]))
	}else{
		cor_rprrc = newcor_rprrc[ordss,]
	}
	write.csv(cor_rprrc,row.names=F,file=paste0(file_name,"cor_rprrc.csv"),quote=F)
	first_gene <- cor_rprrc[1,1]
	system(paste("Rscript scatter_plotdep.R",first_gene,file_name,sep =" "))
	
	#######################################################
	out=cbind(rep(NA,length(partnergs)),rep(NA,length(partnergs)),rep(NA,length(partnergs)))
	idp1=which(slpair[,1]==pivotg)
	idp2=which(slpair[,2]==pivotg)
	for(m in 1:length(partnergs)){
		idp12=which(slpair[,2]==partnergs[m])
		idp21=which(slpair[,1]==partnergs[m])
		xx=na.omit(union(intersect(idp1,idp12),intersect(idp2,idp21)))
		if(length(xx)>0){
			out[m,]=slpair[xx[1],3:5]
		}
	}
	outres=cbind(rep(pivotg,length(partnergs)),partnergs,out)
	colnames(outres)=c("Pivort Gene","Partner Gene","PubmedID","EvidenceType","SL type")
	
	if(length(partnergs)>1){
		ordss=order(outres[,3])
		outres=outres[ordss,]
	}
	write.table(outres,file=paste0(file_name,"/synthetic_lethality.txt"),quote=F,row.names=F,sep='\t')

}

proc.time() - ptm