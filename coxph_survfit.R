#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(dplyr))

##################### 
commands=commandArgs(trailingOnly=T)
file_name = commands[7]
gene_or_cell = "cell"
if (length(commands) == 8){
    gene_or_cell = commands[8]}
# file_name_sub = sub(pattern="/",replacement="", file_name)
# if(!file.exists(file_name_sub))	dir.create(path=file_name_sub)
source("funcs.R")

datalist=loaddata(commands[4],commands[5])
expma=datalist[[1]]
clinical=datalist[[2]]

NN = dim(expma)[1]
MM = dim(expma)[2]

cellss=as.numeric(read.table(paste0(file_name,"cell_exp.csv"))[[1]])

##################### survival interaction
colnames(clinical)=gsub("\\."," ",colnames(clinical))
sur2=c("Survival Time","Death")
surv = Surv(as.numeric(clinical[,match(sur2[1],colnames(clinical))]),as.numeric(clinical[,match(sur2[2],colnames(clinical))]))

#2 groups
cellsh=rep(1,length(cellss))
cellsh[which(cellss>median(cellss))]=2
if(gene_or_cell == "gene" | gene_or_cell == "cell" | gene_or_cell =="DrugGene"){
	cellss2=expma[match(commands[3],rownames(expma)),]
}else if(commands[3]=="multi"){
	drugbank=read.table("drugbank.txt",sep='\t')
	targets=paste(as.character(drugbank[which(as.character(drugbank[,1])==commands[3]),2]),collapse=',')
	cellss2=cellcal(targets,"none",expma,MM)
}else{
	cellss2=cellcal(commands[3],"none",expma,MM)
}
######################################
diff=sum(abs(cellss2-cellss))
if(diff<0.1){
	cellss2=cellss2[sample(1:MM,MM)]
}
######################################

cellsh2=rep(1,MM)
cellsh2[which(cellss2>median(cellss2))]=2

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
resultcell_pz = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]

if(!file.exists(paste0(file_name,"lable3")))	dir.create(path=paste0(file_name,"lable3"))
if(!file.exists(paste0(file_name,"lable3/pic1")))	dir.create(path=paste0(file_name,"lable3/pic1"))
survout(surv,cellsh2,paste0(file_name,"lable3/pic1/"),"Partner Gene Low|Partner Gene High")
write.table(c("Total Cohort","Partner Gene Low|Partner Gene High",resultcell_pz),file=paste0(file_name,"lable3/pic1/pic.txt"),row.names=F,quote=F)
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
resultcell_pz = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]
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
write.table(c("Overall Survival of 4 Groups",labels4,resultcell_pz),file=paste0(file_name,"lable1/pic.txt"),row.names=F,quote=F)
res=cbind(class4,clinical)
colnames(res)=c("Variable",colnames(clinical))
write.table(res,file=paste0(file_name,"lable1/download.txt"),sep='\t',quote=F)
write.table(c(),file=paste0(file_name,"end.txt"),row.names=F,col.names=F)
