#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(dplyr))
######################################################################################################################
#  8 parameters:  1.positive genes of pivot cell (separated by a comma)
#                 2.negative genes of pivot cell (separated by a comma)
#                 3.positive genes of partner cell (separated by a comma)
#                 4.negative genes of partner cell (separated by a comma)
#                 5.gene expression data file (eg. "data/RNAseqLIHC.Rdata")
#                 6.clinical data file (eg. "data/ClinicalLIHC.Rdata")
#                 7.cofounding factors (eg. "Age" or c("Age","TNM staging"))
#                 8.project name

######################################################################################################################

##################### expression level
commands=commandArgs(trailingOnly=T)
file_name = commands[8]
# file_name_sub = sub(pattern="/",replacement="", file_name)
# if(!file.exists(file_name_sub)) dir.create(path=file_name_sub)

source("funcs.R")
datalist=loaddata(commands[5],commands[6])
expma=datalist[[1]]
clinical=datalist[[2]]
NN = dim(expma)[1]
MM = dim(expma)[2]

cellss1=cellcal(commands[1],commands[2],expma,MM)
cellss2=cellcal(commands[3],commands[4],expma,MM)

##################### expression correlation
cor_rprr=c()
corres=cor.test(cellss1,cellss2,method="pearson")
cor_rprr[1]=corres$estimate 
cor_rprr[2]=corres$p.value
lm.reg=lm(cellss1~cellss2)
cor_rprr[3]=summary(lm.reg)$r.squared

write.table(t(cor_rprr),row.names=F,col.names=F,file=paste0(commands[8],"cell_corr.csv"),quote=F)
write.table(t(cellss1),row.names=F,col.names=F,file=paste0(commands[8],"cell_corr.csv"),quote=F,append=TRUE)
write.table(t(cellss2),row.names=F,col.names=F,file=paste0(commands[8],"cell_corr.csv"),quote=F,append=TRUE)
res=cbind(colnames(expma),cellss1,cellss2)
colnames(res)=c("Patient ID","Pivot cell","Partner cell")
write.table(res,file=paste0(commands[8],"plot.csv"),sep='\t',quote=F,row.names=F)

##################### survival interaction
colnames(clinical)=gsub("\\."," ",colnames(clinical))
sur2=c("Survival Time","Death")
surv = Surv(as.numeric(clinical[,match(sur2[1],colnames(clinical))]),as.numeric(clinical[,match(sur2[2],colnames(clinical))]))

cellsh1=rep(1,length(cellss1))
cellsh1[which(cellss1>median(cellss1))]=2
cellsh2=rep(1,length(cellss2))
cellsh2[which(cellss2>median(cellss2))]=2

Interaction=cellsh1 * cellsh2
B3=cofac_vars(commands[7],clinical,cbind(cellsh1,cellsh2,Interaction))
coxout=coxcal(surv,B3)
result_interaction = coxout["Interaction", c("z", "Pr(>|z|)")]
#write.table(t(as.matrix(result_pg)),file=paste0(commands[8],"xls.txt"),sep='\t',quote=F,row.names=F)
names(result_interaction)=c("z.Score","p.value")
write.csv(t(result_interaction),file=paste0(commands[8],"result_interaction.csv"),row.names=F,quote=F)

B1 = cofac_vars(commands[7],clinical,cellsh1)
coxout=coxcal(surv,B1)
resultcell_pz = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]

if(!file.exists(paste0(commands[8],"lable2")))	dir.create(path=paste0(commands[8],"lable2"))
if(!file.exists(paste0(commands[8],"lable2/pic1")))	dir.create(path=paste0(commands[8],"lable2/pic1"))
survout(surv,cellsh1,paste0(commands[8],"lable2/pic1/"),"Pivot Cell Low|Pivot Cell High")
write.table(c("Total Cohort","Pivot Cell Low|Pivot Cell High",resultcell_pz),file=paste0(commands[8],"lable2/pic1/pic.txt"),row.names=F,quote=F)
res=cbind(cellss1,clinical)
colnames(res)=c("Variable",colnames(clinical))
write.table(res,file=paste0(commands[8],"lable2/pic1/download.txt"),sep='\t',quote=F)

############################################
idgh=which(cellsh2==2)
coxout=coxcal(surv[idgh,],as.data.frame(B1[idgh,]))
resultcell_pz = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]

if(!file.exists(paste0(commands[8],"lable2")))	dir.create(path=paste0(commands[8],"lable2"))
if(!file.exists(paste0(commands[8],"lable2/pic2")))	dir.create(path=paste0(commands[8],"lable2/pic2"))
survout(surv[idgh,],cellsh1[idgh],paste0(commands[8],"lable2/pic2/"),"Pivot Cell Low|Pivot Cell High")
write.table(c("Partner Cell High Subcohort","Pivot Cell Low|Pivot Cell High",resultcell_pz),file=paste0(commands[8],"lable2/pic2/pic.txt"),row.names=F,quote=F)
res=cbind(cellss1[idgh],clinical[idgh,])
colnames(res)=c("Variable",colnames(clinical))
write.table(res,file=paste0(commands[8],"lable2/pic2/download.txt"),,sep='\t',quote=F)

###########################################
idgl=which(cellsh2==1)
coxout=coxcal(surv[idgl,],as.data.frame(B1[idgl,]))
resultcell_pz = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]

if(!file.exists(paste0(commands[8],"lable2")))	dir.create(path=paste0(commands[8],"lable2"))
if(!file.exists(paste0(commands[8],"lable2/pic3")))	dir.create(path=paste0(commands[8],"lable2/pic3"))
survout(surv[idgl,],cellsh1[idgl],paste0(commands[8],"lable2/pic3/"),"Pivot Cell Low|Pivot Cell High")
write.table(c("Partner Cell Low Subcohort","Pivot Cell Low|Pivot Cell High",resultcell_pz),file=paste0(commands[8],"lable2/pic3/pic.txt"),row.names=F,quote=F)
res=cbind(cellss1[idgl],clinical[idgl,])
colnames(res)=c("Variable",colnames(clinical))
write.table(res,file=paste0(commands[8],"lable2/pic3/download.txt"),sep='\t',quote=F)

############################################
B2 = cofac_vars(commands[7],clinical,cellsh2)
coxout=coxcal(surv,B2)
resultcell_pz = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]

if(!file.exists(paste0(commands[8],"lable3")))	dir.create(path=paste0(commands[8],"lable3"))
if(!file.exists(paste0(commands[8],"lable3/pic1")))	dir.create(path=paste0(commands[8],"lable3/pic1"))
survout(surv,cellsh2,paste0(commands[8],"lable3/pic1/"),"Partner Cell Low|Partner Cell High")
write.table(c("Total Cohort","Partner Cell Low|Partner Cell High",resultcell_pz),file=paste0(commands[8],"lable3/pic1/pic.txt"),row.names=F,quote=F)
res=cbind(cellss2,clinical)
colnames(res)=c("Variable",colnames(clinical))
write.table(res,file=paste0(commands[8],"lable3/pic1/download.txt"),sep='\t',quote=F)

############################################
idgh=which(cellsh1==2)
coxout=coxcal(surv[idgh,],as.data.frame(B2[idgh,]))
resultcell_pz = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]

if(!file.exists(paste0(commands[8],"lable3/pic2")))	dir.create(path=paste0(commands[8],"lable3/pic2"))
survout(surv[idgh,],cellsh2[idgh],paste0(commands[8],"lable3/pic2/"),"Partner Cell Low|Partner Cell High")
write.table(c("Pivot Cell High Subcohort","Partner Cell Low|Partner Cell High",resultcell_pz),file=paste0(commands[8],"lable3/pic2/pic.txt"),row.names=F,quote=F)
res=cbind(cellss2[idgh],clinical[idgh,])
colnames(res)=c("Variable",colnames(clinical))
write.table(res,file=paste0(commands[8],"lable3/pic2/download.txt"),,sep='\t',quote=F)

###########################################
idgl=which(cellsh1==1)
coxout=coxcal(surv[idgl,],as.data.frame(B2[idgl,]))
resultcell_pz = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]

if(!file.exists(paste0(commands[8],"lable3/pic3")))	dir.create(path=paste0(commands[8],"lable3/pic3"))
survout(surv[idgl,],cellsh2[idgl],paste0(commands[8],"lable3/pic3/"),"Partner Cell Low|Partner Cell High")
write.table(c("Pivot Cell Low Subcohort","Partner Cell Low|Partner Cell High",resultcell_pz),file=paste0(commands[8],"lable3/pic3/pic.txt"),row.names=F,quote=F)
res=cbind(cellss2[idgl],clinical[idgl,])
colnames(res)=c("Variable",colnames(clinical))
write.table(res,file=paste0(commands[8],"lable3/pic3/download.txt"),sep='\t',quote=F)

##############################################
class4=rep(1,MM)
class4[intersect(which(cellsh2==2),which(cellsh1==1))]=2
class4[intersect(which(cellsh2==1),which(cellsh1==2))]=3
class4[intersect(which(cellsh2==2),which(cellsh1==2))]=4
B2 = cofac_vars(commands[7],clinical,class4)
coxout=coxcal(surv,B2)
resultcell_pz = coxout[dim(coxout)[1], c("exp(coef)","z", "Pr(>|z|)")]

labels4="Pivot Cell Low & Partner Cell Low|Pivot Cell Low & Partner Cell High|Pivot Cell High & Partner Cell Low|Pivot Cell High & Partner Cell High"
if(!file.exists(paste0(commands[8],"lable1")))	dir.create(path=paste0(commands[8],"lable1"))
survout(surv,class4,paste0(commands[8],"lable1/"),labels4)
write.table(c("Overall Survival of 4 Groups",labels4,resultcell_pz),file=paste0(commands[8],"lable1/pic.txt"),row.names=F,quote=F)
res=cbind(class4,clinical)
colnames(res)=c("Variable",colnames(clinical))
write.table(res,file=paste0(commands[8],"lable1/download.txt"),sep='\t',quote=F)

write.table(c(),file=paste0(commands[8],"end.txt"),row.names=F,col.names=F)