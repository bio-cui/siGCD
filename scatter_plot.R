#!/usr/bin/env Rscript

commands=commandArgs(trailingOnly=T)
file_name = commands[4]
# file_name_sub = sub(pattern="/",replacement="", file_name)
# if(!file.exists(file_name_sub))	dir.create(path=file_name_sub)

source("funcs.R")
cellss=read.table(paste0(commands[4],"cell_exp.csv"))[[1]]
#print(commands[3])
datalist=loaddata(commands[2],commands[3])

expma=datalist[[1]]
MM = dim(expma)[2]
gene_or_cell=commands[6]
if(gene_or_cell == "gene" | gene_or_cell == "cell" | gene_or_cell =="DrugGene"){
	geneintid = match(commands[1],rownames(expma))
	res=cbind(colnames(expma),cellss,expma[geneintid,])
}else if(commands[5]=="multi"){
	drugbank=read.table("drugbank.txt",sep='\t')
	targets=paste(as.character(drugbank[which(as.character(drugbank[,1])==commands[1]),2]),collapse=',')
	cellss2=cellcal(targets,"none",expma,MM)
	res=cbind(colnames(expma),cellss,cellss2)
}else{
	cellss2=cellcal(commands[1],"none",expma,MM)
	res=cbind(colnames(expma),cellss,cellss2)
}

colnames(res)=c("Patient ID","Pivot cell",commands[1])
write.table(res,file=paste0(commands[4],"plot.csv"),sep='\t',row.names=F,quote=F)
