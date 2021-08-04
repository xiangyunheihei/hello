# hello
PERP
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma", version = "3.8")

library("limma")

setwd("C:\\Users\\windows\\Desktop\\1")             
gene="PERP"                                                           
normalNum=59                                                         
tumorNum=535                                                         

rt=read.table("symbol.txt",sep="\t",header=T,check.names=F)         
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
uniq=rbind(ID=colnames(data),data)
write.table(uniq,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)       
Type=c(rep("Normal",normalNum),rep("Tumor",tumorNum))
single=cbind(ID=colnames(data),expression=data[gene,],Type)
colnames(single)=c("ID",gene,"Type")
write.table(single,file="singleGene.txt",sep="\t",quote=F,row.names=F)

#install.packages("beeswarm")
setwd("C:\\Users\\windows\\Desktop\\1")        
inputFile="singleGene.txt"                                     
yMin=0                     
yMax=800                     
ySeg=yMax*0.94

library(limma)
library(beeswarm)

rt=read.table(inputFile,sep="\t",header=T,row.names=1,check.names=F)
geneName=colnames(rt)[1]
labels=c("Normal","Tumor")
colnames(rt)=c("expression","Type")
wilcoxTest<-wilcox.test(expression ~ Type, data=rt)
wilcoxP=wilcoxTest$p.value
pvalue=signif(wilcoxP,4)
pval=0
if(pvalue<0.001){
     pval=signif(pvalue,4)
     pval=format(pval, scientific = TRUE)
}else{
     pval=round(pvalue,3)
}

outFile=paste(geneName,".pdf",sep="")
pdf(file=outFile,width=7,height=5)
par(mar = c(4,7,3,3))
boxplot(expression ~ Type, data = rt,names=labels,
     ylab = paste(geneName," expression",sep=""),
     cex.main=1.5, cex.lab=1.3, cex.axis=1.2,ylim=c(yMin,yMax),outline = FALSE)
beeswarm(expression ~ Type, data = rt, col = c("blue","red"),lwd=0.1,
     pch = 16, add = TRUE, corral="wrap")
segments(1,ySeg, 2,ySeg);segments(1,ySeg, 1,ySeg*0.96);segments(2,ySeg, 2,ySeg*0.96)
text(1.5,ySeg*1.05,labels=paste("p=",pval,sep=""),cex=1.2)
dev.off()

#install.packages("survival")

setwd("C:\\Users\\windows\\Desktop\\1")
gene="PERP"

library(survival)
rt=read.table("survival.txt",header=T,sep="\t",check.names=F)
rt$futime=rt$futime/365                                    

a=rt[,gene]<=median(rt[,gene])
diff=survdiff(Surv(futime, fustat) ~a,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
     pValue=signif(pValue,4)
     pValue=format(pValue, scientific = TRUE)
}else{
     pValue=round(pValue,3)
}

fit <- survfit(Surv(futime, fustat) ~ a, data = rt)

pdf(file="survival.pdf",
    width=6,
    height=6)
plot(fit, 
     lwd=2,
     col=c("red","blue"),
     xlab="Time (year)",
     mark.time=T,
     ylab="Survival rate",
     main=paste(gene,"(p=", pValue ,")",sep="") )
legend("topright", 
     c("High expression","Low expression"), 
     lwd=2, 
     col=c("red","blue"))
dev.off()
summary(fit)

inputFile="singleGeneClinical.txt"                                 
setwd("C:\\Users\\windows\\Desktop\\1")      
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
clinical="stage"                                                      
gene="PERP"                                                          

xlabel=vector()
tab1=table(rt[,clinical])
labelNum=length(tab1)
dotCol=c(2,3)
if(labelNum==3){
	dotCol=c(2,3,4)
}
if(labelNum==4){
	dotCol=c(2,3,4,5)
}
if(labelNum>4){
	dotCol=rainbow(labelNum)
}
for(i in 1:labelNum){
  xlabel=c(xlabel,names(tab1[i]) )
}

outTab=data.frame()

rt1=rbind(expression=rt[,gene],clinical=rt[,clinical])
rt1=as.matrix(t(rt1))
if(labelNum==2){
    wilcoxTest<-wilcox.test(expression ~ clinical, data=rt1)
}else{
    wilcoxTest<-kruskal.test(expression ~ clinical, data = rt1)}
pValue=wilcoxTest$p.value
outTab=rbind(outTab,cbind(gene=gene,pVal=pValue))
pval=0
if(pValue<0.001){
    pval=signif(pValue,4)
    pval=format(pval, scientific = TRUE)
}else{
    pval=round(pValue,3)
}
  
b = boxplot(expression ~ clinical, data = rt1,outline = FALSE, plot=F) 
yMin=min(b$stats)
yMax = max(b$stats/5+b$stats)
ySeg = max(b$stats/10+b$stats)
ySeg2 = max(b$stats/12+b$stats)
n = ncol(b$stats)

outPdf=paste0(gene,".",clinical,".pdf")
pdf(file=outPdf,
    width=9,
    height=6,)
par(mar = c(4,7,3,3))
boxplot(expression ~ clinical, data = rt1,names=xlabel,
    ylab = paste0(gene," expression"),col=dotCol,
    cex.main=1.6, cex.lab=1.4, cex.axis=1.3,ylim=c(yMin,yMax),outline = FALSE)
segments(1,ySeg, n,ySeg);
segments(1,ySeg, 1,ySeg2)
segments(n,ySeg, n,ySeg2)
text((1+n)/2,ySeg,labels=paste0("p=",pval),cex=1.5,pos=3)
dev.off()

inputFile="singleGeneClinical.txt"                                     
setwd("C:\\Users\\windows\\Desktop\\1")              
rt=read.table(inputFile,sep="\t",header=T,check.names=F)                 
clinical="stage"
gene="PERP"
data=rt[,c(clinical,gene)]
colnames(data)=c("clinical","gene")
data=cbind(data,y=ifelse(rt[,gene]>median(rt[,gene]),1,0))
logit=glm(y~clinical,family=binomial(link='logit'),data=data)
summ=summary(logit)
conf=confint(logit,level = 0.95)
cbind(OR=exp(summ$coefficients[2,"Estimate"]),
      OR.95L=exp(conf[2,1]),
      OR.95H=exp(conf[2,2]),
      p=summ$coefficients[2,"Pr(>|z|)"])
#install.packages('survival')

setwd("C:\\Users\\windows\\Desktop\\1")
library(survival)
rt=read.table("coxInput.txt",header=T,sep="\t",check.names=F,row.names=1)


outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
 cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
 coxSummary = summary(cox)
 coxP=coxSummary$coefficients[,"Pr(>|z|)"]
 outTab=rbind(outTab,
              cbind(id=i,
              HR=coxSummary$conf.int[,"exp(coef)"],
              HR.95L=coxSummary$conf.int[,"lower .95"],
              HR.95H=coxSummary$conf.int[,"upper .95"],
              pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
              )
}
write.table(outTab,file="uniCox.xls",sep="\t",row.names=F,quote=F)

#install.packages('survival')
#install.packages("survminer")

library(survival)
library(survminer)
setwd("C:\\Users\\windows\\Desktop\\1")
rt=read.table("coxInput.txt",header=T,sep="\t",check.names=F,row.names=1)


multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
             HR=multiCoxSum$conf.int[,"exp(coef)"],
             HR.95L=multiCoxSum$conf.int[,"lower .95"],
             HR.95H=multiCoxSum$conf.int[,"upper .95"],
             pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)

pdf(file="forest.pdf",
       width = 7,             
       height = 6,         
       )
ggforest(multiCox,
         main = "Hazard ratio",
         cpositions = c(0.02,0.22, 0.4), 
         fontsize = 0.7, 
         refLabel = "reference", 
         noDigits = 2)
dev.off()

#install.packages("ggplot2")

library(plyr)
library(ggplot2)
library(grid)
library(gridExtra)

setwd("C:\\Users\\windows\\Desktop\\1")          
files=grep(".xls",dir(),value=T)                                    
data = lapply(files,read.delim)                                        
names(data) = files

dataSet = ldply(data, data.frame)
dataSet$pathway = gsub(".xls","",dataSet$.id)                            

gseaCol=c("#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
pGsea=ggplot(dataSet,aes(x=RANK.IN.GENE.LIST,y=RUNNING.ES,colour=pathway,group=pathway))+
  geom_line(size = 1.5) + scale_color_manual(values = gseaCol[1:nrow(dataSet)]) +   
  labs(x = "", y = "Enrichment Score", title = "") + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0),limits = c(min(dataSet$RUNNING.ES - 0.02), max(dataSet$RUNNING.ES + 0.02))) +   
  theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_blank()) + theme(axis.line = element_line(colour = "black")) + theme(axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + 
  geom_hline(yintercept = 0) +   theme(legend.position = c(0,0),legend.justification = c(0,0)) + #legend注释的位值
  guides(colour = guide_legend(title = NULL)) + theme(legend.background = element_blank()) + theme(legend.key = element_blank())+theme(legend.key.size=unit(0.5,'cm'))
pGene=ggplot(dataSet,aes(RANK.IN.GENE.LIST,pathway,colour=pathway))+geom_tile()+
  scale_color_manual(values = gseaCol[1:nrow(dataSet)]) + 
  labs(x = "high expression<----------->low expression", y = "", title = "") + 
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +  
  theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_blank()) + theme(axis.line = element_line(colour = "black"))+
  theme(axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())+ guides(color=FALSE)

gGsea = ggplot_gtable(ggplot_build(pGsea))
gGene = ggplot_gtable(ggplot_build(pGene))
maxWidth = grid::unit.pmax(gGsea$widths, gGene$widths)
gGsea$widths = as.list(maxWidth)
gGene$widths = as.list(maxWidth)
dev.off()

pdf('multipleGSEA.pdf',     
     width=7,                
     height=5.5)             
par(mar=c(5,5,2,5))
grid.arrange(arrangeGrob(gGsea,gGene,nrow=2,heights=c(.8,.3)))
dev.off()
