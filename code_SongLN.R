########################################################################################
# Consensus clustering

library(limma)
library(ConsensusClusterPlus)
expFile="ssGSEAresult.txt"      
workDir="G:\\BioData\\mir_SLN_final_20220509\\2.a.ssGSEAcluster"     
setwd(workDir)      

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

maxK=9
results=ConsensusClusterPlus(data,
              maxK=maxK,
              reps=50,
              pItem=0.8,
              pFeature=1,
              title=workDir,
              clusterAlg="km",
              distance="euclidean",
              seed=123456,
              plot="pdf")

clusterNum=3     
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("ssGSEAcluster")
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(cluster$ssGSEAcluster))
cluster$ssGSEAcluster=letter[match(cluster$ssGSEAcluster, uniqClu)]
clusterOut=rbind(ID=colnames(cluster), cluster)
write.table(clusterOut, file="ssGSEAcluster.txt", sep="\t", quote=F, col.names=F)



########################################################################################
# Corelation analysis between m6A regulators and immune-related genes

library(limma)
corFilter=0.5          
pvalueFilter=0.001     
setwd("E:\\Science\\TCGA_LUAD_m6ATME_Final\\20.coExpression")  

rt=read.table("ssGSEAdiffGeneExp_n59t500.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.1,]

group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
ssGSEA=data[,group==0]
conNum=length(group[group==1])      
treatNum=length(group[group==0])   
sampleType=c(rep(1,conNum), rep(2,treatNum))

rt1=read.table("m6aGeneExp_n59t500.txt", header=T, sep="\t", check.names=F)
rt1=as.matrix(rt1)
rownames(rt1)=rt1[,1]
exp1=rt1[,2:ncol(rt1)]
dimnames1=list(rownames(exp1),colnames(exp1))
m6A=matrix(as.numeric(as.matrix(exp1)), nrow=nrow(exp1), dimnames=dimnames1)
m6A=avereps(m6A)
m6A=m6A[rowMeans(m6A)>0.1,]

group=sapply(strsplit(colnames(m6A),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
m6A=m6A[,group==0]

outTab=data.frame()
for(i in row.names(ssGSEA)){
	if(sd(ssGSEA[i,])>0.1){
		test=wilcox.test(data[i,] ~ sampleType)
		if(test$p.value<0.05){
			for(j in row.names(m6A)){
				x=as.numeric(ssGSEA[i,])
				y=as.numeric(m6A[j,])
				corT=cor.test(x,y)
				cor=corT$estimate
				pvalue=corT$p.value
				if((cor>corFilter) & (pvalue<pvalueFilter)){
					outTab=rbind(outTab,cbind(m6A=j,ssGSEA=i,cor,pvalue,Regulation="postive"))
				}
				if((cor< -corFilter) & (pvalue<pvalueFilter)){
					outTab=rbind(outTab,cbind(m6A=j,ssGSEA=i,cor,pvalue,Regulation="negative"))
				}
			}
		}
	}
}

write.table(file="net.network.txt",outTab,sep="\t",quote=F,row.names=F)
lncNode=data.frame(Node=unique(as.vector(outTab[,"ssGSEA"])), Type="ssGSEA")
mrnaNode=data.frame(Node=unique(as.vector(outTab[,"m6A"])), Type="m6A")
nodeOut=rbind(lncNode, mrnaNode)
write.table(nodeOut, file="net.node.txt", sep="\t", quote=F, row.names=F)

m6assGSEA=unique(as.vector(outTab[,"ssGSEA"]))
m6assGSEAexp=data[m6assGSEA,]
m6assGSEAexp=rbind(ID=colnames(m6assGSEAexp), m6assGSEAexp)
write.table(m6assGSEAexp,file="mirGeneExp_n59t500.txt",sep="\t",quote=F,col.names=F)

########################################################################################
# PCA score

expFile="mirGeneExp.txt"     
setwd("E:\\Science\\TCGA_LUAD_m6ATME_Final\\22.PCAscore")  

data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)

pca=prcomp(data, scale=TRUE)
value=predict(pca)
PCAscore=value[,1]+value[,2]
PCAscore=as.data.frame(PCAscore)
scoreOut=rbind(id=colnames(PCAscore), PCAscore)
write.table(scoreOut, file="PCAscore.txt", sep="\t", quote=F, col.names=F)

#########################################################################################
#K-M curve (grouped by surv_cutpoint)

library(survival)
library(survminer)
scoreFile="PCAscore.txt"    
cliFile="time.txt"      
setwd("E:\\Science\\TCGA_LUAD_m6ATME_Final\\23.PCAscoresur")   

score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
sampleType=gsub("(.*?)\\_.*", "\\1", row.names(score))
score=cbind(score, sampleType)
rownames(score)=gsub("(.*?)\\_(.*?)", "\\2", rownames(score))
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/365

sameSample=intersect(row.names(score), row.names(cli))
data=cbind(cli[sameSample,], score[sameSample,])

res.cut=surv_cutpoint(data, time="futime", event="fustat", variables=c("PCAscore"))
cutoff=as.numeric(res.cut$cutpoint[1])
print(cutoff)
Type=ifelse(data[,"PCAscore"]<=cutoff, "Low", "High")
data$group=Type
outTab=rbind(id=colnames(data), data)
write.table(outTab, file="PCAscore.group.txt", sep="\t", quote=F, col.names=F)

data$group=factor(data$group, levels=c("Low", "High"))
diff=survdiff(Surv(futime, fustat) ~ group, data = data)
length=length(levels(factor(data[,"group"])))
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = data)
#print(surv_median(fit))

bioCol=c("#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
			       data=data,
			       conf.int=F,
			       pval=pValue,
			       pval.size=6,
			       legend.title="PCAscore",
			       legend.labs=levels(factor(data[,"group"])),
			       legend = c(0.8, 0.8),
			       font.legend=12,
			       xlab="Time(years)",
			       break.time.by = 1,
			       palette = bioCol,
			       surv.median.line = "hv",
			       risk.table=T,
			       cumevents=F,
			       risk.table.height=.25)

pdf(file="survival.pdf", onefile = FALSE, width=5, height=5)
print(surPlot)
dev.off()


#########################################################################################
# Differential Analysis (Boxplot)

library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)

setwd("G:\\BioData\\m6ATME_SLN\\m6ATME_SLN\\4.a-d.boxplot")   
ssgseaScore=read.table("cibersort_t500.txt", header=T, sep="\t", check.names=F, row.names=1)
ssgseaScore[1:3,1:3]
dim(ssgseaScore)


cluster=read.table("mirCluster.txt", header=T, sep="\t", check.names=F, row.names=1)
cluster[1:3,]
dim(cluster)
class(cluster)

ssgseaScore=t(ssgseaScore)
ssgseaScore[1:3,1:3]
dim(ssgseaScore)

sameSample=intersect(row.names(ssgseaScore), row.names(cluster))
length(sameSample)
ssgseaScore=ssgseaScore[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
scoreCluster=cbind(ssgseaScore, cluster)
dim(scoreCluster)
scoreCluster[1:3,]

data=melt(scoreCluster, id.vars=c("mirCluster"))
colnames(data)=c("mirCluster", "Immune", "Fraction")
class(data)
data[1:3,]
data

bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data[,"mirCluster"])))]
p=ggboxplot(data, x="Immune", y="Fraction",color= "mirCluster",
     ylab="immune infiltration",
     xlab="",
     legend.title="mirCluster",
     palette=bioCol)
p=p+rotate_x_text(50)
pdf(file="boxplot_cibersort.pdf", width=10, height=7.5)                   
p+stat_compare_means(aes(group=mirCluster),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()

#########################################################################################
# Analysis of drug sensitivity

library(limma)
library(ggpubr)
library(pRRophetic)
library(ggplot2)
set.seed(12345)
expFile="mRNAsymbol_n59t500.txt"     
riskFile="cluster.txt"      
drug="ZM.447439"         
setwd("G:\\BioData\\m6ATME_SLN\\3-420220209\\04.pRRophetic")   

rt = read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
data=data[,group==0]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)
data=t(data)

senstivity=pRRopheticPredict(data, drug, selection=1)
senstivity=senstivity[senstivity!="NaN"]
#senstivity[senstivity>quantile(senstivity,0.99)]=quantile(senstivity,0.99)

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

sameSample=intersect(row.names(risk), names(senstivity))
risk=risk[sameSample, "risk",drop=F]
senstivity=senstivity[sameSample]
rt=cbind(risk, senstivity)

rt$risk=factor(rt$risk, levels=c("A", "B", "C"))
type=levels(factor(rt[,"risk"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

boxplot=ggboxplot(rt, x="risk", y="senstivity", fill="risk",
                  xlab="mirCluster",
                  ylab=paste0(drug, " senstivity (IC50)"),
                  legend.title="mirCluster",
                  palette=c("blue", "orange", "red")
)+ 
  stat_compare_means(comparisons=my_comparisons)
pdf(file=paste0(drug, ".pdf"), width=5, height=4.5)
print(boxplot)
dev.off()

