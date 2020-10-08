setwd('/Users/divyadesadla/Desktop/CMU/CMUFall19/Data Science for Biological Sciences/Homework1')
tbrData <- read.delim("TbrMappedCountMatrix.txt",row.names="Gene")

tbrDataNorm <- tbrData
for (i in 1:dim(tbrData)[2])
{
  tbrDataNorm[,i] <- tbrData[,i]/sum(tbrData[,i])
}

tbrDataTrimmed <- tbrDataNorm[((tbrData$Control1 > 5)+(tbrData$Control2 > 5)+(tbrData$Control3 > 5)+
                                   +                                    (tbrData$Tbr1 > 5)+(tbrData$Tbr2 > 5)+(tbrData$Tbr3 > 5))>2,]
tbrDataGenes <- tbrDataTrimmed[substr(row.names(tbrDataTrimmed),1,3)=="PMI",]

tbrDataMatrix <- sapply(tbrDataGenes,cbind)

tbrDataMatrix_transposed = t(tbrDataMatrix)

kmeans2 <- kmeans(tbrDataMatrix_transposed,centers=2,iter.max=10,nstart=25)
hist(kmeans2$cluster)

kmeans2$cluster

kmeans2_1 <- kmeans(tbrDataMatrix_transposed,centers=2,iter.max=10,nstart=25)
kmeans2_2 <- kmeans(tbrDataMatrix_transposed,centers=2,iter.max=10,nstart=25)
kmeans2_3 <- kmeans(tbrDataMatrix_transposed,centers=2,iter.max=10,nstart=25)
kmeans2_4 <- kmeans(tbrDataMatrix_transposed,centers=2,iter.max=10,nstart=25)
kmeans2_5 <- kmeans(tbrDataMatrix_transposed,centers=2,iter.max=10,nstart=25)

kmeans2_1$cluster
kmeans2_2$cluster
kmeans2_3$cluster
kmeans2_4$cluster
kmeans2_5$cluster

kmeans3_1 <- kmeans(tbrDataMatrix_transposed,centers=3,iter.max=10,nstart=25)
kmeans3_2 <- kmeans(tbrDataMatrix_transposed,centers=3,iter.max=10,nstart=25)
kmeans3_3 <- kmeans(tbrDataMatrix_transposed,centers=3,iter.max=10,nstart=25)
kmeans3_4 <- kmeans(tbrDataMatrix_transposed,centers=3,iter.max=10,nstart=25)
kmeans3_5 <- kmeans(tbrDataMatrix_transposed,centers=3,iter.max=10,nstart=25)

kmeans3_1$cluster
kmeans3_2$cluster
kmeans3_3$cluster
kmeans3_4$cluster
kmeans3_5$cluster

library("edgeR")
tbrData <- read.delim("TbrMappedCountMatrix.txt",row.names="Gene")
tbrDataTrimmed <- subset(tbrData,(Control1>5)+(Control2>5)+(Control3>5)+(Tbr1>5)+(Tbr2>5)+(Tbr3>5)>2)
dim(tbrDataTrimmed)
group <- factor(c(1,2,3,1,2,3))
design <- model.matrix(~group)
dge <- DGEList(counts=tbrDataTrimmed,group=group)
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge,design)
fit <- glmFit(dge,design)
lrt <- glmLRT(fit,coef=2)
edgeResults <- topTags(lrt,n=10000,p.value=0.005/dim(tbrData)[1])
dim(edgeResults)
edgeResults

