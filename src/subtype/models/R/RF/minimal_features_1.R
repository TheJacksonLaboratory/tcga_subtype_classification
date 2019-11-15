library(xtable)
library(Biobase)
library(class)
library(randomForest)
library(xtable)
library(party)
library(caret)
library(e1071)
library(gbm)

set.seed= 1024

source("../auxilary_functions.R")
load("../../../../data/raw/combined_data.rdata")

eset.1 <- combined.eset[,combined.eset$sample.type %in% c("01","03")]
eset.2 <- combined.eset[,combined.eset$tumor.type =="SKCM" & combined.eset$sample.type == "06"]
eset <- Biobase::combine(eset.1,eset.2)

tumor.types <- unique(eset$tumor.type)

exp <- Biobase::exprs(eset)
exp.log2 <- log2(exp +1)
r.max <- apply(exp.log2,1,max)
r.sd  <- apply(exp.log2,1,sd)
exp.sel.all <- exp.log2[r.max > 8 & r.sd > 1,]
p.anns <- featureData(eset)@data
probe.anns <- p.anns[rownames(exp.sel.all),]
probe.anns.sel <- probe.anns[probe.anns$symbol != "?",]
exp.sel <- exp.sel.all[rownames(probe.anns.sel),]
rownames(exp.sel) <- probe.anns.sel$symbol
true.label <- factor(eset$tumor.type)
all.dat <- data.frame(t(exp.sel) ,tumor.type = eset$tumor.type)
exp.t = data.frame(t(exp.sel))

load("final_model.rdata")


dat <- data.frame(symbol = rownames(final.model$importance),final.model$importance)
dat.ord <- dat[order(dat$MeanDecreaseGini,decreasing=T),]

num = 50 
err.rate =numeric(length = num)
for(i in 1:num)
{
	len = 5*i 
	features <- rownames(dat.ord)[1:len]
	train.dat <- data.frame(exp.t[,features], tumor.type = eset$tumor.type)	
	temp.model <- randomForest(tumor.type  ~ . , data = train.dat ,ntree= 100)
	err.rate[i] <- apply(temp.model$err.rate,2,sum)[1]
	print(i)
}




save(temp.model, file = "model_with_250_features.rdata")

