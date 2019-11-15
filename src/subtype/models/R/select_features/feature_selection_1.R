library(Biobase)
library(class)
library(randomForest)
library(xtable)
library(party)
library(caret)
library(e1071)
library(gbm)
library(sparsediscrim)

set.seed= 1024

source("../auxilary_functions.R")
load("../../../../data/raw/combined_data.rdata")

eset.1 <- combined.eset[,combined.eset$sample.type %in% c("01","03")]
eset.2 <- combined.eset[,combined.eset$tumor.type =="SKCM" & combined.eset$sample.type == "06"]
eset <- Biobase::combine(eset.1,eset.2)

tumor.types <- unique(eset$tumor.type)

exp <- exprs(eset)
exp.log2 <- log2(exp +1)
r.max <- apply(exp.log2,1,max)
r.sd  <- apply(exp.log2,1,sd)
exp.sel <- exp.log2[r.max > 8 & r.sd > 1,]

true.label <- eset$tumor.type
all.dat <- data.frame(t(exp.sel) ,tumor.type = eset$tumor.type)

features <- select_features_for_training(all.dat)


