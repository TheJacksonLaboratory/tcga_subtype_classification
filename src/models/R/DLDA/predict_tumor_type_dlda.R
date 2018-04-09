library(xtable)
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

K = 3 
flds <- createFolds(1:nrow(all.dat), k = K, list = TRUE, returnTrain = FALSE)
DLDA.predicted.label = character(length = nrow(all.dat))
cl.idx <- which(colnames(all.dat) == "tumor.type")

for(i in 1:K)
{
	test.idx <- flds[[i]]
	train.all <- all.dat[-test.idx,]
	features <- select_features_for_training(train.all)
	test.all <- all.dat[test.idx,-c(cl.idx)]
	train <- train.all[,c(features,"tumor.type")]
	test  <- test.all[,features]
	t.cl.idx <- which(colnames(train) == "tumor.type")
        train.exp <- train[,-c(t.cl.idx)]
        train.cl <- train[,t.cl.idx]
	train.exp <- train[,-c(t.cl.idx)]
	train.cl <- train[,t.cl.idx]
	dlda_out = dlda(x = train.exp,y = train.cl)
	dlda.pred <- predict(dlda_out, test)$class
	DLDA.predicted.label[test.idx] <- as.character(dlda.pred)
}

DLDA.predicted.label = factor(DLDA.predicted.label)
DLDA.res <- confusionMatrix(DLDA.predicted.label,true.label)

write.csv(DLDA.res[[2]],"contingency_table_DLDA.csv")
write.csv(DLDA.res[[4]][,1:7],"performance_metrics_DLDA.csv")


