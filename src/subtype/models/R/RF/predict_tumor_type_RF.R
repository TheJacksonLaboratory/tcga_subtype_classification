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

source("src/models/auxilary_functions.R")
load("data/raw/combined_data.rdata")

eset.1 <- combined.eset[,combined.eset$sample.type %in% c("01","03")]
eset.2 <- combined.eset[,combined.eset$tumor.type =="SKCM" & combined.eset$sample.type == "06"]
eset <- Biobase::combine(eset.1,eset.2)

tumor.types <- unique(eset$tumor.type)

exp <- Biobase::exprs(eset)
exp.log2 <- log2(exp +1)
r.max <- apply(exp.log2,1,max)
r.sd  <- apply(exp.log2,1,sd)
exp.sel <- exp.log2[r.max > 8 & r.sd > 1,]

true.label <- factor(eset$tumor.type)
all.dat <- data.frame(t(exp.sel) ,tumor.type = eset$tumor.type)

K = 3
flds <- createFolds(1:nrow(all.dat), k = K, list = TRUE, returnTrain = FALSE)
RF.predicted.label = character(length = nrow(all.dat))
cl.idx <- which(colnames(all.dat) == "tumor.type")

for(i in 1:K)
{
	test.idx <- flds[[i]]
	train.all <- all.dat[-test.idx,]
	features <- select_features_for_training(train.all)
	test.all <- all.dat[test.idx,-c(cl.idx)]
	train <- train.all[,c(features,"tumor.type")]
	test  <- test.all[,features]
	model <- randomForest(tumor.type  ~ ., data = train, ntree = 1000, mtry = 31)
	RF.pred <- predict(model, newdata = test)
	RF.predicted.label[test.idx] <- as.character(RF.pred)
}

RF.predicted.label = factor(RF.predicted.label)
RF.res <- confusionMatrix(RF.predicted.label, true.label)
save(RF.res, file = "models/primary_type/final_model_cv.rdata")


# final model
set.seed= 1024

exp.sel.all <- exp.log2[r.max > 8 & r.sd > 1,]
p.anns <- featureData(eset)@data
probe.anns <- p.anns[rownames(exp.sel.all),]
probe.anns.sel <- probe.anns[probe.anns$symbol != "?",]
exp.sel <- exp.sel.all[rownames(probe.anns.sel),]
rownames(exp.sel) <- probe.anns.sel$symbol
true.label <- factor(eset$tumor.type)
all.dat <- data.frame(t(exp.sel) ,tumor.type = eset$tumor.type)

final.model <- randomForest(tumor.type  ~ ., data = all.dat, ntree = 1000, mtry = 31)
save(final.model, file = "models/primary_type/final_model.rdata")
