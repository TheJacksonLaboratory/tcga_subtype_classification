library(xtable)
library(Biobase)
library(class)
library(randomForest)
library(party)
library(caret)
library(e1071)
library(gbm)

set.seed= 1024

load("data/raw/combined_data_with_subtype_anns.rdata")

source("src/models/R/auxilary_functions.R")
sample.anns <- pData(eset)
tab <- table(sample.anns$tumor.type,sample.anns$subtype)
tab.bin <- tab
tab.bin[tab>0] <- 1
tab.sel <- tab[apply(tab.bin,1,sum) > 1,]
tumors.all <- rownames(tab.sel)
tumors = setdiff(tumors.all,c("GBM","UCEC"))
probe.anns = featureData(eset)@data
RF.res <- list()
model.list <- list()

for(i in 1:length(tumors))
{
	subtypes = table(sample.anns[sample.anns$tumor.type == tumors[i],"subtype"])
	subtype.names <- setdiff(names(subtypes[subtypes > 25]),"")
	if(length(subtype.names) > 1)
	{
		eset.sel <- eset[,eset$tumor.type == tumors[i] & eset$subtype %in% subtype.names]
		K= 3
		exp.sel <- Biobase::exprs(eset.sel)
		exp.log2 <- log2(exp.sel +1)
		r.max <- apply(exp.log2,1,max)
		r.sd  <- apply(exp.log2,1,sd)
		exp.sel <- exp.log2[r.max > 8 & r.sd > 1,]
		p.anns.sel <- probe.anns[rownames(exp.sel),]
		p.anns.uniq <- p.anns.sel[!duplicated(p.anns.sel$symbol),]
		exp.uniq <- exp.sel[rownames(p.anns.uniq),]
		rownames(exp.uniq) <- p.anns.uniq$symbol
		true.label <- factor(eset.sel$subtype)
		all.dat <- data.frame(t(exp.uniq) ,tumor.type = eset.sel$subtype)
		flds <- createFolds(1:nrow(all.dat), k = K, list = TRUE, returnTrain = FALSE)
		RF.predicted.label = character(length = nrow(all.dat))
		cl.idx <- which(colnames(all.dat) == "tumor.type")
		for(n in 1:K)
		{
      			test.idx <- flds[[n]]
        		train.all <- all.dat[-test.idx,]
        		features <- select_features_for_training(train.all,num = 100)
        		test.all <- all.dat[test.idx,-c(cl.idx)]
        		train <- train.all[,c(features,"tumor.type")]
        		test  <- test.all[,features]
        		model <- randomForest(tumor.type  ~ . , data = train,ntree= 1000,mtry=31)
        		RF.pred <- predict(model, newdata = test)
        		RF.predicted.label[test.idx] <- as.character(RF.pred)
		}
		RF.predicted = factor(RF.predicted.label)
		RF.res[[i]] <- confusionMatrix(RF.predicted,true.label)
        	model.list[[i]] <- randomForest(tumor.type  ~ . , data = all.dat,ntree= 1000,mtry=31)
	}
}

save(RF.res,model.list,file="models/subtype/subtype_predictor_results.rdata")

#write.csv(RF.res[[2]],"contingency_table_RF.csv")
#write.csv(RF.res[[4]],"preformance_metrics_RF.csv")

