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

subtypes = table(sample.anns[sample.anns$tumor.type == "OV","subtype"])
subtype.names <- setdiff(names(subtypes[subtypes > 25]),"")
eset.sel <- eset[,eset$tumor.type == "OV" & eset$subtype %in% subtype.names]
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
RF.res <- confusionMatrix(RF.predicted,true.label)

features.all <- select_features_for_training(all.dat,num = 200)
features <- intersect(features.all,rownames(exp.uniq))
exp.sel <- exp.uniq[features,]
all.dat <- data.frame(t(exp.sel) ,tumor.type = eset.sel$subtype)
ovarian.model <- randomForest(tumor.type  ~ . , data = all.dat ,ntree= 1000,mtry=31)


load("data/external/rma_normalized_data_nov_2012.rdata")
library(hgu133plus2.db)
aocs.exp <- Biobase::exprs(eset)
aocs.sample.anns <- pData(eset)
gene.symbol <- as.character(mget(rownames(aocs.exp),hgu133plus2SYMBOL,ifnotfound=NA))
aocs.p.anns <- data.frame(probeid = rownames(aocs.exp), symbol = gene.symbol,stringsAsFactors=F)
r.mu <- apply(aocs.exp,1,mean)
o = order(r.mu,decreasing=T)
aocs.exp.ord <- aocs.exp[o,]
aocs.p.anns.ord <- aocs.p.anns[o,]
idx <- duplicated(aocs.p.anns.ord$symbol)
aocs.uniq <- aocs.exp.ord[! idx,]
aocs.p.anns.uniq <- aocs.p.anns.ord[!idx,]
rownames(aocs.uniq) <- aocs.p.anns.uniq$symbol
features.common <- intersect(features,rownames(aocs.uniq))

exp.t <- t(exp.uniq[features.common,])
exp.scaled = scale(exp.t)
all.dat <- data.frame(exp.scaled ,tumor.type = eset.sel$subtype)
ovarian.model <- randomForest(tumor.type  ~ . , data = all.dat ,ntree= 1000,mtry=31)

aocs.t <- t(aocs.uniq[features.common,])
aocs.scaled = scale(aocs.t)
aocs.pred <- predict(ovarian.model, newdata = aocs.scaled)

aocs.true.labels <- aocs.sample.anns$MolSubtype
idx <- aocs.true.labels %in% c(1,2,4,5)
aocs.true.labels[!idx] <- NA
aocs.true.labels[aocs.sample.anns$MolSubtype == 1] <- "Mesenchymal"
aocs.true.labels[aocs.sample.anns$MolSubtype == 2] <- "Immunoreactive"
aocs.true.labels[aocs.sample.anns$MolSubtype == 4] <- "Differentiated"
aocs.true.labels[aocs.sample.anns$MolSubtype == 5] <- "Proliferative"

lev <- c("Mesenchymal", "Immunoreactive", "Differentiated","Proliferative")
aocs.pred <- factor(aocs.pred,levels = lev)
aocs.true <- factor(aocs.true.labels ,levels = lev)

ovarian.cm <- confusionMatrix(aocs.pred,aocs.true)

save(ovarian.cm, file = "models/subtype_ext_val/subtype_classifier_results_from_AOCS_samples.rdata")

