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

subtypes = table(sample.anns[sample.anns$tumor.type == "BRCA","subtype"])
subtype.names <- setdiff(names(subtypes[subtypes > 25]),"")
eset.sel <- eset[,eset$tumor.type == "BRCA" & eset$subtype %in% subtype.names]
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
breast.model <- randomForest(tumor.type  ~ . , data = all.dat ,ntree= 1000,mtry=31)


load("data/external/metabric_combined.rdata")
#met.eset <- metabric.combined.eset[, metabric.combined.eset$Pam50Subtype %in% c("Basal","Her2","LumA","LumB")]
met.eset <- metabric.combined.eset
met.exp <- Biobase::exprs(met.eset)
met.p.anns <- featureData(met.eset)@data
met.sample.anns <- pData(met.eset)
r.mu <- apply(met.exp,1,mean)
o = order(r.mu,decreasing=T)
met.exp.ord <- met.exp[o,]
met.p.anns.ord <- met.p.anns[o,]
idx <- duplicated(met.p.anns.ord$Symbol)
met.uniq <- met.exp.ord[! idx,]
met.p.anns.uniq <- met.p.anns.ord[!idx,]
rownames(met.uniq) <- met.p.anns.uniq$Symbol
features.common <- intersect(features,rownames(met.uniq))

exp.t <- t(exp.uniq[features.common,])
exp.scaled = scale(exp.t)
all.dat <- data.frame(exp.scaled ,tumor.type = eset.sel$subtype)
breast.model <- randomForest(tumor.type  ~ . , data = all.dat ,ntree= 1000,mtry=31)

met.t <- t(met.uniq[features.common,])
met.scaled = scale(met.t)
met.pred <- predict(breast.model, newdata = met.scaled)

met.true.labels <- met.sample.anns$Pam50Subtype
idx <- met.true.labels %in% c("Basal","Her2","LumA","LumB")
met.true.labels[!idx] <- NA
met.true.labels[met.sample.anns$Pam50Subtype == "Basal"] <- "Basal-like"
met.true.labels[met.sample.anns$Pam50Subtype == "Her2"] <- "Her2 enriched"
met.true.labels[met.sample.anns$Pam50Subtype == "LumA"] <- "Luminal A"
met.true.labels[met.sample.anns$Pam50Subtype == "LumB"] <- "Luminal B"
lev <- c("Basal-like", "Her2 enriched", "Luminal A","Luminal B")
met.pred <- factor(met.pred,levels = lev)
met.true <- factor(met.true.labels ,levels = lev)

breast.cm <- confusionMatrix(met.pred,met.true)

save(breast.cm, file = "models/subtype_ext_val/subtype_classifier_results_from_metabric_samples.rdata")

