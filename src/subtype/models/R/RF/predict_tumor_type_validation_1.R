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

tumor.type <- unique(eset$tumor.type)
exp <- Biobase::exprs(eset)
exp.log2 <- log2(exp +1)
r.max <- apply(exp.log2,1,max)
r.sd  <- apply(exp.log2,1,sd)
exp.sel <- exp.log2[r.max > 8 & r.sd > 1,]
all.dat <- data.frame(t(exp.sel) ,tumor.type = eset$tumor.type)

features <- select_features_for_training(all.dat,num = 100)
load("data/external/gse2109.rdata")
gse.p21.anns <- featureData(gse2109.eset)@data
gse.p21.entid <- paste("X",gse.p21.anns$ENTREZ_GENE_ID,sep="")
load("data/external/gse18549.rdata")
gse.p18.anns <- featureData(gse18549.eset)@data
gse.p18.entid <- paste("X",gse.p18.anns$ENTREZ_GENE_ID,sep="")
load("data/external/gse19750.rdata")

load("data/external/pdx_eset.rdata")
pdx.anns <- featureData(pdx.eset)@data
pdx.entid <- paste("X",pdx.anns$ENTREZ_ID,sep="")

features.common = intersect(gse.p18.entid,features)
features.common = intersect(gse.p21.entid,features.common)
features.common = intersect(pdx.entid,features.common)

gse.p21.anns.sel <- gse.p21.anns[gse.p21.entid %in% features.common,]
gse.21.exp <- Biobase::exprs(gse2109.eset[rownames(gse.p21.anns.sel),])
r.mu <- apply(gse.21.exp,1,mean)
o <- order(r.mu,decreasing=T)
gse.21.ord <- gse.21.exp[o,]
gse.p21.anns.ord <- gse.p21.anns.sel[o,]
gse.21.uniq <- gse.p21.anns.ord[!duplicated(gse.p21.anns.ord$ENTREZ_GENE_ID),]
gse.21.exp.uniq <- gse.21.ord[rownames(gse.21.uniq),]
rownames(gse.21.exp.uniq) <- paste("X",gse.21.uniq$ENTREZ_GENE_ID,sep="")

gse.p18.anns.sel <- gse.p18.anns[gse.p18.entid %in% features.common,]
gse.18.exp <- Biobase::exprs(gse18549.eset[rownames(gse.p18.anns.sel),])
r.mu <- apply(gse.18.exp,1,mean)
o <- order(r.mu,decreasing=T)
gse.18.ord <- gse.18.exp[o,]
gse.p18.anns.ord <- gse.p18.anns.sel[o,]
gse.18.uniq <- gse.p18.anns.ord[!duplicated(gse.p18.anns.ord$ENTREZ_GENE_ID),]
gse.18.exp.uniq <- gse.18.ord[rownames(gse.18.uniq),]
rownames(gse.18.exp.uniq) <- paste("X",gse.18.uniq$ENTREZ_GENE_ID,sep="")


gse18 <- gse.18.exp.uniq
gse21 <- log2(gse.21.exp.uniq[rownames(gse18),])

train.dat <- all.dat[,features.common]
train.scaled <- scale(train.dat)
all.sel.dat <-  data.frame(train.scaled, tumor.type = eset$tumor.type)
final.model <- randomForest(tumor.type  ~ . , data = all.sel.dat,ntree= 1000,mtry=31)

test.dat <- scale(t(gse21))
gse2109.predicted.type <- predict(final.model, newdata = test.dat)
true.type = as.character(gse2109.eset$source_name_ch1)
redifined.type = true.type
redifined.type[] <- "NotSure"
redifined.type[grep("vary",true.type)] <- "OV"
redifined.type[grep("mentum",true.type)] <- "OV"
redifined.type[grep("allopian",true.type)] <- "OV"
redifined.type[grep("ladder",true.type)] <- "BLCA"
redifined.type[grep("reast",true.type)] <- "BRCA"
redifined.type[grep("olon",true.type)] <- "COAD"
redifined.type[grep("idney",true.type)] <- "KIRC/KIRP/KIRH"
redifined.type[grep("rain",true.type)] <- "LGG/GBM"
redifined.type[grep("iver",true.type)] <- "LIHC"
redifined.type[grep("ung",true.type)] <- "LUAD/LUSC"
redifined.type[grep("rostate",true.type)] <- "PRAD"
redifined.type[grep("Testic",true.type)] <- "TGCT"
redifined.type[grep("Thyroid",true.type)] <- "THCA"

gse2109.predicted.type = as.character(gse2109.predicted.type)
gse2109.predicted.type[gse2109.predicted.type %in% c("KIRP","KIRC","KIRH")] <- "KIRC/KIRP/KIRH"
gse2109.predicted.type[gse2109.predicted.type %in% c("LUSC","LUAD")] <- "LUAD/LUSC"
gse2109.predicted.type[gse2109.predicted.type %in% c("LGG","GBM")] <- "LGG/GBM"

res <- data.frame(predicted = gse2109.predicted.type, true.type = redifined.type, stringsAsFactors=FALSE)
res.sel <- res[res$true.type != "NotSure",]
lev = unique(res.sel$predicted)
gse2109.predicted.type  <- factor(res.sel$predicted , levels = lev)
gse2109.true.type  <- factor(res.sel$true.type, levels=lev)
gse2109.cm  <-  confusionMatrix(gse2109.predicted.type,gse2109.true.type)

test.dat <- scale(t(gse18))
gse18549.predicted <- as.character(predict(final.model, newdata = test.dat))
gse18549.true <- gse18549.eset$primary
gse18549.true[] <- NA
gse18549.true[grep("Ovary", gse18549.eset$primary)] <- "OV"
gse18549.true[grep("Ovarian", gse18549.eset$primary)] <- "OV"
gse18549.true[grep("Breast", gse18549.eset$primary)] <- "BRCA"
gse18549.true[grep("Prostate", gse18549.eset$primary)] <- "PRAD"
gse18549.true[grep("Renal", gse18549.eset$primary)] <- "KIRC/KIRP/KIRH"
gse18549.true[grep("Lung", gse18549.eset$primary)] <- "LUAD/LUSC"
gse18549.true[grep("Colon", gse18549.eset$primary)] <- "COAD"
gse18549.true[grep("Rectal", gse18549.eset$primary)] <- "READ"

gse18549.predicted[gse18549.predicted %in% c("KICH","KIRC","KIRP")] <- "KIRC/KIRP/KIRH"
gse18549.predicted[gse18549.predicted %in% c("LUSC","LUAD") ] <- "LUAD/LUSC"
x <- data.frame(predicted.label = gse18549.predicted, true.label= gse18549.true,stringsAsFactors=FALSE)
x.c <- x[!is.na(x$true.label),]
lev.all <- unique(gse18549.predicted)
x.all <- factor(x.c$predicted.label, levels = lev.all)
y.all <- factor(x.c$true.label, levels=lev.all)
gse18549.cm <-  confusionMatrix(x.all,y.all)

pdx.anns.sel <- pdx.anns[pdx.entid %in% features.common,]
pdx.exp <- Biobase::exprs(pdx.eset[rownames(pdx.anns.sel),])
r.mu <- apply(pdx.exp,1,mean)
o <- order(r.mu,decreasing=T)
pdx.ord <- pdx.exp[o,]
pdx.anns.ord <- pdx.anns.sel[o,]
pdx.uniq <- pdx.anns.ord[!duplicated(pdx.anns.ord$ENTREZ_ID),]
pdx.exp.uniq <- pdx.ord[rownames(pdx.uniq),]
rownames(pdx.exp.uniq) <- paste("X",pdx.uniq$ENTREZ_ID,sep="")
test.dat <- scale(t(log2(pdx.exp.uniq+1)))
pdx.predicted <- as.character(predict(final.model, newdata = test.dat))
#pdx.predicted[pdx.predicted %in% c("LAML","DLBC")] <- "LAML/DLBC"
pdx.predicted[pdx.predicted %in% c("LGG","GBM")] <- "LGG/GBM"
pdx.predicted[pdx.predicted %in% c("KIRC","KIRP","KIRH")] <- "KIRC/KIRP/KIRH"
pdx.predicted[pdx.predicted %in% c("LUAD","LUSC")] <- "LUAD/LUSC"

pdx.true.labels = as.character(pdx.eset$primary_site)
pdx.true.labels[pdx.true.labels == "ampulla"] <- NA
pdx.true.labels[pdx.true.labels == "Bone marrow"] <- "LAML"
pdx.true.labels[grep("B-[c/C]ell",pdx.eset$clinical_diagnosis)] <- "DLBC"
pdx.true.labels[pdx.true.labels == "Bladder"] <- "BLCA"
pdx.true.labels[pdx.true.labels == "Brain"] <- "LGG/GBM"
pdx.true.labels[pdx.true.labels == "Breast"] <- "BRCA"
pdx.true.labels[pdx.true.labels == "Colon"] <- "COAD"
pdx.true.labels[pdx.true.labels == "Kidney"] <- "KIRC/KIRP/KIRH"
pdx.true.labels[pdx.true.labels == "Lung"] <- "LUAD/LUSC"
pdx.true.labels[pdx.true.labels == "Ovary"] <- "OV"
pdx.true.labels[pdx.true.labels == "Pancreas"] <- "PAAD"
pdx.true.labels[pdx.true.labels == "Rectum"] <- "READ"
pdx.true.labels[pdx.true.labels == "Skin"] <- "SKCM"
pdx.true.labels[pdx.true.labels == "SoftTissue"] <- "SARC"
pdx.true.labels[pdx.true.labels == "Bone"] <- "SARC"

pdx.levels  <- unique(pdx.predicted)
pdx.predicted  <- factor(pdx.predicted, levels = pdx.levels)
pdx.true  <- factor(pdx.true.labels, levels=pdx.levels)
pdx.cm <- confusionMatrix(pdx.predicted,pdx.true)


save(pdx.cm, gse18549.cm, gse2109.cm, file ="models/primary_type_ext_val.rdata")

