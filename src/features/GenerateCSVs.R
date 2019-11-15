print('External Validation')
library(xtable)
library(Biobase)
library(class)
library(randomForest)
library(xtable)
library(party)
library(caret)
library(e1071)
library(gbm)

set.seed(1024)
load("../DataPool/combined_data.rdata")
eset.1 <- combined.eset[,combined.eset$sample.type %in% c("01","03")]
#eset.2 <- combined.eset[,combined.eset$tumor.type =="SKCM" & combined.eset$sample.type == "06"]
eset <- Biobase::combine(eset.1)

tumor.type <- unique(eset$tumor.type)
exp <- Biobase::exprs(eset)
exp.log2 <- log2(exp +1)
exp.sel <- exp.log2
all.dat <- data.frame(t(exp.sel) ,tumor.type = eset$tumor.type)
features.common <- colnames(all.dat)

train.dat <- all.dat
train.scaled <- train.dat
tp <- as.character(train.dat$tumor.type)
tp[tp %in% c('COAD','READ')] <- 'COADREAD'
all.sel.dat <-  data.frame(train.scaled)
all.sel.dat$tumor.type <- tp
write.csv(all.sel.dat, file='../DataPool/tcga_training.csv')

load("../DataPool/pdx_eset.rdata")
pdx.anns <- featureData(pdx.eset)@data
pdx.entid <- paste("X",pdx.anns$ENTREZ_ID,sep="")

pdx.anns.sel <- pdx.anns[pdx.entid %in% features.common,]
pdx.exp <- Biobase::exprs(pdx.eset[rownames(pdx.anns.sel),])
r.mu <- apply(pdx.exp,1,mean)
o <- order(r.mu,decreasing=T)
pdx.ord <- pdx.exp[o,]
pdx.anns.ord <- pdx.anns.sel[o,]
pdx.uniq <- pdx.anns.ord[!duplicated(pdx.anns.ord$ENTREZ_ID),]
pdx.exp.uniq <- pdx.ord[rownames(pdx.uniq),]
rownames(pdx.exp.uniq) <- paste("X",pdx.uniq$ENTREZ_ID,sep="")
test.dat <- t(log2(pdx.exp.uniq+1))

pdx.true.labels = as.character(pdx.eset$primary_site)
pdx.true.labels[pdx.true.labels == "ampulla"] <- NA
pdx.true.labels[pdx.true.labels == "ampulla of Vater"] <- NA
pdx.true.labels[pdx.true.labels == "Bone marrow"] <- "LAML"
pdx.true.labels[grep("B-[c/C]ell",pdx.eset$clinical_diagnosis)] <- "DLBC"
pdx.true.labels[pdx.true.labels == "Bladder"] <- "BLCA"
pdx.true.labels[pdx.true.labels == "Brain"] <- "LGG/GBM"
pdx.true.labels[pdx.true.labels == "Breast"] <- "BRCA"
pdx.true.labels[pdx.true.labels == "Colon"] <- "COADREAD"
pdx.true.labels[pdx.true.labels == "Kidney"] <- "KIRC/KIRP/KIRH"
pdx.true.labels[pdx.true.labels == "Lung"] <- "LUAD/LUSC"
pdx.true.labels[pdx.true.labels == "Ovary"] <- "OV"
pdx.true.labels[pdx.true.labels == "Pancreas"] <- "PAAD"
pdx.true.labels[pdx.true.labels == "Rectum"] <- "COADREAD"
pdx.true.labels[pdx.true.labels == "Skin"] <- "SKCM"
pdx.true.labels[pdx.true.labels == "SoftTissue"] <- "SARC"
pdx.true.labels[pdx.true.labels == "Soft tissue"] <- "SARC"
pdx.true.labels[pdx.true.labels == "Bone"] <- "SARC"

output <- cbind(test.dat, tumor.type=pdx.true.labels)
write.csv(output, file='../DataPool/ExternalDataPdx.csv')

#################clinical data
eset.backup <- eset #in case eset is overwrite by clinical data object
load("../DataPool/clinical_lab_samples.rdata")
clinical.eset <- eset
eset <- eset.backup
clinical.anns <- rownames(featureData(clinical.eset)@data)
clinical.entid <- paste("X",gsub(".*\\|","",clinical.anns),sep="")
clinical.anns.sel <- clinical.anns[clinical.entid %in% features.common]
clinical.exp <- Biobase::exprs(clinical.eset[clinical.anns.sel,])
r.mu <- apply(clinical.exp,1,mean)
o <- order(r.mu,decreasing=T)
clinical.ord <- clinical.exp[o,]
clinical.anns.ord <- clinical.anns.sel[o]
clinical.uniq <- clinical.anns.ord[!duplicated(clinical.anns.ord)]
clinical.exp.uniq <- clinical.ord[clinical.uniq,]
rownames(clinical.exp.uniq) <- paste("X",gsub(".*\\|","",clinical.uniq),sep="")
test.dat <- t(log2(clinical.exp.uniq+1))
test.dat[is.na(test.dat)] <- 0

true_labels <- read.csv("../DataPool/clinical_true_labels.csv", stringsAsFactors=FALSE)
clinical.true.labels <- true_labels$true.label
names(clinical.true.labels) <- true_labels$sampleid
clinical.true.labels <- clinical.true.labels[rownames(test.dat)]
clinical.true.labels[clinical.true.labels %in% c('COAD','READ')] <- 'COADREAD'
output <- cbind(test.dat, tumor.type=clinical.true.labels)
write.csv(output, file='../DataPool/ExternalDataClinical.csv')

######################tcga metastatic
eset.tcga.meta <- combined.eset[,combined.eset$sample.type == "06"]
tumor.type.meta <- unique(eset.tcga.meta$tumor.type)
exp.meta <- Biobase::exprs(eset.tcga.meta)
exp.log2.meta <- log2(exp.meta +1)
exp.sel.meta <- exp.log2.meta
all.dat.meta <- data.frame(t(exp.sel.meta),tumor.type = eset.tcga.meta$tumor.type)

all.dat.meta <- all.dat.meta[,colnames(all.dat.meta) %in% features.common]
test.dat <- all.dat.meta
test.dat[is.na(test.dat)] <- 0
meta.tcga.true.labels <- eset.tcga.meta$tumor.type
meta.tcga.true.labels[meta.tcga.true.labels %in% c('COAD','READ')] <- 'COADREAD'

output <- cbind(test.dat, tumor.type=meta.tcga.true.labels)
write.csv(output, file='../DataPool/ExternalDataMeta.csv')

load("../DataPool/features_10.rdata")
load("../DataPool/stad_coad_model.rdata")
write.csv(c(coad_stad_features, features), file='../DataPool/randomForestColumns_254.csv')

