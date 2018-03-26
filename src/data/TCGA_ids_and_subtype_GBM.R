pron.cases <- scan("data/raw/cases_proneural.txt",what="character")
pron.ids <- pron.cases[grep("TCGA", pron.cases)]
pron.dat <- data.frame(sampleid = pron.ids, subtype="Proneural",cancertype="Glioblastoma")

neu.cases <- scan("data/raw/cases_neural.txt",what="character")
neu.ids <- neu.cases[grep("TCGA", neu.cases)]
neu.dat <- data.frame(sampleid = neu.ids, subtype="Neural",cancertype="Glioblastoma")

clas.cases <- scan("data/raw/gbm_cases_classical.txt",what="character")
clas.ids <- clas.cases[grep("TCGA", clas.cases)]
clas.dat <- data.frame(sampleid = clas.ids, subtype="Classical",cancertype="Glioblastoma")

mesen.cases <- scan("data/raw/gbm_cases_mesenchymal.txt",what="character")
mesen.ids <- mesen.cases[grep("TCGA", mesen.cases)]
mesen.dat <- data.frame(sampleid=mesen.ids, subtype="Mesenchymal",cancertype="Glioblastoma")

all_GBM.dat <-rbind(pron.dat, neu.dat, clas.dat, mesen.dat)

all_GBM.mat <- as.matrix(all_GBM.dat)
all_GBM.TCGAids <- all_GBM.mat[grep("TCGA", all_GBM.mat)]
GBM.sampleid.mat2 <- substr(all_GBM.TCGAids, 1, 15)
GBM.sampleid.mat3 <- substr(all_GBM.TCGAids, 1, 12)
GBM.sampleid.dat <- as.data.frame(GBM.sampleid.mat2, GBM.sampleid.mat3)

rownames(GBM.sampleid.dat) -> GBM.sampleid.dat$sampleid
rownames(GBM.sampleid.dat)=NULL
colnames(GBM.sampleid.dat) <- c("sampleid", "patientid")
all_GBM.dat <- merge(all_GBM.dat, GBM.sampleid.dat, by="sampleid", all.x=FALSE)
 write.table(all_GBM.dat, file="data/interim/all_GBM.dat", sep="\t")
