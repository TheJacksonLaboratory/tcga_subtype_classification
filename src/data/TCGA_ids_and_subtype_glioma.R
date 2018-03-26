lgr1.cases <-scan("data/raw/cases_LGr1.txt",what="character")
lgr1.ids <- lgr1.cases[grep("TCGA", lgr1.cases)]
lgr1.dat <- data.frame(patientid=lgr1.ids, subtype="most upregulate TERT", cancertype="Glioma")

lgr2.cases <- scan("data/raw/cases_LGr2.txt",what="character")
lgr2.ids <- lgr2.cases[grep("TCGA", lgr2.cases)]
lgr2.dat <- data.frame(patientid=lgr2.ids, subtype="some upregulate TERT", cancertype="Glioma")

lgr3.cases <- scan("data/raw/cases_LGr3.txt",what="character")
lgr3.ids <- lgr3.cases[grep("TCGA", lgr3.cases)]
lgr3.dat <- data.frame(patientid=lgr3.ids, subtype="TP53 mutations", cancertype="Glioma")

lgr4.cases <- scan("data/raw/cases_LGr4.txt", what="character")
lgr4.ids <- lgr4.cases[grep("TCGA", lgr4.cases)]
lgr4.dat <- data.frame(patientid=lgr4.ids, subtype="IDH wild-type", cancertype="Glioma")

all_glioma <- rbind(lgr1.dat, lgr2.dat, lgr3.dat, lgr4.dat)
all_glioma[order(all_glioma$patientid), ]


glioma.sampleid.cases <- scan("data/raw/glioma_gdc_manifest.2018-01-31T19_50_13.151970.txt",what="character")
glioma.sampleid.ids <- glioma.sampleid.cases[grep("TCGA-.*gdc_realn", glioma.sampleid.cases)]
glioma.sampleid.mat <-matrix(glioma.sampleid.ids)
glioma.sampleid.mat2 <- substr(glioma.sampleid.mat, 6, 20)
glioma.sampleid.mat3 <- substr(glioma.sampleid.mat, 6, 17)
glioma.sampleid.dat <- as.data.frame(glioma.sampleid.mat2, glioma.sampleid.mat3)

rownames(glioma.sampleid.dat) -> glioma.sampleid.dat$sampleid
rownames(glioma.sampleid.dat)=NULL
colnames(glioma.sampleid.dat) <- c("sampleid", "patientid")
all_glioma.dat <- merge(all_glioma, glioma.sampleid.dat, by="patientid", all.x=FALSE)
write.table(all_glioma.dat, file="data/interim/all_glioma.dat", sep="\t")
