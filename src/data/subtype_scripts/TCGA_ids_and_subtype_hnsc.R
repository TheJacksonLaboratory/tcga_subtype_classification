atyp.cases <-scan("data/raw/cases_atypical.txt",what="character")
atyp.ids <-atyp.cases[grep("TCGA", atyp.cases)]
atyp.dat <- data.frame(patientid=atyp.ids, subtype="Atypical",cancertype="Head and Neck")

basal_hnsc.cases <- scan("data/raw/hnsc_cases_basal.txt",what="character")
basal_hnsc.ids <- basal_hnsc.cases[grep("TCGA", basal_hnsc.cases)]
basal_hnsc.dat <- data.frame(patientid=basal_hnsc.ids, subtype="Basal",cancertype="Head and Neck")

class_hnsc.cases <- scan("data/raw/hnsc_cases_classical.txt", what="character")
class_hnsc.ids <- class_hnsc.cases[grep("TCGA", class_hnsc.cases)]
class_hnsc.dat <- data.frame(patientid=class_hnsc.ids, subtype="Classical",cancertype="Head and Neck")

mesen_hnsc.cases <- scan("data/raw/hnsc_cases_mesenchymal.txt", what="character")
mesen_hnsc.ids <- mesen_hnsc.cases[grep("TCGA", mesen_hnsc.cases)]
mesen_hnsc.dat <- data.frame(patientid=mesen_hnsc.ids, subtype="Mesenchymal",cancertype="Head and Neck")

all_hnsc <- rbind(atyp.dat, basal_hnsc.dat, class_hnsc.dat, mesen_hnsc.dat)
all_hnsc$patientid <- gsub("\\.", "-", all_hnsc$patientid)
all_hnsc[order(all_hnsc$patientid), ]

hnsc_sampleid.cases <- scan("data/raw/hnsc_data_clinical_sample.txt",what="character")
hnsc_sampleid.ids <- hnsc_sampleid.cases[grep("TCGA", hnsc_sampleid.cases)]
hnsc_sampleid.mat <- matrix(hnsc_sampleid.ids)
hnsc_sampleid.mat <- hnsc_sampleid.mat[order(hnsc_sampleid.mat), ]
n <- length(hnsc_sampleid.mat)
hnsc_sampleid2 <- hnsc_sampleid.mat[seq(1,n,2)]
hnsc_sampleid2.mat <- matrix(hnsc_sampleid2)
hnsc_sampleid3.mat <- substr(hnsc_sampleid2.mat, 1, 12)
hnsc_sampleid.dat <- as.data.frame(hnsc_sampleid2.mat, hnsc_sampleid3.mat)
rownames(hnsc_sampleid.dat) -> hnsc_sampleid.dat$patientid
rownames(hnsc_sampleid.dat)=NULL
colnames(hnsc_sampleid.dat) <-c("sampleid", "patientid")
hnsc_sampleid.dat[order(hnsc_sampleid.dat$patientid), ]


all_hnsc.dat <- merge(all_hnsc, hnsc_sampleid.dat, by="patientid", all.x=FALSE)
write.table(all_hnsc.dat, file="data/interim/all_hnsc.dat", sep="\t")
