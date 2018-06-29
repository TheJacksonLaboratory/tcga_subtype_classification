CIN.cases <-scan("data/raw/cases_cin.txt",what="character")
CIN.ids <- CIN.cases[grep("TCGA", CIN.cases)]
CIN.dat <- data.frame(sampleid=CIN.ids, subtype="Chromosomal instability",cancertype="Gastric")

EBV.cases <- scan("data/raw/cases_ebv.txt",what="character")
EBV.ids <- EBV.cases[grep("TCGA", EBV.cases)]
EBV.dat <- data.frame(sampleid=EBV.ids, subtype="Epstein-Barr virus positive",cancertype="Gastric")

GS.cases <- scan("data/raw/cases_gs.txt",what="character")
GS.ids <- GS.cases[grep("TCGA",GS.cases)]
GS.dat <- data.frame(sampleid=GS.ids, subtype="Genomically stable",cancertype="Gastric")

MSI.cases <- scan("data/raw/stad_cases_msi.txt",what="character")
MSI.ids <- MSI.cases[grep("TCGA", MSI.cases)]
MSI.dat <- data.frame(sampleid=MSI.ids, subtype="Microsatellite instability",cancertype="Gastric")

all_stad.dat <- rbind(CIN.dat, EBV.dat, GS.dat, MSI.dat)


all_stad.mat <- as.matrix(all_stad.dat)
all_stad.TCGAids <- all_stad.mat[grep("TCGA", all_stad.mat)]
stad.sampleid.mat2 <- substr(all_stad.TCGAids, 1, 15)
stad.sampleid.mat3 <- substr(all_stad.TCGAids, 1, 12)
stad.sampleid.dat <- as.data.frame(stad.sampleid.mat2, stad.sampleid.mat3)

rownames(stad.sampleid.dat) -> stad.sampleid.dat$sampleid
rownames(stad.sampleid.dat)=NULL
colnames(stad.sampleid.dat) <- c("sampleid", "patientid")
all_stad.dat <- merge(all_stad.dat, stad.sampleid.dat, by="sampleid", all.x=FALSE)
write.table(all_stad.dat, file="data/interim/all_stad.dat", sep="\t")
