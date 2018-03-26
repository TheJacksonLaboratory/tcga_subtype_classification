BRAF.cases <- scan("data/raw/cases_BRAF.txt",what="character")
BRAF.ids <- BRAF.cases[grep("TCGA",BRAF.cases)]
BRAF.dat <- data.frame(sampleid=BRAF.ids, subtype="BRAF",cancertype="Cutaneous melanoma")

RAS.cases <- scan("data/raw/cases_RAS.txt",what="character")
RAS.ids <- RAS.cases[grep("TCGA", RAS.cases)]
RAS.dat <- data.frame(sampleid=RAS.ids, subtype="RAS",cancertype="Cutaneous melanoma")

NF1.cases <- scan("data/raw/cases_NF1.txt",what="character")
NF1.ids <- NF1.cases[grep("TCGA", NF1.cases)]
NF1.dat <- data.frame(sampleid=NF1.ids, subtype="NF1",cancertype="Cutaneous melanoma")

TWT.cases <- scan("data/raw/cases_Triple_WT.txt",what="character")
TWT.ids <- TWT.cases[grep("TCGA", TWT.cases)]
TWT.dat <- data.frame(sampleid=TWT.ids, subtype="Triple wild-type",cancertype="Cutaneous melanoma")

all_skcm.dat <- rbind(BRAF.dat, RAS.dat, NF1.dat, TWT.dat)

all_skcm.mat <- as.matrix(all_skcm.dat)
all_skcm.TCGAids <- all_skcm.mat[grep("TCGA", all_skcm.mat)]
skcm.sampleid.mat2 <- substr(all_skcm.TCGAids, 1, 15)
skcm.sampleid.mat3 <- substr(all_skcm.TCGAids, 1, 12)
skcm.sampleid.dat <- as.data.frame(skcm.sampleid.mat2, skcm.sampleid.mat3)

rownames(skcm.sampleid.dat) -> skcm.sampleid.dat$sampleid
rownames(skcm.sampleid.dat)=NULL
colnames(skcm.sampleid.dat) <- c("sampleid", "patientid")
all_skcm.dat <- merge(all_skcm.dat, skcm.sampleid.dat, by="sampleid", all.x=FALSE)
write.table(all_skcm.dat, file="data/interim/all_skcm.dat", sep="\t")
