POLE_ucec.cases <- scan("data/raw/cases_pole.txt", what="character")
POLE_ucec.ids <- POLE_ucec.cases[grep("TCGA", POLE_ucec.cases)]
POLE_ucec.dat <- data.frame(sampleid=POLE_ucec.ids, subtype="POLE", cancertype="Endometrial")

MSI_ucec.cases <- scan("data/raw/ucec_cases_msi.txt", what="character")
MSI_ucec.ids <- MSI_ucec.cases[grep("TCGA", MSI_ucec.cases)]
MSI_ucec.dat <- data.frame(sampleid=MSI_ucec.ids, subtype="Microsatellite instability hypermutated",cancertype="Endometrial")

CN_low_ucec.cases <- scan("data/raw/cases_cnlow.txt",what="character")
CN_low_ucec.ids <- CN_low_ucec.cases[grep("TCGA", CN_low_ucec.cases)]
CN_low_ucec.dat <- data.frame(sampleid=CN_low_ucec.ids, subtype="Copy number low", cancertype="Endometrial")

CN_high_ucec.cases <- scan("data/raw/cases_cnhigh.txt",what="character")
CN_high_ucec.ids <- CN_high_ucec.cases[grep("TCGA", CN_high_ucec.cases)]
CN_high_ucec.dat <- data.frame(sampleid=CN_high_ucec.ids, subtype="Copy number high", cancertype="Endometrial")

all_ucec.dat <-rbind(POLE_ucec.dat, MSI_ucec.dat,CN_low_ucec.dat,CN_high_ucec.dat)


all_ucec.mat <- as.matrix(all_ucec.dat)
all_ucec.TCGAids <- all_ucec.mat[grep("TCGA", all_ucec.mat)]
ucec.sampleid.mat2 <- substr(all_ucec.TCGAids, 1, 15)
ucec.sampleid.mat3 <- substr(all_ucec.TCGAids, 1, 12)
ucec.sampleid.dat <- as.data.frame(ucec.sampleid.mat2, ucec.sampleid.mat3)

rownames(ucec.sampleid.dat) -> ucec.sampleid.dat$sampleid
rownames(ucec.sampleid.dat)=NULL
colnames(ucec.sampleid.dat) <- c("sampleid", "patientid")
all_ucec.dat <- merge(all_ucec.dat, ucec.sampleid.dat, by="sampleid", all.x=FALSE)
write.table(all_ucec.dat, file="data/interim/all_ucec.dat", sep="\t")
