bron.cases <- scan("data/raw/cases_bronchioid.txt",what="character")
bron.ids <- bron.cases[grep("TCGA", bron.cases)]
bron.dat <- data.frame(sampleid=bron.ids, subtype="Terminal respiratory unit",cancertype="Lung adenocarcinoma")

mag.cases <- scan("data/raw/cases_magnoid.txt", what="character")
mag.ids <- mag.cases[grep("TCGA", mag.cases)]
mag.dat <- data.frame(sampleid=mag.ids, subtype="Proximal-proliferative",cancertype="Lung adenocarcinoma")

squ.cases <- scan("data/raw/cases_squamoid.txt", what="character")
squ.ids <- squ.cases[grep("TCGA", squ.cases)]
squ.dat <- data.frame(sampleid=squ.ids, subtype="Proximal-proliferative",cancertype="Lung adenocarcinoma")

all_luad.dat <- rbind(bron.dat, mag.dat, squ.dat)


all_luad.mat <- as.matrix(all_luad.dat)
all_luad.TCGAids <- all_luad.mat[grep("TCGA", all_luad.mat)]
luad.sampleid.mat2 <- substr(all_luad.TCGAids, 1, 15)
luad.sampleid.mat3 <- substr(all_luad.TCGAids, 1, 12)
luad.sampleid.dat <- as.data.frame(luad.sampleid.mat2, luad.sampleid.mat3)

rownames(luad.sampleid.dat) -> luad.sampleid.dat$sampleid
rownames(luad.sampleid.dat)=NULL
colnames(luad.sampleid.dat) <- c("sampleid", "patientid")
all_luad.dat <- merge(all_luad.dat, luad.sampleid.dat, by="sampleid", all.x=FALSE)
write.table(all_luad.dat, file="data/interim/all_luad.dat", sep="\t")
