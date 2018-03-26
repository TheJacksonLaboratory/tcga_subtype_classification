bas_lus.cases <- scan("data/raw/lusc_cases_basal.txt",what="character")
bas_lus.ids <- bas_lus.cases[grep("TCGA", bas_lus.cases)]
bas_lus.dat <- data.frame(sampleid=bas_lus.ids, subtype="Basal",cancertype="Lung squamous")

class_lus.cases <- scan("data/raw/lusc_cases_classical.txt", what="character")
class_lus.ids <- class_lus.cases[grep("TCGA", class_lus.cases)]
class_lus.dat <- data.frame(sampleid=class_lus.ids, subtype="Classical",cancertype="Lung squamous")

prim.cases <- scan("data/raw/cases_primitive.txt",what="character")
prim.ids <- prim.cases[grep("TCGA", prim.cases)]
prim.dat <- data.frame(sampleid=prim.ids, subtype="Primitive",cancertype="Lung squamous")

sec.cases <- scan("data/raw/cases_secretory.txt",what="character")
sec.ids <- sec.cases[grep("TCGA", sec.cases)]
sec.dat <- data.frame(sampleid=sec.ids, subtype="Secretory",cancertype="Lung squamous")

all_lusc.dat <- rbind(bas_lus.dat, class_lus.dat, prim.dat, sec.dat)


all_lusc.mat <- as.matrix(all_lusc.dat)
all_lusc.TCGAids <- all_lusc.mat[grep("TCGA", all_lusc.mat)]
lusc.sampleid.mat2 <- substr(all_lusc.TCGAids, 1, 15)
lusc.sampleid.mat3 <- substr(all_lusc.TCGAids, 1, 12)
lusc.sampleid.dat <- as.data.frame(lusc.sampleid.mat2, lusc.sampleid.mat3)

rownames(lusc.sampleid.dat) -> lusc.sampleid.dat$sampleid
rownames(lusc.sampleid.dat)=NULL
colnames(lusc.sampleid.dat) <- c("sampleid", "patientid")
all_lusc.dat <- merge(all_lusc.dat, lusc.sampleid.dat, by="sampleid", all.x=FALSE)
write.table(all_lusc.dat, file="data/interim/all_lusc.dat", sep="\t")
