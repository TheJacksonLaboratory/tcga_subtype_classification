luma.cases <-scan("data/raw/cases_luma.txt",what="character")
luma.ids <-luma.cases[grep("TCGA",luma.cases)]
luma.dat <- data.frame(sampleid = luma.ids, subtype="Luminal A",cancertype="Breast")

her2e.cases <-scan("data/raw/cases_her2.txt",what="character")
her2e.ids <- her2e.cases[grep("TCGA",her2e.cases)]
her2e.dat <- data.frame(sampleid = her2e.ids, subtype="Her2 enriched",cancertype="Breast")

bas.cases <-scan("data/raw/brca_cases_basal.txt",what="character")
bas.ids <- bas.cases[grep("TCGA",bas.cases)]
bas.dat <- data.frame(sampleid = bas.ids, subtype="Basal-like",cancertype="Breast")

lumb.cases <-scan("data/raw/cases_lumb.txt",what="character")
lumb.ids <- lumb.cases[grep("TCGA", lumb.cases)]
lumb.dat <- data.frame(sampleid = lumb.ids, subtype="Luminal B",cancertype="Breast")

all_brca.dat <- rbind(luma.dat,lumb.dat,her2e.dat,bas.dat)


all_brca.mat <- as.matrix(all_brca.dat)
all_brca.TCGAids <- all_brca.mat[grep("TCGA", all_brca.mat)]
brca.sampleid.mat2 <- substr(all_brca.TCGAids, 1, 15)
brca.sampleid.mat3 <- substr(all_brca.TCGAids, 1, 12)
brca.sampleid.dat <- as.data.frame(brca.sampleid.mat2, brca.sampleid.mat3)

rownames(brca.sampleid.dat) -> brca.sampleid.dat$sampleid
rownames(brca.sampleid.dat)=NULL
colnames(brca.sampleid.dat) <- c("sampleid", "patientid")
all_brca.dat <- merge(all_brca.dat, brca.sampleid.dat, by="sampleid", all.x=FALSE)
write.table(all_brca.dat, file="data/interim/all_brca.dat", sep="\t")
