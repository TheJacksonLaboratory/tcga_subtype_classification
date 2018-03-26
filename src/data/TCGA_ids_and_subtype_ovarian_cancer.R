diff.cases <- scan("data/raw/cases_expression_cluster1.txt",what="character")
diff.ids <- diff.cases[grep("TCGA",diff.cases)]
imm.cases <- scan("data/raw/cases_expression_cluster2.txt",what="character")
imm.ids <- imm.cases[grep("TCGA",imm.cases)]
diff.dat <- data.frame(sampleid = diff.ids,subtype="Differentiated",cancertype="Ovarian")
imm.dat <- data.frame(sampleid = imm.ids,subtype ="Immunoreactive",cancertype="Ovarian")
all.res <- rbind(diff.dat,imm.dat)
mes.cases <- scan("data/raw/cases_expression_cluster3.txt", what="character")
mes.ids <- mes.cases[grep("TCGA", mes.cases)]
mes.dat <- data.frame(sampleid=mes.ids,subtype="Mesenchymal",cancertype="Ovarian")
pro.cases <- scan("data/raw/cases_expression_cluster4.txt", what="character")
pro.ids <-pro.cases[grep("TCGA",pro.cases)]
pro.dat <-data.frame(sampleid=pro.ids,subtype="Proliferative",cancertype="Ovarian")
all_ov.dat <- rbind(diff.dat,imm.dat,mes.dat,pro.dat)


all_ov.mat <- as.matrix(all_ov.dat)
all_ov.TCGAids <- all_ov.mat[grep("TCGA", all_ov.mat)]
ov.sampleid.mat2 <- substr(all_ov.TCGAids, 1, 15)
ov.sampleid.mat3 <- substr(all_ov.TCGAids, 1, 12)
ov.sampleid.dat <- as.data.frame(ov.sampleid.mat2, ov.sampleid.mat3)

rownames(ov.sampleid.dat) -> ov.sampleid.dat$sampleid
rownames(ov.sampleid.dat)=NULL
colnames(ov.sampleid.dat) <- c("sampleid", "patientid")
all_ov.dat <- merge(all_ov.dat, ov.sampleid.dat, by="sampleid", all.x=FALSE)
write.table(all_ov.dat, file="data/interim/all_ov.dat", sep="\t")
