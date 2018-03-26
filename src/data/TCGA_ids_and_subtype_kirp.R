c1_kirp.cases <- scan("data/raw/cases_mRNA_cluster_1.txt",what="character")
c1_kirp.ids <- c1_kirp.cases[grep("TCGA", c1_kirp.cases)]
c1_kirp.dat <- data.frame(patientid=c1_kirp.ids, subtype="type I papillary and MET mutations",cancertype="Papillary renal")

c2_kirp.cases <- scan("data/raw/cases_mRNA_cluster_2.txt",what="character")
c2_kirp.ids <- c2_kirp.cases[grep("TCGA", c2_kirp.cases)]
c2_kirp.dat <- data.frame(patientid=c2_kirp.ids, subtype="Mostly type II papillary or unclassified",cancertype="Papillary renal")

c3_kirp.cases <- scan("data/raw/cases_mRNA_cluster_3.txt", what="character")
c3_kirp.ids <- c3_kirp.cases[grep("TCGA", c3_kirp.cases)]
c3_kirp.dat <- data.frame(patientid=c3_kirp.ids, subtype="Best survival",cancertype="Papillary renal")

all_kirp <-rbind(c1_kirp.dat, c2_kirp.dat, c3_kirp.dat)
all_kirp_order <- all_kirp[order(all_kirp$patientid), ]

kirp.sampleid.cases <- scan("data/raw/kirp_cases_all.txt",what="character")
kirp.sampleid.ids <- kirp.sampleid.cases[grep("TCGA", kirp.sampleid.cases)]
kirp.sampleid.ids.mat <-matrix(kirp.sampleid.ids)
kirp.sampleid.ids.mat[order(kirp.sampleid.ids.mat),]
kirp.sampleid.ids.uniq.mat <- unique(kirp.sampleid.ids.mat)


kirp.sampleid.ids.mat2 <- substr(kirp.sampleid.ids.uniq.mat, 1, 12)
kirp.patientids.mat <- substr(kirp.sampleid.ids.uniq.mat, 1, 15)
kirp_patient_sample_ids.dat <- as.data.frame(kirp.patientids.mat, kirp.sampleid.ids.mat2)

rownames(kirp_patient_sample_ids.dat) -> kirp_patient_sample_ids.dat$patientid
rownames(kirp_patient_sample_ids.dat)=NULL
colnames(kirp_patient_sample_ids.dat) <- c("sampleid", "patientid")


all_kirp.dat <- merge(all_kirp_order, kirp_patient_sample_ids.dat, by="patientid", all.x=FALSE)
write.table(all_kirp.dat, file="data/interim/all_kirp.dat", sep="\t")
