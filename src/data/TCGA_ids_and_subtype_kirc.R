m1.cases <- scan("data/raw/cases_m1.txt",what="character")
m1.ids <- m1.cases[grep("TCGA", m1.cases)]
m1.dat <- data.frame(patientid=m1.ids, subtype="Chromatin remodeling processes",cancertype="CLear cell renal")

m2.cases <- scan("data/raw/cases_m2.txt",what="character")
m2.ids <- m2.cases[grep("TCGA", m2.cases)]
m2.dat <- data.frame(patientid=m2.ids, subtype="High expression of miR-21",cancertype="CLear cell renal")

m3.cases <- scan("data/raw/cases_m3.txt",what="character")
m3.ids <- m3.cases[grep("TCGA", m3.cases)]
m3.dat <- data.frame(patientid=m3.ids, subtype="PTEN mutations and CDKN2A deletions",cancertype="CLear cell renal")

m4.cases <- scan("data/raw/cases_m4.txt",what="character")
m4.ids <- m4.cases[grep("TCGA", m4.cases)]
m4.dat <- data.frame(patientid=m4.ids, subtype="BAP1 and mTOR mutations",cancertype="CLear cell renal")

all_kirc <- rbind(m1.dat, m2.dat, m3.dat, m4.dat)
all_kirc[order(all_kirc$patientid), ]

sampleid.cases <- scan("data/raw/kirc_data_clinical_sample.txt",what="character")
sampleid.ids <- sampleid.cases[grep("TCGA", sampleid.cases)]
sampleid.dat <- data.frame(sampleid=sampleid.ids)

sample_id_patient_id <- as.matrix(sampleid.dat)
n <- length(sample_id_patient_id)
sample_id <-sample_id_patient_id[seq(2,n,2)]
patient_id <- sample_id_patient_id[seq(1,n,2)]
patient_id_sample_id.dat <- data.frame(patientid = patient_id, sampleid=sample_id)
patient_id_sample_id.dat[order(patient_id_sample_id.dat$patientid),]

all_kirc.dat <- merge(all_kirc, patient_id_sample_id.dat, by="patientid", all.x=FALSE)
write.table(all_kirc.dat, file="data/interim/all_kirc.dat", sep="\t")
