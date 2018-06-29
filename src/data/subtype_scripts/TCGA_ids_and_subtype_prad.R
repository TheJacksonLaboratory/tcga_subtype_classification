ERG.cases <- scan("data/raw/cases_ERG.txt",what="character")
ERG.ids <- ERG.cases[grep("TCGA", ERG.cases)]
ERG.dat <- data.frame(patientid=ERG.ids, subtype="ERG",cancertype="Prostate")

ETV1.cases <- scan("data/raw/cases_ETV1.txt",what="character")
ETV1.ids <- ETV1.cases[grep("TCGA", ETV1.cases)]
ETV1.dat <- data.frame(patientid=ETV1.ids, subtype="ETV1",cancertype="Prostate")

ETV4.cases <- scan("data/raw/cases_ETV4.txt",what="character")
ETV4.ids <- ETV4.cases[grep("TCGA", ETV4.cases)]
ETV4.dat <- data.frame(patientid=ETV4.ids, subtype="ETV4",cancertype="Prostate")

FLI1.cases <-scan("data/raw/cases_FLI1.txt",what="character")
FLI1.ids <- FLI1.cases[grep("TCGA", FLI1.cases)]
FLI1.dat <- data.frame(patientid=FLI1.ids, subtype="FLI1",cancertype="Prostate")

FOXA1.cases <- scan("data/raw/cases_FOXA1.txt",what="character")
FOXA1.ids <- FOXA1.cases[grep("TCGA", FOXA1.cases)]
FOXA1.dat <- data.frame(patientid=FOXA1.ids, subtype="FOXA1",cancertype="Prostate")

IDH1.cases <- scan("data/raw/cases_IDH1.txt",what="character")
IDH1.ids <- IDH1.cases[grep("TCGA", IDH1.cases)]
IDH1.dat <- data.frame(patientid=IDH1.ids, subtype="IDH1",cancertype="Prostate")

SPOP.cases <-scan("data/raw/cases_SPOP.txt",what="character")
SPOP.ids <- SPOP.cases[grep("TCGA", SPOP.cases)]
SPOP.dat <- data.frame(patientid=SPOP.ids, subtype="SPOP",cancertype="Prostate")

all_prad <- rbind(ERG.dat, ETV1.dat, ETV4.dat, FLI1.dat, FOXA1.dat, IDH1.dat, SPOP.dat)

prad.sampleid.cases <- scan("data/raw/prad_data_clinical_sample.txt",what="character")
prad.sampleid.ids <- prad.sampleid.cases[grep("TCGA", prad.sampleid.cases)]
prad.sampleid.dat <- data.frame(sampleid=prad.sampleid.ids)

prad.sample_id_patient_id <- as.matrix(prad.sampleid.dat)
n <- length(prad.sample_id_patient_id)
prad.patient_id <-prad.sample_id_patient_id[seq(2,n,2)]
prad.sample_id <- prad.sample_id_patient_id[seq(1,n,2)]
prad.patient_id_sample_id.dat <- data.frame(patientid = prad.patient_id, sampleid=prad.sample_id)
prad.patient_id_sample_id.dat[order(prad.patient_id_sample_id.dat$patientid),]

all_prad.dat <- merge(all_prad, prad.patient_id_sample_id.dat, by="patientid", all.x=FALSE)
write.table(all_prad.dat, file="data/interim/all_prad.dat", sep="\t")
