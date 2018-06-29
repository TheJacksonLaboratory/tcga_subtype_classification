#!/usr/bin/env Rscript

source("src/data/subtype_scripts/TCGA_ids_and_subtype_GBM.R")
source("src/data/subtype_scripts/TCGA_ids_and_subtype_glioma.R")
source("src/data/subtype_scripts/TCGA_ids_and_subtype_hnsc.R")
source("src/data/subtype_scripts/TCGA_ids_and_subtype_kirc.R")
source("src/data/subtype_scripts/TCGA_ids_and_subtype_kirp.R")
source("src/data/subtype_scripts/TCGA_ids_and_subtype_luad.R")
source("src/data/subtype_scripts/TCGA_ids_and_subtype_lusc.R")
source("src/data/subtype_scripts/TCGA_ids_and_subtype_ovarian_cancer.R")
source("src/data/subtype_scripts/TCGA_ids_and_subtype_prad.R")
source("src/data/subtype_scripts/TCGA_ids_and_subtype_skcm.R")
source("src/data/subtype_scripts/TCGA_ids_and_subtype_stad.R")
source("src/data/subtype_scripts/TCGA_ids_and_subtype_ucec.R")
source("src/data/subtype_scripts/TCGA_ids_and_subtypes_brca.R")

prelim_all_subtypes.dat <- rbind(all_glioma.dat, all_GBM.dat, all_hnsc.dat,
                                 all_kirc.dat, all_kirp.dat, all_luad.dat,
                                 all_lusc.dat, all_ov.dat, all_prad.dat,
                                 all_skcm.dat, all_stad.dat, all_brca.dat,
                                 all_ucec.dat)

prelim_all_subtypes.dat <- prelim_all_subtypes.dat[c("sampleid", "patientid",
                                                     "subtype", "cancertype")]


write.table(prelim_all_subtypes.dat,
            file="data/processed/preliminary_TCGA_ids_and_subtypes_for_cancer_subtypes.tab",
            sep="\t")

write.table(prelim_all_subtypes.dat,
            file="data/processed/preliminary_TCGA_ids_and_subtypes_for_cancer_subtypes.csv",
            sep=",")
