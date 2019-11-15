#!/usr/bin/env Rscript

print("Generating RF models and results.")
print("Training primary type RF.")

# saves final_model.rdata
# saves confusion_matrices.rdata (pdx, gse18549, gse2109)
print("Training primary type RF.")
source("src/models/R/RF/predict_tumor_type_RF.R")
print("Validating primary type RF.")
source("src/models/R/RF/predict_tumor_type_validation.R")

# saves subtype_predictor_results.rdata
# saves subtype_classifier_results_from_metabric_samples.rdata
# saves subtype_classifier_results_from_AOCS_samples.rdata
print("Training subtype RF.")
source("src/models/R/RF/predict_tumor_subtype_RF.R")
print("Validating subtype RF.")
source("src/models/R/RF/predict_breast_tumor_subtype_RF.R")
source("src/models/R/RF/predict_ovarian_tumor_subtype_RF.R")


