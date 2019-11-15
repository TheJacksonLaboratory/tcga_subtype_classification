#!/usr/bin/env Rscript

write_outputs <- function(model, outdir, prefix) {
  tab <- model$table
  tab.fname <- paste(prefix, "_contingency_table.csv", sep = "")
  tab.path <- paste(outdir, tab.fname, sep = "/")
  write.csv(tab, tab.path)

  pm <- model$byClass
  pm.fname <- paste(prefix, "_performance_metrics.csv", sep = "")
  pm.path <- paste(outdir, pm.fname, sep = "/")
  write.csv(pm, pm.path)
}

print("Extracting metrics and contingency tables from models.")
print("Extracting primary type results.")

rel.dir <- 'models/primary_type'
load(paste(rel.dir, 'final_model.rdata', sep = '/'))
write_outputs(rel.dir, final.model, 'primary')
rm(list=ls())

rel.dir <- 'models/primary_type_ext_val'
load(paste(rel.dir, 'primary_ext_val_models.rdata', sep = '/'))
write_outputs(rel.dir, pdx.cm, 'pdx')
write_outputs(rel.dir, gse18549.cm, 'gse18549')
write_outputs(rel.dir, gse2109.cm, 'gse2109')
rm(list=ls())

print("Extracting subtype results.")
rel.dir <- 'models/subtype'
load(paste(rel.dir, 'subtype_predictor_results.rdata', sep = '/'))
tumors = names(RF.res)
for(i in 1:length(tumors))
{
  write_outputs(rel.dir, RF.res[[i]], tumors[i])
}
rm(list=ls())

rel.dir <- 'models/subtype_ext_val'
load(paste(rel.dir, 'subtype_classifier_results_from_metabric_samples.rdata', sep = '/'))
write_outputs(rel.dir, breast.cm, "breast")

load(paste(rel.dir, 'subtype_classifier_results_from_AOCS_samples.rdata', sep = '/'))
write_outputs(rel.dir, ovarian.cm, "ovarian")
