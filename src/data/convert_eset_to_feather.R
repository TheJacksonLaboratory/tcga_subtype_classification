#!/usr/bin/env Rscript

library(feather)

# loads combined.eset
load('../../data/raw/combined_data.rdata')

head(combined.eset@phenoData)

samples <- phenoData(combined.eset)@data
samples$index <- rownames(samples)
rownames(samples) <- NULL

features <- featureData(combined.eset)@data
features$index <- rownames(features)
rownames(features) <- NULL

expr <- as.data.frame(exprs(combined.eset))
expr$index <- rownames(expr)
rownames(expr) <- NULL

write_feather(expr, '../../data/raw/combined_expr.feather')
write_feather(samples, '../../data/raw/combined_samp.feather')
write_feather(features, '../../data/raw/combined_feat.feather')
