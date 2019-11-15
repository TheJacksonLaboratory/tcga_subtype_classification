load("../visualization/cup_project/subtype_predictor_results.rdata")
output <- NULL
max_length <- 5926
for(i in 1:length(model.list)){
  cancerType <- names(model.list)[i]
  m <- model.list[[i]]
  col <- rownames(m$importance)
  length(col) <- max_length
  names(col) <- cancerType
  output <- cbind(output, col)
}

colnames(output) <- names(model.list)
output <- output[-1,]#first row is X
write.csv(output, file='subtype_classifier_features.csv', na = "",
          row.names = FALSE)

