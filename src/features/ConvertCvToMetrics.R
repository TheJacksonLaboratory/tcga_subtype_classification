library(caret)
#Inception
#fileName <- "Confusion_matrix_CV_InceptionNet1d_diverse_DEG_Feature_40_deep_orderChromo_scan.csv"
#fileName <- "Confusion_matrix_Meta_InceptionNet1d_diverse_DEG_Feature_40_deep_orderChromo_scan.csv"
#cnn 1d
#fileName <- "Confusion_matrix_CV_CNN1d_DEG_Feature_40_smote_OrderChromo.csv"
#fileName <- "Confusion_matrix_Meta_CNN1d_DEG_Feature_40_smote_OrderChromo.csv"
#resnet
#fileName <- "Confusion_matrix_CV_ResNet_DEG_Feature_1000_smote_orderChromo_Best.csv"
#fileName <- "Confusion_matrix_Meta_ResNet_DEG_Feature_1000_smote_orderChromo_Best.csv"

Confusion_matrix <- read.csv(fileName, row.names=1, stringsAsFactors=FALSE)
#row as true, column a prediction
true_values <- c()
predict_values <- c()
for(i in 1:nrow(Confusion_matrix)){
  for(j in 1:ncol(Confusion_matrix)){
    v <- as.numeric(Confusion_matrix[i,j])
    if (v > 0){
      true_values <- c(true_values,rep(rownames(Confusion_matrix)[i], v))
      predict_values <- c(predict_values, rep(colnames(Confusion_matrix)[j], v))
    }
  }
}

levels <- rownames(Confusion_matrix)
predict_values  <- factor(predict_values, levels = levels)
true_values  <- factor(true_values, levels=levels)
cm <- confusionMatrix(predict_values, true_values)
write.csv(cm[[2]],paste0("contingency_table_",fileName))
write.csv(cm[[4]],paste0("preformance_metrics_",fileName))
