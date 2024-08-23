#!/usr/bin/env Rscript 

options(stringsAsFactors = F) # important !
#Rscript random_forest_tissue_classification_cv.R --feature_matrix path_to_feature_matrix.tsv --annotation path_to_annotation.tsv --trees 500 --mtry 10 --folds 5 --scaling standard --output output_directory/ --model_output saved_model.rds

# Load necessary libraries
if(!require(randomForest)) install.packages("randomForest")
if(!require(caret)) install.packages("caret")
if(!require(argparse)) install.packages("argparse")
if(!require(scales)) install.packages("scales")

library(randomForest)
library(caret)
library(argparse)
library(scales)
library(data.table)

# Define command-line arguments using argparse
parser <- ArgumentParser(description='Random Forest Tissue Classification with K-fold Cross-Validation')
parser$add_argument('-f', '--feature_matrix', type='character', required=TRUE,
                    help='Path to the feature count matrix (TSV file, features as rows, samples as columns)')
parser$add_argument('-a', '--annotation', type='character', required=TRUE,
                    help='Path to the sample-tissue annotation table (TSV file with "sample" and "group" columns)')
parser$add_argument('-s', '--scaling', type='character', choices=c("none", "standard"), default="none",
                    help='Scaling method for feature matrix [choices: none, standard; default=none]')
parser$add_argument('-o', '--output', type='character', default='output/',
                    help='Directory for output files [default=output/]')
parser$add_argument('-mi', '--model_input', type='character', default='random_forest_model.rds',
                    help='Filename for reading the trained model [default=random_forest_model.rds]')
parser$add_argument('-fsel', '--feature_selection', type='character', #default=5000,
                    help='path of selected feature')

# Parse command line arguments
args <- parser$parse_args()
for(i in 1:length(args)){
  v.name <- names(args)[i]
  assign(v.name,args[[i]])
  message(paste0(v.name,": ",get(v.name)))
}


# # test
# feature_matrix="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE163534/call_peak_all/count_matrix/cfpeakCNN_b5_d50_p1_CPM.txt"
# annotation="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/GSE163534/sample_table_filter2.txt"
# output="/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/GSE163534/TOO/"
# scaling="standard"
# model_input="/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/TCGA_small16/TOO_GSE163534/random_forest_model.rds"
# feature_selection="/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/output/TCGA_small16/TOO_GSE163534/feature_importance.csv"

# run
# Load the feature matrix and annotation table (assuming TSV format)
# feature_matrix <- read.csv(feature_matrix, sep="\t", row.names=1, header = T, check.names = F, stringsAsFactors = F)
# annotation <- read.csv(annotation, sep="\t", header = T, check.names = F, stringsAsFactors = F)
feature_matrix <- fread(feature_matrix, sep="\t", header=TRUE, check.names = F, stringsAsFactors = F, data.table = F)
rownames(feature_matrix) <- feature_matrix[,1]
rownames(feature_matrix) <- unlist(sapply(strsplit(rownames(feature_matrix),"|",fixed = T),"[",4))
feature_matrix[,1] <- NULL
annotation <- fread(annotation, sep="\t", header=TRUE, check.names = F, stringsAsFactors = F, data.table = F)
feature_selection <- fread(feature_selection, sep=",", header=TRUE, quote = "\"", check.names = F, stringsAsFactors = F, data.table = F)
table(rownames(feature_matrix) %in% feature_selection$V1)
feature_matrix <- feature_matrix[feature_selection$V1,]
feature_matrix[1:3,1:3]

rownames(annotation) <- annotation$sample
sid <- intersect(colnames(feature_matrix),annotation$sample)
feature_matrix <- feature_matrix[,sid]
annotation <- annotation[sid,]


# # Remove features with zero standard deviation
# feature_sd <- apply(feature_matrix, 1, sd)
# non_zero_sd_features <- feature_sd != 0
# feature_matrix <- feature_matrix[non_zero_sd_features, ]
# dim(feature_matrix)
# feature_matrix[1:3,1:3]
# # tmp1 <- (feature_matrix[,"TCGA-AX-A1CI-11A-11R-A136-13_mirna_gdc_realn"])
# 
# # Check if the feature matrix is empty after removing zero SD features
# if (nrow(feature_matrix) == 0) {
#   stop("Feature matrix is empty after removing features with zero standard deviation. Please check your data.")
# }


feature_matrix_row_names <- rownames(feature_matrix)
feature_matrix_col_names <- colnames(feature_matrix)
# # Apply feature selection based on variance
# feature_variance <- apply(feature_matrix, 1, sd)
# top_features <- order(feature_variance, decreasing=TRUE)[1:feature_selection]
# feature_matrix <- feature_matrix[top_features, ]
dim(feature_matrix)


# Apply scaling method if specified
if (scaling == "standard") {
  message(scaling," scaling")
  feature_matrix <- t(apply(feature_matrix, 1, scale))
}
dim(feature_matrix)
# tmp2 <- (feature_matrix[,"TCGA-AX-A1CI-11A-11R-A136-13_mirna_gdc_realn"])
# cor.test(tmp1,tmp2,method = "spearman")
# sd(tmp1)

# Transpose the feature matrix so that rows are samples and columns are features
feature_matrix <- t(feature_matrix)
#dplyr::mutate(across(, ~as.numeric(.)))
feature_matrix <- apply(feature_matrix,2,as.numeric)
colnames(feature_matrix) <- feature_matrix_row_names
rownames(feature_matrix) <- feature_matrix_col_names

# Ensure that the sample names in the annotation match the feature matrix
#rownames(annotation) <- annotation$sample
annotation <- annotation[rownames(feature_matrix),]
# summary(feature_matrix["TCGA-AX-A1CI-11A-11R-A136-13_mirna_gdc_realn",])
# rownames(feature_matrix)


# Combine feature matrix with tissue labels
data <- as.data.frame(cbind(Tissue = annotation$group, feature_matrix))
# Convert all columns (except the first one) to numeric
data[,2:ncol(data)] <- data.frame(lapply(data[,2:ncol(data)], function(x) if(is.character(x)) as.numeric(as.character(x)) else x))
data[1:3,1:5]
# str(data)
print(unique(data$Tissue))

# Convert the Tissue column to a factor
data$Tissue <- as.factor(data$Tissue)
#table(is.na(data))
#table(is.nan(data))
data[is.na(data)] <- 0

# Set up cross-validation
set.seed(123)  # for reproducibility
# folds2 <- createFolds(data$Tissue, k = folds, list = TRUE)
# 
# # Initialize a variable to store overall confusion matrix results
overall_conf_matrix <- NULL
# 
# # Perform K-fold cross-validation
# for(i in seq_along(folds2)) {
#   # i <- 1
#   message(i)
#   trainIndex <- folds2[[i]]
#   
#   # Split into training and test data
#   trainData <- data[trainIndex, ]
#   testData <- data[-trainIndex, ]
  
  # # Train Random Forest model
  # rf_model <- randomForest(Tissue ~ ., data = trainData, importance = TRUE, 
  #                          ntree = trees, mtry = mtry)
  # # trainData[,"ENST00000385300_____1_60_81_+|pri_miRNA|peak_6065|peak_6065|ENST00000385300_____1|60|81|pri_miRNA"] <- as.numeric(as.character( trainData[,"ENST00000385300_____1_60_81_+|pri_miRNA|peak_6065|peak_6065|ENST00000385300_____1|60|81|pri_miRNA"]))
testData <- data
rf_model <- readRDS(model_input)
  # Predict on test data
  rf_predictions <- predict(rf_model, newdata = testData)
  
  # Evaluate the model
  conf_matrix <- confusionMatrix(rf_predictions, testData$Tissue)
  # print(conf_matrix)
  
  # Sum the confusion matrix for overall performance
  if (is.null(overall_conf_matrix)) {
    overall_conf_matrix <- conf_matrix$table
  } else {
    overall_conf_matrix <- overall_conf_matrix + conf_matrix$table
  }
# }

# Calculate overall accuracy
overall_accuracy <- sum(diag(overall_conf_matrix)) / sum(overall_conf_matrix)
print(paste("Overall Accuracy: ", overall_accuracy))

# Save the overall confusion matrix and model from the last fold
if (!dir.exists(output)) {
  dir.create(output)
}

write.csv(as.data.frame(overall_conf_matrix), file=paste0(output, "overall_confusion_matrix.csv"))
#saveRDS(rf_model, file=paste0(output, model_output))

# Optional: Feature importance from the last model
importance_df <- as.data.frame(importance(rf_model))
write.csv(importance_df, file=paste0(output, "feature_importance.csv"))

# Plot and save feature importance from the last model
png(filename = paste0(output, "feature_importance_plot.png"))
varImpPlot(rf_model)
dev.off()





# Calculate average confusion matrix
average_conf_matrix <- overall_conf_matrix 
print("Average Confusion Matrix:")
print(average_conf_matrix)
# Save the average confusion matrix
if (!dir.exists(output)) {
  dir.create(output)
}
write.csv(as.data.frame(average_conf_matrix), file=paste0(output, "average_confusion_matrix.csv"))

# Prepare data for ggplot2
conf_matrix_df <- as.data.frame(as.table(average_conf_matrix))
colnames(conf_matrix_df) <- c("Prediction", "Reference", "Frequency")

# Plot confusion matrix with ggplot2
ggplot(conf_matrix_df, aes(x = Prediction, y = Reference, fill = Frequency)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  labs(title = "Confusion Matrix", x = "Predicted Label", y = "True Label") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
ggsave(filename = paste0(output, "average_confusion_matrix_plot.png"), width = 5, height = 4)
