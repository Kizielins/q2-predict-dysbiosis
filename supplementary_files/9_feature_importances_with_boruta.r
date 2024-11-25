library(Boruta)
library(randomForest)
library(AUC)
library(gplots)
library(ROCR)
library(pROC)
library(xgboost)


# Load the data

# Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if a cohort argument is provided
if (length(args) < 1) {
    stop("Please provide a cohort argument.")
}

# Define the cohort variable from command line argument
cohort <- args[1]

input_filename <- paste0("boruta_results/", cohort, "_input.txt")
data <- read.table(input_filename, sep='\t', header=TRUE)
#head(data)
# Ensure 'binary' column is treated as a factor (since it's the target)
data$category <- as.factor(data$category)
    
    # Run Boruta on the cohort-specific data without the 'cohort' column
boruta_result <- Boruta(category ~ ., data = data, getImp = getImpXgboost) #getImp = getImpRfRaw
    
    # Display results
print(boruta_result)
    
    # Save importance history for the cohort
output_filename <- paste0("boruta_results/", cohort, ".csv")
write.csv(boruta_result$ImpHistory, output_filename, row.names = FALSE)
    
    # Save feature rankings (mean importance, decision) for the cohort
feature_ranking <- attStats(boruta_result)
ranking_output_filename <- paste0("boruta_results/feature_ranking_", cohort, ".csv")
write.csv(feature_ranking, ranking_output_filename, row.names = TRUE)

