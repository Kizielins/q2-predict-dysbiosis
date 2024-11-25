X<- read.delim(file = "Data/calculated_parameters_for_merged.txt", sep = "\t")
metadata<- read.delim(file= "Data/merged_metadata.txt", sep="\t")
library(dplyr)

metadata$cohort %>% unique() %>% sort() -> coh_names
lapply(coh_names, function(coh) c(table(metadata[metadata$cohort==coh,]$category), total=sum(metadata$cohort==coh)) )-> class_distr
names(class_distr)<- coh_names
print("if name of the category of samples appears below dataset name it might be problematic to include it in the test")
print("if total appears, dataset is too small overall")
lapply(class_distr, function(x) names(x)[ x<10 ]) #

lapply(class_distr, function(x) min(x)>30) %>% unlist () -> ok_data
class_distr[ok_data] %>% names() -> ok_cohorts
test_meta<- metadata[ metadata$cohort %in% ok_cohorts, ]
test_X<- X[ X$Sample %in% test_meta$sample , ]
rownames(test_X)<- test_X$Sample
rownames(test_meta)<- test_meta$sample
test_X[ rownames(test_meta), ]-> test_X
all(rownames(test_meta)==rownames(test_X))
test_meta$cohort %>% unique() %>% sort()-> test_cohorts

test_subsets<- lapply(test_cohorts, function(cname) test_meta$cohort==cname  )
names(test_subsets)<- test_cohorts
test_Y<- test_meta$category
names(test_Y)<- rownames(test_meta)
test_X<- test_X[,-(1:2)]

binary_Y<- test_Y
for (i in seq_along(test_subsets))
    {
      cohort<- test_subsets[[i]]
      print(sprintf("old class distribution for %s", names(test_subsets)[[i]] ))
      print(table(test_Y[ cohort ]))
      old<- test_Y[ cohort ]
      classNames<- unique(old)
      if (length(classNames)>2)
          {
           new<- old
           new[ old!= "healthy" ] = paste0( classNames[classNames!= "healthy" ], collapse = "_or_" )
           binary_Y[ cohort ]<- new 
          }
       print(sprintf("new class distribution for %s", names(test_subsets)[[i]] ))
        print(table(binary_Y[cohort]))
      }

n_estimators = 500
n_splits = 3   #not used yet
n_repeats = 30 #not used yet
library(randomForest)
set.seed(42)
subset_rfs<- lapply(test_subsets, function(cohort_subset)
        randomForest(x = test_X[cohort_subset,  ], 
                     y = as.factor(binary_Y[cohort_subset]),
                     ntree = n_estimators,importance = TRUE)
                    )
names(subset_rfs)<- test_cohorts

library(pROC)
suppressMessages(
lapply(seq_along(subset_rfs), function(i) auc(response = binary_Y[ test_subsets[[i]] ],  
                                              predictor=subset_rfs[[i]]$votes[,1] ) ) -> aucs
       )
names(aucs)<- test_cohorts
library(Boruta)

lapply(test_subsets, function(cohort_subset)
    Boruta(x = test_X[cohort_subset,] , y= as.factor(binary_Y[ cohort_subset ]), getImp=getImpLegacyRfRaw )
    )-> list_of_borutas

mean_imps_boruta<- lapply(list_of_borutas, function(rokita)
    rokita$ImpHistory %>% as.matrix() %>% matrixStats::colMeans2()
    )
imps_RF<- lapply(subset_rfs, function(rf)
                rf$importance[,"MeanDecreaseAccuracy"] )

cohort_imp_dfs<-list()
for (i in seq_along(mean_imps_boruta))
    {    rbind(mean_imps_boruta[[i]][1:6],
          imps_RF[[i]] ) %>% as.data.frame() -> df_cohort
          df_cohort$model= c("boruta","RF")
          df_cohort$AUC_health_pred= rep(aucs[[i]],2)
          cohort_imp_dfs[[i]]<-df_cohort
     }
names(cohort_imp_dfs)<- names(imps_RF)

do.call(rbind,
        lapply(cohort_imp_dfs, function(chunk)
    {
        chunk[1,1:6] %>% as.numeric() %>% is.finite() -> confirmed
        chunk[2,c(1:6,8)]
    }
    ))-> num_heatmap
rownames(num_heatmap)<- paste0(rownames(num_heatmap),"_AUC=", round(num_heatmap$AUC_health_pred,2))
num_heatmap<- num_heatmap[, - ncol(num_heatmap)]
#rank_heatmap<- t(apply(num_heatmap,1,rank))
#
rank_heatmap <- t(apply(num_heatmap, 1, function(x) rank(-x)))

library(ggplot2)
library(reshape2)
#melt(rank_heatmap)-> heatmap_flat
#colnames(heatmap_flat)=c("x","y","fill")

# Save the plot as PNG
#png("feature_ranks_rf_importance.png")

#options(repr.plot.height=13)
#ggplot(heatmap_flat, aes(x=y,y=x, fill=fill)) + geom_tile() +
# scale_fill_gradient(low = "#8bde4b", high = "#bf3043") +
#geom_text(aes(label = fill, size = 2)) +
#theme(text=element_text(size=16)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#ggtitle("Feature ranks per cohort")

#dev.off()

# Calculate distance matrix and perform hierarchical clustering

# Calculate distance matrix and perform hierarchical clustering on rows
dist_matrix <- dist(rank_heatmap, method = "euclidean")  # Use Euclidean distance for rows
row_clustering <- hclust(dist_matrix, method = "ward.D2")  # Hierarchical clustering

# Reorder the rows of rank_heatmap based on clustering results
rank_heatmap <- rank_heatmap[row_clustering$order, ]  # Reorder rows by the clustering order
write.table(rank_heatmap, file = "/mnt/c/Users/kizie/OneDrive - Uniwersytet Jagielloński/Health_index/Clean_analysis/clustered_rank_heatmap.txt", sep = "\t", col.names = NA, quote = FALSE)

# Reshape the matrix for ggplot
melt(rank_heatmap) -> heatmap_flat
colnames(heatmap_flat) = c("x", "y", "fill")

# Save the clustered heatmap as PNG
#png("clustered_feature_ranks_rf_importance.png", width = 1000, height = 1300, res = 150)

options(repr.plot.height = 13)
ggplot(heatmap_flat, aes(x = y, y = x, fill = fill)) +
  geom_tile() +
  scale_fill_gradient(high = "#bf3043", low = "#8bde4b") + 
  geom_text(aes(label = fill), size = 2) +
  theme(text = element_text(size = 16)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ggtitle("Clustered Feature Ranks (Higher Importance = Lower Rank) per Cohort")

#dev.off()

