# Atlas clustering dataset
# create random forest classificator

# load libraries
library(caTools)
library(rpart)
library(rpart.plot)
library(randomForest)
library(rattle)
library(caret)
library(patchwork)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(scales)
library(RColorBrewer)
library(pals)
library(ggsci)

atlas <- readRDS("atlas_proportion_dataset_clustered.rds")
atlas$cluster <- as.factor(atlas$cluster_kmeans_k6)
atlas$cluster_kmeans_k6 <- NULL
atlas$tsne_y <- NULL
atlas$tsne_x <- NULL
colnames(atlas) <- make.names(colnames(atlas))

# first we split the dataset into training and test
split <- sample.split(atlas$cluster, SplitRatio = 0.75)
training_set <- subset(atlas[, c(11:ncol(atlas))], split == TRUE)
test_set <- subset(atlas[, c(11:ncol(atlas))], split == FALSE)

# then train the RF model
rf <- randomForest(cluster ~ ., ntree = 1000, data = as.data.frame.data.frame(training_set), importance = TRUE, proximity = TRUE)
plot(rf)

# fit model on test subset
predicted_cluster <- predict(rf, test_set[1:25])
# evaluate with the confusion matrix
confusionMatrix(
  data = predicted_cluster,
  reference = test_set$cluster
)

# Variable Importance
varImpPlot(rf, type = 1, sort = TRUE, n.var = 15, main = "Top 15 - Variable Importance")
varImpPlot(rf, type = 2, sort = TRUE, n.var = 15, main = "Top 15 - Variable Importance")

df <- as.data.frame(rf$importance[, c("MeanDecreaseAccuracy", "MeanDecreaseGini")])
df$cell_type <- rownames(df)
df$cell_type <- gsub("\\.", " ", df$cell_type)
colnames(df) <- c("Mean Decrease Accuracy", "Mean Decrease Gini", "cell_type")
df <- pivot_longer(df, cols = c("Mean Decrease Accuracy", "Mean Decrease Gini"), names_to = "measure")
df <- df[order(-df$value), ]
df <- df[df$cell_type %in% df$cell_type[1:15], ]
df$cell_type <- as.factor(df$cell_type)
df$cell_type <- factor(df$cell_type, levels = unique(df$cell_type[order(df$value)]))

ggplot(df, aes(x = cell_type, y = value)) +
  geom_point(col = "tomato2", size = 5) + # Draw points
  coord_flip() +
  facet_grid(. ~ measure, scales = "free", labeller = label_wrap_gen(width = 10)) +
  theme_minimal() +
  labs(x = "", y = "") +
  theme(text = element_text(size = 25)) +
  theme(panel.background = element_rect(fill = NA, color = "black"))
