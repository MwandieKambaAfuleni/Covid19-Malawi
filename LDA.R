
library(MASS)
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyverse)
library(caret)
theme_set(theme_classic())

#============================================================================
# Linear Discriminant analysis - dates separating the waves
#============================================================================

## Data extracetd from patient metadata for sequences 
dd<-read.csv("./data/lda.csv", header = TRUE)
dd<-dd %>% mutate(date =dmy(CollectionDate)) 
dd <- data.frame(transform(dd, day = as.numeric(date)))
ggplot(dd, aes(x=day, y=Variant)) + geom_point(aes(color = Variant))
dd <- data.frame("day" = dd$day, "Variant"=dd$Variant)
scaledPred <- scale(dd[1])
apply(scaledPred, 2, mean)
apply(scaledPred, 2, sd) 


#make this example reproducible
set.seed(1)
#Use 70% of dataset as training set and remaining 30% as testing set
sample <- sample(c(TRUE, FALSE), nrow(dd), replace=TRUE, prob=c(0.7,0.3))
train <- dd[sample, ]
test <- dd[!sample, ] 
#fit LDA model
model <- lda(Variant~., data=train)
#view model output
model

##### separating dates:
mean(model$means[4,],model$means[1,]) #
mean(model$means[1,],model$means[2,])
mean(model$means[2,],model$means[3,])


#=====================================







#use LDA model to make predictions on test data
predicted <- predict(model, test)
#view predicted class for first six observations in test set
head(predicted$class)
#view posterior probabilities for first six observations in test set
head(predicted$posterior)
#view linear discriminants for first six observations in test set
head(predicted$x)
#find accuracy of model
mean(predicted$class==test$Variant)   
#define data to plot
lda_plot <- cbind(train, predict(model)$x)
#create plot
ggplot(lda_plot, aes(day, LD1)) +   #LD2, Variant, day
  geom_point(aes(color = Variant))




#============================================================================





















