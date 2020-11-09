# This script contains wrapper functions to build a SVM model and predict new
# data.

require(e1071)

build_model <- function(data, labels) {
  
  if (sum(is.na(data)) > 0) {
    stop('There are missing values in the data')
  }
  
  model <- svm(x = t(data), y = as.factor(labels))
  return(model)
}

predict_new_data <- function(new_data, model) {
  
  predictions <- predict(model, t(new_data))
  return(predictions)
}

performance <- function(labels, predictions) {
  
  print(table(predictions, labels))
  
}



