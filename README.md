# Random Forest
Course materials for Oncode masterclass on machine learning 202
SVM vs. Random forest comparision
Pre-processing
Split into labels classes/features
* How many features 
* What are the classes/labels?
* Type of data (Categorical, numerical)
* Is there any missing data?

Imputation vs. missing data

Scaling/Normalisation
* RPKM, TMM

Training/test split | k-folds | leave-one-out

Parameter tuning, pruning
What is a good feature? Splitting criterion selection 

RandomForest(Train_X, Train_y)
Predict(X_test)
ROC-curve (accuracy vs. sensitivity)
Feature-importance 

Pick top-features, compare to original paper

Pitfalls - drawbacks of random forest

Advanced excercise
Out-of baggings error
* Select number of trees which minimises OOB error and re-run
* 
