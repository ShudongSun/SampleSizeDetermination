% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_fitting.R
\name{get_predictions}
\alias{get_predictions}
\title{Run and get predicted results}
\usage{
get_predictions(train_data, test_data, model = "randomforest", func = NULL)
}
\arguments{
\item{train_data}{a dataframe including training data x and the responding variable y as "class"}

\item{test_data}{a dataframe including test data x and the responding variable y as "class"}

\item{model}{base classification model.
\itemize{
\item svm: Support Vector Machines. \code{\link[e1071]{svm}} in \code{e1071} package
\item randomforest: Random Forest. \code{\link[randomForest]{randomForest}} in \code{randomForest} package
\item tree: Classificatin Tree. \code{\link[tree]{tree}} in \code{tree} package
\item self: You can use your self-defined function. You need to pass your self-defined function via the "func" parameter.
}}

\item{func}{If you set "model" to "self", you have to pass your self-defined model function. This function should be able to take "x_train" and "y_train" as the first two inputs to train the model and then take "x_test" as the third input and return the predicted scores of x_test data. For example,
\preformatted{library(e1071)

predict_model <- function(x_train, y_train, x_test){
     data_trainxy<-data.frame(x_train,y_train=as.factor(y_train))
     fit_svm<-svm(y_train~.,data=data_trainxy,probability=TRUE)
     pred_svm <- predict(fit_svm, x_test, probability=TRUE,decision.values = TRUE)
     p_svm=as.data.frame(attr(pred_svm, "probabilities"))$"1"
     return(p_svm)
}

result = get_p_result(data_list=data_list, model=c("self","randomforest"), func=predict_model)}}
}
\value{
the scores predicted by models
}
\description{
Run statistical model and get AUC results.
}
