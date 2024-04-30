#' @title Sample Size Determination
#' @description build the model and determine the sample size based on elbow method of ARI and AMI.
#'
#' @param x A data frame or a matrix of predictors to be fitted. When using "pilot" mode, you need to use this parameter to input the predictor variable of pilot data; When using "true" mode, you need to use this parameter to input the predictor variable of training data.
#' @param y A response vector. When using "pilot" mode, you need to use this parameter to input the predictor variable of pilot data. When using "true" mode, you need to use this parameter to input the predictor variable of training data.
#' @param model base classification model.
#' \itemize{
#' \item svm: Support Vector Machines. \code{\link[e1071]{svm}} in \code{e1071} package
#' \item randomforest: Random Forest. \code{\link[randomForest]{randomForest}} in \code{randomForest} package
#' \item tree: Classificatin Tree. \code{\link[tree]{tree}} in \code{tree} package
#' \item self: You can use your self-defined function. You need to pass your self-defined function via the "func" parameter.
#' }
#'
#' @param func If you set "model" to "self", you have to pass your self-defined model function. This function should be able to take *train_data_x* and *train_data_y* as the first two inputs to train the model and then take *test_data_x* as the third input and return the predicted scores of *test_data_x* data. For example,
#'
#' \preformatted{library(e1071)
#'
#' predict_model <- function(train_data_x, train_data_y, test_data_x){
#'     train_data = data.frame(train_data_x, as.factor(train_data_y))
#'     names(train_data)[length(train_data)] = "class"
#'     fit_svm<-svm(class~.,data=train_data,probability=TRUE)
#'     pred <- predict(fit_svm, test_data_x)
#'     return(pred)
#' }
#'
#' elbow = ssd(x_pilot, y_pilot, model="self", func=predict_model)
#' }
#'
#' @param index the index used to judge the quality of the model.
#' \itemize{
#' \item classification error: what we usually mean.
#' \item ARI: Adjusted Rand Index.
#' \item AMI: Adjusted Mutual Information.
#' \item self: You can use your self-defined metric function. You need to pass your self-defined function via the "metric_func" parameter.
#' }
#'
#' @param metric_func If you set "index" to "self", you have to pass your self-defined metric function. This function should be able to take *prediction* and *true_label* as the first two inputs to calculate the metric. For example,
#'
#' \preformatted{
#'
#' classification_accuracy <- function(prediction, true_label){
#'     correct <- sum(true_label == prediction)
#'     acc <- correct/length(prediction)
#'     return(acc)
#' }
#'
#' elbow = ssd(x_pilot, y_pilot, index="self", metric_func=classification_accuracy)
#' }
#'
#' @param metric_name If you set "index" to "self", you can pass your metric name here to show in the plot.
#'
#' @param n_train_list a series of numbers which represent the generated training data size for each class.
#' @param n_test the generated test data size for each class.
#' @param mode The default value is "pilot", which means that we will see the input data(x and y) as the pilot data(usually small number), and use them to generate simulated data (corresponding to 'n_train_list'), and then get the curve and elbow point; "true" means that we will sample the data corresponding to 'n_train_list' from input data (x and y), and then get the curve and elbow point.
#' @param num_repeat The number of times to repeat for each training data size in order to get a average value of index to reduce the effect of randomness. The default value is 30.
#' @param print_progress_bar Whether to print the progress bar and elapsed time. The default value is TRUE.
#' @param n.cores The number of nodes to be forked when using multi-core parallel computing. If not being set(n.cores=NULL), \code{n.cores <- parallel::detectCores() - 1} would be used.
#' @param test_x A data frame or a matrix of predictors. When using "true" mode, you need to use this parameter to input the predictor variable of test data.
#' @param test_y A response vector. When using "true" mode, you need to use this parameter to input the predictor variable of test data.
#'
#' @return Return the average classification errors/ARI/AMI of the repeat results for the selected mode.
#' @import mclust MASS aricode progress parallel foreach doParallel ggplot2 minpack.lm
#' @importFrom randomForest randomForest
#' @importFrom stats cov nls pnorm predict qnorm quantile
#'
#' @examples
#' pilot_data <- read.csv(system.file("extdata", "data_pbmc68k_pilot_18pc.csv", package = "SSD"),row.names=1)
#'
#' x_pilot = pilot_data[,-length(pilot_data)]
#' y_pilot = pilot_data[,length(pilot_data)]
#' table(pilot_data$phenoid)
#'
#' result_pilot = ssd(x_pilot, y_pilot)
#'
#' ### use true data:
#'
#' train_data <- read.csv(system.file("extdata", "data_pbmc68k_train_23pc.csv", package = "SSD"),row.names=1)
#' test_data <- read.csv(system.file("extdata", "data_pbmc68k_test_23pc.csv", package = "SSD"),row.names=1)
#'
#' x_true_train = train_data[,-length(train_data)]
#' y_true_train = train_data[,length(train_data)]
#' table(train_data$phenoid)
#' x_true_test = test_data[,-length(test_data)]
#' y_true_test = test_data[,length(test_data)]
#' table(test_data$phenoid)
#'
#' result_true = ssd(x=x_true_train, y=y_true_train, mode="true", test_x=x_true_test, test_y=y_true_test)
#'
#' @export
ssd <- function(x, y, model="randomforest", func=NULL, index="ARI", n_train_list = seq(from=30,to=600,by=30), n_test = 300, mode="pilot", num_repeat=30, print_progress_bar=TRUE, n.cores=NULL, test_x=NULL, test_y=NULL, metric_func=NULL, metric_name=NULL) {

  if(!model %in% c( "svm", "randomforest", "tree","self")){
    stop('\nmodel \'',model, '\' cannot be found')
  }

  if(!index %in% c("ARI","AMI","classification error","self")){
    stop('\nindex \'',index, '\' cannot be found')
  }

  n_iter <- length(n_train_list)+1
  if(print_progress_bar==TRUE){
    pb <- progress_bar$new(format = ":percent [:bar] [Elapsed time: :elapsedfull] :message",
                           total = n_iter,
                           complete = "=",   # Completion bar character
                           incomplete = "-", # Incomplete bar character
                           current = ">",    # Current bar character
                           clear = FALSE,    # If TRUE, clears the bar when finish
                           width = 100,      # Width of the progress bar
                           show_after = 0)
    pb$tick(0,tokens = list(message="[Initializing...]                 "))
  }

  if(is.null(n.cores)){
    n.cores <- parallel::detectCores() - 1
  }
  my.cluster <- parallel::makeCluster(
    n.cores,
    type = "PSOCK"
  )
  doParallel::registerDoParallel(cl = my.cluster)


  if(mode=="pilot"){
    pilot_data = cbind(x,y)
    names(pilot_data)[length(names(pilot_data))] = "class"
    num_class = length(table(y))
    num_PC = length(colnames(x))

    pilot_data_rp = pilot_data
    means = matrix(0, nrow = num_class, ncol = num_PC)
    covs = vector(mode = "list", length = num_class)
    for(i in 1:num_class){
      class_i_ids = which(pilot_data$class == names(table(pilot_data$class))[i])
      pilot_i_data = pilot_data[class_i_ids,]
      for(j in 1:num_PC){
        pilot_data_i_j = pilot_i_data[,j]
        temp = qnorm(rank(pilot_data_i_j)/(length(pilot_data_i_j) + 1)) # rank
        pilot_data_i_j_rp = as.numeric(quantile(pilot_data_i_j, probs = pnorm(temp))) # percentile
        pilot_data_rp[class_i_ids,j] = pilot_data_i_j_rp
        # pilot_data_rp[class_i_ids,j] = pilot_data_i_j
        means[i,j] = mean(pilot_data_i_j_rp)
      }
      covs[[i]] = cov(pilot_data[class_i_ids,1:num_PC])
    }
    p_nclass = as.numeric(table(y))
    for(i in 1:num_class){
      if(i==1){
        sigma_hat = (p_nclass[i]-1)*covs[[i]]
      }else{
        sigma_hat = sigma_hat + (p_nclass[i]-1)*covs[[i]]
      }
    }
    sigma_hat = sigma_hat/(sum(p_nclass)-num_class)

    n_train_list = n_train_list
    metrics_syn <- vector(length=length(n_train_list))
    j=1

    for(n_train in n_train_list){

      n_test = n_test

      if(print_progress_bar==TRUE){pb$tick(tokens = list(message = paste("[Training size =",n_train,"is running...]")))}

      result_temp <- foreach(
        k = 1:num_repeat,
        .combine = 'rbind',
        .packages = c('MASS','randomForest','mclust','aricode','progress')
      ) %dopar% {

        for(i in 1:num_class){
          if(i==1){
            train_data = cbind(data.frame(mvrnorm(n_train, means[i,], sigma_hat)),class = factor(i,levels=1:num_class))
            test_data = cbind(data.frame(mvrnorm(n_test, means[i,], sigma_hat)),class = factor(i,levels=1:num_class))
          }else{
            train_data = rbind(train_data,cbind(mvrnorm(n_train, means[i,], sigma_hat),class = factor(i,levels=1:num_class)))
            test_data = rbind(test_data,cbind(mvrnorm(n_test, means[i,], sigma_hat),class = factor(i,levels=1:num_class)))
          }
        }


        pred = get_predictions(train_data=train_data, test_data=test_data, model=model, func=func)


        # Record error rate for one fold (0 or 1 for LOOCV)
        if(index=="classification error"){
          metrics_temp = misclass_err(pred, test_data[,ncol(test_data)])
        }
        if(index=="ARI"){
          metrics_temp = adjustedRandIndex(test_data[,ncol(test_data)], pred)
        }
        if(index=="AMI"){
          metrics_temp = AMI(test_data[,ncol(test_data)], pred)
        }
        if(index=="self"){
          metrics_temp = metric_func(pred, test_data[,ncol(test_data)])
        }
        # ARIs_temp = adjustedRandIndex(test_data[,ncol(test_data)], pred)
        # AMIs_temp = AMI(test_data[,ncol(test_data)], pred)

        return(metrics_temp)

      }

      metrics_syn[j] = mean(result_temp)

      j=j+1
    }
    result = metrics_syn

  }else if(mode=="true"){
    if(is.null(test_x) | is.null(test_y)){
      stop('\nYou must input test_x and test_y parameters when using this mode.')
    }

    train_data_input = cbind(x,y)
    names(train_data_input)[length(names(train_data_input))] = "class"
    num_class_train = length(table(y))
    num_PC_train = length(colnames(x))

    test_data_input = cbind(test_x, test_y)
    names(test_data_input)[length(names(test_data_input))] = "class"
    num_class_test = length(table(test_y))
    num_PC_test = length(colnames(test_x))

    num_class = num_class_train
    if(num_PC_train != num_PC_test){
      stop('\n# of PCs in training data must be matched with # of PCs in test data.')
    }

    n_train_list = n_train_list
    metrics_true <- vector(length=length(n_train_list))
    j=1
    for(n_train in n_train_list){
      # print(j)
      n_test = n_test

      if(print_progress_bar==TRUE){pb$tick(tokens = list(message = paste("[Training size =",n_train,"is running...]")))}

      result_temp <- foreach(
        k = 1:num_repeat,
        .combine = 'rbind',
        .packages = c('MASS','randomForest','mclust','aricode','progress')
      ) %dopar% {

        # do here
        for(i in 1:num_class){
          class_i_ids_train = which(train_data_input$class == names(table(train_data_input$class))[i])
          if(length(class_i_ids_train)>=n_train){
            train_i_ids = sample(class_i_ids_train, n_train)
            train_i_data = train_data_input[train_i_ids,]
          }else{
            train_i_ids = sample(class_i_ids_train,n_train,replace=TRUE)
            train_i_data = train_data_input[train_i_ids,]
          }

          class_i_ids_test = which(test_data_input$class == names(table(test_data_input$class))[i])
          if(length(class_i_ids_test)>=n_test){
            test_i_ids = sample(class_i_ids_test, n_test)
            test_i_data = test_data_input[test_i_ids,]
          }else{
            test_i_ids = sample(class_i_ids_test,n_test,replace=TRUE)
            test_i_data = test_data_input[test_i_ids,]
          }

          if(i == 1){
            train_data = train_i_data
            test_data = test_i_data
          }else{
            train_data = rbind(train_data, train_i_data)
            test_data = rbind(test_data, test_i_data)
          }
        }
        # table(train_data$class)
        # table(test_data$class)

        pred = get_predictions(train_data=train_data, test_data=test_data, model=model, func=func)

        # Record error rate for one fold (0 or 1 for LOOCV)
        if(index=="classification error"){
          metrics_temp = misclass_err(pred, test_data[,ncol(test_data)])
        }
        if(index=="ARI"){
          metrics_temp = adjustedRandIndex(test_data[,ncol(test_data)], pred)
        }
        if(index=="AMI"){
          metrics_temp = AMI(test_data[,ncol(test_data)], pred)
        }
        if(index=="self"){
          metrics_temp = metric_func(pred, test_data[,ncol(test_data)])
        }
        # ARIs_temp = adjustedRandIndex(test_data[,ncol(test_data)], pred)
        # AMIs_temp = AMI(test_data[,ncol(test_data)], pred)

        return(metrics_temp)

      }

      metrics_true[j] = mean(result_temp)

      j=j+1
    }
    result = metrics_true
  }

  parallel::stopCluster(cl = my.cluster)


  if(mode=="pilot"){results_calculated = metrics_syn}
  if(mode=="true"){results_calculated = metrics_true}

  metric_name_pass = index

  if(index=="self"){
    if(is.null(metric_name)){
      metric_name_pass = "Self-difined Metric"
    }else{
      metric_name_pass = metric_name
    }
  }

  plot_fit_ipl(results_calculated, n_train_list, mode, metric_name_pass, model)


  # if(index %in% c("ARI","AMI")){
  #   result = findPC(sdev = -pred_index_syn,number = c(length(nn)),method = 'all',figure = T)
  # }else if(index == "classification error"){
  #   result = findPC(sdev = pred_index_syn,number = c(length(nn)),method = 'all',figure = T)
  # }
  if(print_progress_bar==TRUE){pb$tick(tokens = list(message="[Finished]                         "))}

  return(result)
}



#' @title calculate classification error
#'
#' @param pred predicted values by the model.
#' @param actual actual values.
#'
#' @return Return the classification error
#' @export
misclass_err <- function(pred, actual){
  incorrect <- sum(!actual == pred)
  correct <- sum(actual == pred)
  err <- 1 - correct/length(pred)

  return(err)
}


#' @importFrom randomForest randomForest
#' @importFrom e1071 svm
#' @importFrom tree tree
#' @title Run and get predicted results
#' @description  Run statistical model and get AUC results.
#' @param train_data a dataframe including training data x and the responding variable y as "class"
#' @param test_data a dataframe including test data x and the responding variable y as "class"
#' @param model base classification model.
#' \itemize{
#' \item svm: Support Vector Machines. \code{\link[e1071]{svm}} in \code{e1071} package
#' \item randomforest: Random Forest. \code{\link[randomForest]{randomForest}} in \code{randomForest} package
#' \item tree: Classificatin Tree. \code{\link[tree]{tree}} in \code{tree} package
#' \item self: You can use your self-defined function. You need to pass your self-defined function via the "func" parameter.
#' }
#' @param func If you set "model" to "self", you have to pass your self-defined model function. This function should be able to take "x_train" and "y_train" as the first two inputs to train the model and then take "x_test" as the third input and return the predicted scores of x_test data. For example,
#' \preformatted{library(e1071)
#'
#' predict_model <- function(x_train, y_train, x_test){
#'      data_trainxy<-data.frame(x_train,y_train=as.factor(y_train))
#'      fit_svm<-svm(y_train~.,data=data_trainxy,probability=TRUE)
#'      pred_svm <- predict(fit_svm, x_test, probability=TRUE,decision.values = TRUE)
#'      p_svm=as.data.frame(attr(pred_svm, "probabilities"))$"1"
#'      return(p_svm)
#' }
#'
#' result = get_p_result(data_list=data_list, model=c("self","randomforest"), func=predict_model)}
#'
#'
#' @return the scores predicted by models
#' @export
get_predictions <- function(train_data, test_data, model="randomforest", func=NULL)
{
  train_x = train_data[,-length(test_data)]
  train_y = train_data$class
  test_x = test_data[,-length(test_data)]

  ###SVM
  if(model == "svm"){
    fit_svm<-svm(as.factor(class)~.,data=train_data,probability=TRUE)
    pred <- predict(fit_svm, test_x)
  }

  ###RF
  if(model == "randomforest"){
    # randomforrest
    rf <- randomForest(as.factor(class)~.,data=train_data)
    pred <- predict(rf, test_x, type="class")
  }

  ###Classification
  if(model == "tree"){
    fit_tree = tree(as.factor(class)~ ., data = train_data)
    pred = predict(fit_tree, newdata = test_x, type = 'class')
  }

  ###Self-defined
  if(model == "self"){
    pred = func(train_x, train_y, test_x)
  }
  return(pred)

}


#' @import ggplot2 minpack.lm
#' @title plot and fit the data using inverse power law
#'
#' @param results_calculated average index calculated by ssd function, only input one type (one row out of three) of the three.
#' @param n_train_list the series of training sizes you set to generate *results_calculated*. This parameter mush match with *results_calculated*.
#' @param mode the *mode* parameter you set to generate *results_calculated*.
#' @param index the *index* parameter you set to generate *results_calculated*.
#' @param model the *model* parameter you set to generate *results_calculated*.
#'
#' @return Return the plot
#' @export
plot_fit_ipl <- function(results_calculated, n_train_list, mode, index, model){

  n = n_train_list
  nn=seq(1,max(n)+50,1)
  if(mode=="pilot"){group2=c("synthetic-data")}
  if(mode=="true"){group2=c("true-data")}

  newx = data.frame(n=nn)

  # fit_index = nls(results_calculated ~ a*n^(-b), start = list(a=0.5, b=0.1))
  fit_index = nlsLM(results_calculated ~ a*n^(-b)+c, start = list(a=0.5, b=0.1, c=0.1),
                    control = nls.lm.control(maxiter = 1000))

  pred_index = predict(fit_index,newx)

  size <- Class <- NULL

  index_df <-data.frame(
    index = c(results_calculated),
    size = rep(n,1),
    Class = rep(group2,each=length(results_calculated))
  )
  index_pred_df <-data.frame(
    index = c(pred_index),
    nn = rep(nn,1),
    Class = rep(group2,each=length(nn))
  )

  index_min = min(c(results_calculated))
  index_max = max(c(results_calculated))
  delta_max_min = index_max-index_min

  result_plot = plot(ggplot(data = index_df, aes(x = size, y = index, group=Class,color=Class))+
         geom_point()+
         geom_line(data = index_pred_df, aes(x=nn, y=index, group=Class,color=Class), na.rm=TRUE)+
         labs(x = "Size of training data",y=index,title = paste("Predicted Curve:",model))+
         theme_bw()+theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
         ylim(max(index_min-delta_max_min*0.2,0),index_max+delta_max_min*0.2))

  return(result_plot)
}
