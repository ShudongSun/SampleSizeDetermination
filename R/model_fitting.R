#' @title Sample Size Determination
#' @description build the model and determine the sample size based on elbow method of ARI and AMI.
#'
#' @param x a data frame or a matrix of predictors to be fitted.
#' @param y A response vector.
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
#' }
#'
#' @param n_train_list a series of numbers which represent the generated training data size for each class.
#' @param n_test the generated test data size for each class.
#' @param mode The default value is "pilot", which means that we will see the input data(x and y) as the pilot data(usually small number), and use them to generate simulated data (corresponding to 'n_train_list'), and then get the curve and elbow point; "true" means that we will sample the data corresponding to 'n_train_list' from input data (x and y), and then get the curve and elbow point.
#' @param num_repeat The number of times to repeat for each training data size in order to get a average value of index to reduce the effect of randomness. The default value is 30.
#' @param print_progress_bar Whether to print the progress bar and elapsed time. The default value is TRUE.
#' @param n.cores The number of nodes to be forked when using multi-core parallel computing. If not being set(n.cores=NULL), \code{n.cores <- parallel::detectCores() - 1} would be used.
#'
#'
#' @return Return the average classification errors/ARI/AMI of the repeat results for the selected mode.
#' @import mclust MASS aricode progress parallel foreach doParallel ggplot2 minpack.lm
#' @importFrom randomForest randomForest
#' @importFrom stats cov nls pnorm predict qnorm quantile
#'
#' @examples
#' data_pmbc <- read.csv(system.file("extdata", "data_pmbc_24pc.csv", package = "SSD"),row.names=1)
#' data = data_pmbc
#'
#' ### use pilot data:
#' p = 15
#' table(data$phenoid)
#' num_class = length(table(data$phenoid))
#' num_PC = length(colnames(data))-1
#' for(i in 1:num_class){
#'   class_i_ids = which(data$phenoid == names(table(data$phenoid))[i])
#'   pilot_i_ids = sample(class_i_ids, p)
#'   pilot_i_data = data[pilot_i_ids,]
#'   if(i == 1){
#'     pilot_data = pilot_i_data
#'   }else{
#'     pilot_data = rbind(pilot_data, pilot_i_data)
#'   }
#' }
#' x_pilot = pilot_data[,-length(pilot_data)]
#' y_pilot = pilot_data[,length(pilot_data)]
#' table(pilot_data$phenoid)
#'
#' result_pilot = ssd(x_pilot, y_pilot)
#'
#' ### use true data:
#'
#' x_true = data[,-length(data)]
#' y_true = data[,length(data)]
#'
#' result_true = ssd(x_true, y_true, mode="true")
#'
#' @export
ssd <- function(x, y, model="randomforest", func=NULL, index="ARI", n_train_list = seq(from=30,to=600,by=30), n_test = 300, mode="pilot", num_repeat=30, print_progress_bar=TRUE, n.cores=NULL) {

  if(!model %in% c( "svm", "randomforest", "tree","self")){
    stop('model \'',model, '\' cannot be found')
  }

  if(!index %in% c("ARI","AMI","classification error")){
    stop('index \'',index, '\' cannot be found')
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
    errors_syn <- vector(length=length(n_train_list))
    ARIs_syn <- vector(length=length(n_train_list))
    AMIs_syn <- vector(length=length(n_train_list))
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
        cm = table(test_data[,ncol(test_data)], pred)
        errors_temp = misclass_err(test_data[,ncol(test_data)], pred)
        ARIs_temp = adjustedRandIndex(test_data[,ncol(test_data)], pred)
        AMIs_temp = AMI(test_data[,ncol(test_data)], pred)

        return(cbind(errors_temp,ARIs_temp,AMIs_temp))

      }

      errors_syn[j] = mean(result_temp[,"errors_temp"])
      ARIs_syn[j] = mean(result_temp[,"ARIs_temp"])
      AMIs_syn[j] = mean(result_temp[,"AMIs_temp"])

      j=j+1
    }
    result = rbind(errors_syn,ARIs_syn,AMIs_syn)

  }else if(mode=="true"){
    data = cbind(x,y)
    names(data)[length(names(data))] = "class"
    num_class = length(table(y))
    num_PC = length(colnames(x))

    n_train_list = n_train_list
    errors_true <- vector(length=length(n_train_list))
    ARIs_true <- vector(length=length(n_train_list))
    AMIs_true <- vector(length=length(n_train_list))
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

        for(i in 1:num_class){
          class_i_ids = which(data$class == names(table(data$class))[i])
          if(length(class_i_ids)>=(n_train+n_test)){
            train_test_i_ids = sample(class_i_ids, (n_train+n_test))
            train_i_data = data[train_test_i_ids[1:n_train],]
            test_i_data = data[train_test_i_ids[(n_train+1):(n_train+n_test)],]
          }else{
            r = n_train/((n_train+n_test))
            train_i_ids = sample(class_i_ids[1:ceiling(length(class_i_ids)*r)],n_train,replace=TRUE)
            test_i_ids = sample(class_i_ids[ceiling((length(class_i_ids)*r)+1):length(class_i_ids)],n_test,replace=TRUE)
            train_i_data = data[train_i_ids,]
            test_i_data = data[test_i_ids,]
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
        cm = table(test_data[,ncol(test_data)], pred)

        errors_temp = misclass_err(test_data[,ncol(test_data)], pred)
        ARIs_temp = adjustedRandIndex(test_data[,ncol(test_data)], pred)
        AMIs_temp = AMI(test_data[,ncol(test_data)], pred)

        return(cbind(errors_temp,ARIs_temp,AMIs_temp))

      }

      errors_true[j] = mean(result_temp[,"errors_temp"])
      ARIs_true[j] = mean(result_temp[,"ARIs_temp"])
      AMIs_true[j] = mean(result_temp[,"AMIs_temp"])

      j=j+1
    }
    result = rbind(errors_true,ARIs_true,AMIs_true)
  }

  parallel::stopCluster(cl = my.cluster)

  if(index=="ARI"){
    if(mode=="pilot"){index_used = ARIs_syn}
    if(mode=="true"){index_used = ARIs_true}
  }else if(index=="AMI"){
    if(mode=="pilot"){index_used = AMIs_syn}
    if(mode=="true"){index_used = AMIs_true}
  }else if(index=="classification error"){
    if(mode=="pilot"){index_used = errors_syn}
    if(mode=="true"){index_used = errors_true}
  }

  # print(index_used)


  plot_fit_ipl(index_used, n_train_list, mode, index, model)


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
#' @param index_used average index calculated by ssd function, only input one type (one row out of three) of the three.
#' @param n_train_list the series of training sizes you set to generate *index_used*. This parameter mush match with *index_used*.
#' @param mode the *mode* parameter you set to generate *index_used*.
#' @param index the *index* parameter you set to generate *index_used*.
#' @param model the *model* parameter you set to generate *index_used*.
#'
#' @return Return the plot
#' @export
plot_fit_ipl <- function(index_used, n_train_list, mode, index, model){

  n = n_train_list
  nn=seq(1,max(n)+50,1)
  if(mode=="pilot"){group2=c("synthetic-data")}
  if(mode=="true"){group2=c("true-data")}

  newx = data.frame(n=nn)

  # fit_index = nls(index_used ~ a*n^(-b), start = list(a=0.5, b=0.1))
  fit_index = nlsLM(index_used ~ a*n^(-b)+c, start = list(a=0.5, b=0.1, c=0.1),
                    control = nls.lm.control(maxiter = 1000))

  pred_index = predict(fit_index,newx)

  size <- Class <- NULL

  index_df <-data.frame(
    index = c(index_used),
    size = rep(n,1),
    Class = rep(group2,each=length(index_used))
  )
  index_pred_df <-data.frame(
    index = c(pred_index),
    nn = rep(nn,1),
    Class = rep(group2,each=length(nn))
  )

  index_min = min(c(index_used))
  index_max = max(c(index_used))

  result_plot = plot(ggplot(data = index_df, aes(x = size, y = index, group=Class,color=Class))+
         geom_point()+
         geom_line(data = index_pred_df, aes(x=nn, y=index, group=Class,color=Class), na.rm=TRUE)+
         labs(x = "Size of training data",y=index,title = paste("Predicted Curve:",model))+
         theme_bw()+theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
         ylim(max(index_min-0.02,0),index_max+0.05))

  return(result_plot)
}
