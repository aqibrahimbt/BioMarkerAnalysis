
#' Checks if the minimum sample number is reached and if the two samples are balanced.
#' @param conditions labels of conditions per sample
#' @export
#'
check_input_samples <- function(conditions, lower_limit=5){
  
  min_samples = min(table(conditions))
  
  valid_comparison = FALSE
  balanced = FALSE
  
  # check if we reach the lower sample number limit, otherwise we can't perform a classification
  if(min_samples>=lower_limit){
    valid_comparison = TRUE
  }
  
  if(length(samples_con1)==length(samples_con2)){
    balanced = TRUE
  }
  
  return(list(valid_comparison=valid_comparison, balanced = balanced, min_sample_number = min_samples))
}

# Get list of seeds for recursive feature elimination 
#' @param repeats - For repeated k-fold cross-validation only: the number of complete sets of folds to compute
#' @param resampling - number of cross validation to use
#' @export
set_seeds <- function(repeats, resampling, subset_sizes){
  length.seeds = (repeats*resampling) +1
  n.tune.parameters = length(subset_sizes)+1
  seeds <- vector(mode = "list", length = length.seeds) #length is = (n_repeats*nresampling)+1
  for(i in 1:length.seeds) seeds[[i]]<- sample.int(n=1000, n.tune.parameters) #(number of tuning parameter)
  seeds[[length.seeds]]<-sample.int(1000, 1)#for the last model
  return(seeds)
}


#' Recursive feature elimination (Backwards Feature Selection) incorporating resampling. Returns a class object with selcted features properties
#' http://topepo.github.io/caret/recursive-feature-elimination.html
#' @param countdata - count data
#' @param conditions - conditions to be used for classification
#' @param subset_sizes - subset_sizes for the recursive feature elimination eg: subsetSizes <- c(2, 4, 6, 8). 
#' @param repeats - the cross-validation procedure is repeated n times (2.1 of Algorithm 2)
#' @param number - number of folds for resampling
#' @param method eg: cv or repeatedcv
#' @param metric - metric to use e.g "Accuracy", "RMSE". 
#' @param balanced - specify if the model will build based on balanced input samples (default: FALSE)
#' @param lower_limit - minimum number of samples per class (default: 5) 
#' @import caret
#' @import jsonlite
#' @export
#'
basic_feature_elimination <- function(countdata, conditions , subset_sizes, repeats, number, method, metric, balanced=FALSE, lower_limit = 5){
  
  set.seed(200) # set seed to get always the same seeds for feature elimination
  
  validate = check_input_samples(conditions, lower_limit)
  seeds = set_seeds(repeats, number, subset_sizes)
  ctrl <- rfeControl(functions = rfFuncs,
                     method = method, 
                     number = number, 
                     repeats = repeats, 
                     verbose = FALSE, 
                     returnResamp = "all", 
                     seeds = seeds)
  
  if(!validate$valid_comparison){
    log_warning(paste0("Cannot run since sample class is less than ",lower_limit))
    return(NULL)
  }else if(balanced){
    log_running("Stratified sampling")
    rfProfile <- rfe(countdata, conditions, sizes = subset_sizes, sampsize=c(validate$min_sample_number,validate$min_sample_number), strata=conditions, rfeControl = ctrl, metric = metric)
  }else{
    log_running("Without stratified sampling")
    rfProfile <- rfe(countdata, conditions, sizes = subset_sizes, rfeControl = ctrl, metric = metric)
    
  }
  return(rfProfile)
}

#' Recursive feature elimination (Backwards Feature Selection) incorporating resampling. Returns a class object with selcted features properties
#' http://topepo.github.io/caret/recursive-feature-elimination.html
#' @param countdata - count data
#' @param conditions - conditions to be used for classification
#' @param subset_sizes - subset_sizes for the recursive feature elimination eg: subsetSizes <- c(2, 4, 6, 8). 
#' @param repeats - the cross-validation procedure is repeated n times (2.1 of Algorithm 2)
#' @param number - number of folds for resampling
#' @param method eg: cv or repeatedcv
#' @param metric - metric to use e.g "Accuracy", "RMSE". 
#' @param balanced - specify if the model will build based on balanced input samples (default: FALSE)
#' @param two_step - use a second rfe two improve feature size
#' @export
#'
recursive_feature_elimination <- function(countdata, conditions, method='repeatedcv', repeats = 10, number = 5,  
                                          metric = 'Accuracy' , balanced=F, two_step=F){

  #Run 1
  subset_sizes = c(2:10,seq(from=15,to=100,by = 5),seq(from=110,200,by = 10),seq(from=220,ncol(countdata),by = 20))
  rfProfile <- basic_feature_elimination(countdata, conditions, subset_sizes, repeats, number, method, metric, balanced)
  profiles <- list(as.list(rfProfile))
  #Run 2
  if(two_step){
    subset_sizes_2 <-c(seq(from=max(2,rfProfile$bestSubset-15),to=min(ncol(countdata),rfProfile$bestSubset+15),by = 1))
    rfProfile_2 <- basic_feature_elimination(countdata, conditions, subset_sizes_2, repeats, number, method, metric, balanced)
    profiles[[2]] <- as.list(rfProfile_2)
  }
 
  return(profiles)
}



#'Aggregates the results from the random forest profile from recursive feature elimination step
#'
#'@param profiles - profiles from feature elimination of the dataset
#'@export
#'
aggregate_profiles <- function(profiles){
  aggregate_data = list()
  for (i in 1:length(profiles)){
    rfProfile <- profiles[[i]]
    res<-cbind(rfProfile$resample[,c(1,2)], data.frame(do.call('rbind', strsplit(as.character(rfProfile$resample$Resample),'.',fixed=TRUE))))
    names(res)<-c("Variables","Accuracy","Fold","Repeat")
    aggdata <-aggregate(Accuracy~Variables+Repeat, res, mean)
    aggregate_data[[i]] <-aggregate(Accuracy~Variables,aggdata, FUN = function(x) quantile(x,probs=c(0,0.25,0.5,0.75,1)))
  }
  return(do.call(rbind, aggregate_data))
}


#'Generates log scale of aggregated accuracy data of recursive feature elimination
#'
#'@param opt_data - reduced/optimized feature list
#'@param agg_data - aggregated accuracy data of recursive feature elimination
#'@param rfProfile - random forest model/profile
#'@export
#'
feature_logscale <- function(opt_data, agg_data, rfProfile){
  logScale<-c()
  rfProfile <- rfProfile
  curX<-ncol(as.data.frame(opt_data))
  while(curX>1){
    curDiffs<-agg_data[,1]-curX
    curInd<-which(abs(curDiffs)==min(abs(curDiffs)))[1]
    curElem<-agg_data[curInd,1]
    logScale<-c(logScale,curElem)
    curX<-round(curX/2.0)
  }
  logScale<-sort(unique(c(logScale,rfProfile$bestSubset)))
  data<-agg_data[agg_data[,1]%in%logScale,]
  featNewJson<-list()
  for(i in 1:nrow(data)){
    featNewJson[[i]]<-list("number_of_features"=unbox(data[i,"Variables"]),
                           "quantiles"=list(
                             "0.00"=unbox(data[i,'Accuracy'][1]),
                             "0.25"=unbox(data[i,'Accuracy'][2]), 
                             "0.50"=unbox(data[i,'Accuracy'][3]), 
                             "0.75"=unbox(data[i,'Accuracy'][4]), 
                             "1.00"=unbox(data[i,'Accuracy'][5])))
    if(data[i,"Variables"]==rfProfile$bestSubset){
      featNewJson[[i]][["optimal_number_of_features"]]=unbox(TRUE)
    }else{
      featNewJson[[i]][["optimal_number_of_features"]]=unbox(FALSE)
    }
  }

  return(featNewJson)
}

#' set number of trees based on sample size, due to performance the maximum tree number for
#' classification with more than 100 samples is set to 10,000 otherwise to 100,000
#' @param ntree number of trees
#' @export
set_number_of_trees <- function(ntree){
  if(is.null(ntree)){
    if(length(conditions)>100){
      ntree = 10000
    }else{
      ntree = 100000
    }
  }
  return(ntree)
}

#' Optimize random forest model with cross validation
#' optimization parameter: mtry 
#'
#'@param count_data count data, used for classification
#'@param conditions classification labels 
#'@param model_type eg: 'rf', 'ranger'
#'@param number Either the number of folds or number of resampling iterations
#'@param repeats For repeated k-fold cross-validation only: the number of complete sets of folds to compute
#'@param balanced specify if the model will build based on balanced input samples (default: FALSE)
#'@param ntree number of trees, default 100,000 if sample size >100, set to 10,000
#'@import caret
#'@export
#'

train_model <- function(count_data, conditions, model_type, repeats = 10, number = 5,  balanced=FALSE, ntree=NULL) {

  ntree =  set_number_of_trees(ntree)
  
  validate = check_input_samples(conditions)
  control <- trainControl(method="repeatedcv", number = number, repeats = repeats,returnResamp="all")
  
  if(balanced){
    
    log_running('Running Stratified Sampling')
    rf.model <- train(x = count_data,
                      y = as.factor(conditions),
                      method = model_type, 
                      metric = "Accuracy",
                      trControl = control,
                      verbose = F,
                      strata=as.factor(conditions),
                      sampsize=c(validate$min_sample_number,validate$min_sample_number),
                      ntree = ntree,
                      importance=TRUE,
                      proximity=TRUE,
                      )
    
  }else{
    
    log_running('Running Without Stratified Sampling')
    rf.model <- train(x = count_data,
                      y = as.factor(conditions),
                      method = model_type, 
                      metric = "Accuracy",
                      trControl = control,
                      tuneGrid = data.frame(mtry = 6),
                      savePredictions = TRUE,
                      verbose = F)
  }
  
  return(rf.model)
}

#' Calculate the final random forest model based on the reduced/optimized feature set
#' The final model can be balanced or un-balanced
#' 
#' @param countdata count data matrix
#' @param conditions labels of count data matrix
#' @param balanced specify if the model will build based on balanced input samples (default: FALSE)
#' @param ntree number of trees to use for random forest model
#' @param lower_limit - minimum number of samples per class (default: 5) 
#' @param ntree number of trees, default 100,000 if sample size >100, set to 10,000
#' @import randomForest
#' @export
#' 
get_final_model <- function(countdata, conditions, balanced=FALSE, ntree=NULL, lower_limit=5){
  
  ntree =  set_number_of_trees(ntree)
  mtry = get_mtry(ncol(countdata))
  print(mtry)
  validate = check_input_samples(conditions, lower_limit)
  
  if(!validate$valid_comparison){
    log_warning(paste0("Cannot run since sample class is less than ",lower_limit))
    return(NULL)
  }else if(balanced){
    log_running(paste("balanced sampling with number of samples in each class:", validate$min_sample_number))
    rf.model <- randomForest(countdata,
                             conditions,
                             mtry=mtry,
                             ntree=ntree,
                             importance=TRUE,
                             proximity=TRUE,
                             strata=conditions,
                             sampsize=c(validate$min_sample_number,validate$min_sample_number),
                             keep.inbag = T)
  }else{
    log_running("Unbalanced sampling")
    rf.model <- randomForest(countdata,
                             conditions,
                             mtry=mtry,
                             ntree=ntree,
                             importance=TRUE,
                             proximity=TRUE,
                             keep.inbag = T)
  }
  
  return(rf.model)
}

#' Randomly split count data in train and test set. test = 1/fold; train = 1 - test
#' 
#' @param fold - number of folds for resampling
#' @param countdata - count data matrix
#' @param conditions - labels of count data matrix
#' @param seed - seed to make sampling reproducable 
#' @exprot
#' 
split_train_test <- function(fold, countdata,conditions, seed){
  set.seed(seed)
  sample_number = nrow(countdata)
  partition = trunc(sample_number/fold)

  test_index <- sample(1:sample_number,partition,replace = F)
  train <- countdata[ -test_index,]
  test <- countdata[test_index,]
  
  condition_train <- conditions[-test_index]
  condition_test <- conditions[test_index]
  
  return(list("train"=train,"train_condition"=condition_train,"test"=test,"test_condition"=condition_test))
}

#' Number of variables randomly sampled as candidates at each split.
#'
#' @param nfeature number of maximal features/variables (genes) 
#' @export
#'
get_mtry <- function(nfeature){
  mtry=floor(max(sqrt(nfeature),1))
  return(mtry)
}


#' Get cross validation error based on all models used by cross validation
#' 
#' @param cv cross validation object, including for all repeats predictions, conditions and feature importance values 
#' @export
#' 
cv_error <- function(cv){
  
  # get condition names
  prediction.names = colnames(cv$predictions[[1]])
  
  # unlist predictions based on first condition
  cv.pred = unlist(lapply(cv$predictions, function(x) x[,prediction.names[1]]))
  cv.condition = unlist(cv$conditions)
  
  # label predictions
  cv.pred.label = ifelse(cv.pred > 0.5, prediction.names[1], prediction.names[2])
  
  # create a table with real conditions and predictions
  cv.table = as.data.frame(cv.pred.label, stringsAsFactors = F)
  colnames(cv.table ) = c("prediction")
  cv.table ["condition"] = as.character(cv.condition)
  cv.table ["sample"] = names(cv.pred.label)

  # summarize number of samples and number of correct predictions
  cv.table.errorCV  = cv.table  %>%
    group_by(sample) %>%
    add_count(sample) %>% #count number of samples
    add_count(prediction==condition) %>% #count number of samples where prediction == condition
    summarise(sample_counts=max(n),
              correct_prediction=max(nn))
  
  cv.table.errorCV["cross_validation_error"] = 1-cv.table.errorCV["correct_prediction"]/cv.table.errorCV["sample_counts"]
  return(cv.table.errorCV)
}

#' Get feature importance of all models used by cross validation
#' 
#' @param cv cross validation object, including for all repeats predictions, conditions and feature importance values 
#' @export
#' 
cv_feature_importance <- function(cv){
  cv_inportance = as.data.frame(do.call(rbind, cv$importance))
  cv_inportance["feature_names"] =  rownames(cv_inportance)
  
  cv_importance_tmp = cv_inportance %>% 
    group_by(feature_names) %>%
    summarise(quantile = list(quantile(MeanDecreaseGini)))
  
  cv_importance_quant = as.data.frame(cv_importance_tmp$quantile)
  colnames(cv_importance_quant) = cv_importance_tmp$feature_names
  cv_importance_quant_df = as.data.frame(t(cv_importance_quant))
  
  return(cv_importance_quant_df)
}

#' Write data to json file 
#' 
#' @param dataset - data set name 
#' @param data - arbitray data set
#' @param description - description of input data set
#' @param output_dir - path to output directory 
#' @param file_name - json file name, without extension
#' @import jsonlite
#' @export
#' 
write_to_json <- function(dataset, data, description, output_dir, file_name){
  meta = list(dataset = dataset,  description = description)
  file_name = generate_filename(file_name, output_dir, '.json')
  write_json(list('metadata' = meta, 'data' = data),file_name, pretty = TRUE)
}

#' Random forest model cross validation with n folds and m repeats 
#'
#' @param countdata count data matrix
#' @param conditions labels of count data matrix
#' @param repeats - the cross-validation procedure is repeated n times
#' @param fold - number of folds for resampling
#' @param ntree number of trees, default 100,000 if sample size >100, set to 10,000
#' @import randomForest
#' @export
#' 
cross_validation <- function(fold, repeats, countdata, conditions, ntree=NULL){
  
  ntree =  set_number_of_trees(ntree)
  
  test_predictions <-list()
  test_conditions <-list()
  test_importtance <- list()
  
  mtry = get_mtry(ncol(countdata))
  
  for(i in 1:repeats){
    split = split_train_test(fold,countdata,conditions,i+22)
    rf.model.cv<-randomForest(split$train,split$train_condition,mtry=mtry,importance=TRUE,proximity=TRUE, ntree=ntree)
    test_predictions[[i]]<-predict(rf.model.cv,split$test,type='prob')
    test_conditions[[i]]<-split$test_condition
    test_importtance[[i]] <- rf.model.cv$importance
  }
  return(list("conditions"=test_conditions, "predictions"=test_predictions, "importance"=test_importtance))
}

#' Receiver Operating Characteristic curve (ROC) and area under this curve (AUC)
#' 
#' @param rf_model random forest model 
#' @import ROCR
#' @export
roc <- function(rf_model){
  rf_predictions <- predict(rf_model, type = 'prob')
  pred<-prediction(rf_predictions[,2], conditions)
  roc<-performance(pred, 'tpr', 'fpr')
  auc<-performance(pred, 'auc')
  x <- unlist(slot(roc, "x.values"))
  y <- unlist(slot(roc, "y.values"))
  auc_val = unlist(slot(auc, "y.values"))
  
  return(list("FP"=x,"TP"=y,"AUC"=auc_val))
}

#' Out-of-bag error(OOB) obtained for each independent tree that forms the forest
#' 
#' @param rf_model random forest model 
#' @export
oob <- function(rf_model){
  oob<-rf_model$err.rate
  oob_df = as.data.frame(oob)
  oob_df["ntree"] = rownames(oob_df)
  return(oob_df)
}



#' Preprocessing step to change data into the right format and splits into training and test datasets
#' 
#' @param count_matrix count matrix
#' @param conditions conditions 
#' @param pct - percentage for test and training data
#' @export
preProcessing <- function(count_matrix, conditions, pct){
  input <- count_matrix
  input <- as.data.frame(cbind(conditions, input))
  colnames(input)[1] <- 'Conditions'
  smp_siz = floor(pct*nrow(input))
  train_ind = sample(seq_len(nrow(input)),size = smp_siz)

  train = as.data.frame(input[train_ind,])
  test= as.data.frame(input[-train_ind,])

  return(list("train"=train, "test"=test))
}


#' Kmeans clustering step
#' 
#' @param train training dataset
#' @param no_clusters number of clusters  
#' @export
kmeansclustering <- function(train, no_clusters){
  train = t(train[, -1])
  kclusters = kmeans(train, no_clusters)
  
  return(kclusters)
}



#' SVM Scoring step ( runs svm on the clusters and scores clusters)
#' 
#' @param train training dataset
#' @param kclusters required number of clusters  
#' @export
svm_scoring <- function(train, kclusters){
  cluster_size = length(kclusters$size)
  
  svm_result <- data.frame(matrix(nrow = cluster_size, ncol = 3))
  colnames(svm_result) <- c("cluster", "accuracy", "count_genes")
  
  #train_dump = t(train[, -1])
  
  for (i in 1:cluster_size){
    # cat('Prediction for Cluster', i, '\n')
    
    cluster <- as.data.frame(train[, kclusters$cluster == i])
    cluster <- cbind(train$Conditions, cluster)
    colnames(cluster)[1] <- 'Conditions'
    
    #Basic Model
    svm_model= svm(Conditions~.,data=cluster, type='C-classification', cross=3)
    
    svm_result[i, 1] <- i
    svm_result[i, 2] <- svm_model$tot.accuracy
    svm_result[i, 3] <- ncol(cluster)
  }
  
  return(svm_result)
}


#' RCE step to remove clusters with the lowest accuracy. 
#' 
#' @param svm_result result from svm scoring step
#' @param kclusters required number of clusters  
#' @param train training dataset 
#' @export
rce_step <-function(svm_result, kclusters, train){
  
  subset_res_low <-subset(svm_result, svm_result$accuracy <= quantile(svm_result$accuracy, prob =0.20))
  cluster_nos <- as.vector(subset_res_low$cluster)
  
  train = lapply(cluster_nos, getClusterGenes, kclusters=kclusters, train=train)
  train = as.data.frame(train[length(cluster_nos)])
  
  return(list("train" = train, "cluster_k" = length(cluster_nos)))
}



#' Retrieves gene list given cluster number 
#' 
#' @param cluster_no cluster number
#' @param train training dataset 
#' @param kclusters training dataset 
#' @export
getClusterGenes <- function(cluster_no, train, kclusters){
  
  cluster <- as.data.frame(train[, kclusters$cluster == cluster_no])
  
  col_names = as.list(colnames(cluster))
  myvars <- names(train) %in% col_names
  train <- train[!myvars]
  
  return(train)
}



#' Returns a new dataframe with selected genes
#' 
#' @param train training dataset
#' @param no_clusters number of clusters 
#' @param required_no_clusters required number of clusters
#' @export
rce_svm <- function(train, no_clusters, required_no_clusters){
  
  while (no_clusters > required_no_clusters){
    kclusters = kmeansclustering(train, no_clusters)
    svm_result = svm_scoring(train, kclusters)
    new_train = rce_step(svm_result, kclusters, train)
    
    train = new_train$train
    no_clusters = no_clusters - new_train$cluster_k
    cat("No of clusters ", no_clusters, "\n")
    cat("No of Variables ", ncol(train), "\n")
  }
  return(train[, -1])
}





