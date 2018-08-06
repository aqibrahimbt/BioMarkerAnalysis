
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
#' @param balanced - True or False (For stratified sampling)
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
#' @param balanced - True or False (For stratified sampling)
#' @param two_step - use a second rfe two improve feature size
#' @export
#'
recursive_feature_elimination <- function(countdata, conditions, method='repeatedcv', repeats = 10, number = 5,  metric = 'Accuracy' , balanced=F, two_step=F){

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



#'Aggregates the results from the rfprofiles from feature elimination step
#'@param rfProfile - profiles from feature elimination of the dataset
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


#'Generates log scale of the data
#'
#'@param opt_data   Optimised data
#'@param agg_data - Aggregated data from rfprofile
#'@param rfProfile
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

#' Train model with cross validation 
#'
#'@param count_data count data, used for classification
#'@param conditions classification labels 
#'@param model_type eg: 'rf', 'ranger'
#'@param number Either the number of folds or number of resampling iterations
#'@param repeats For repeated k-fold cross-validation only: the number of complete sets of folds to compute
#'@param balanced Balancing of imput samples 
#'@param ntree number of trees, default 100,000 if sample size >100, set to 10,000
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


get_final_model <- function(countdata, conditions, balanced=FALSE, ntree=NULL, lower_limit=5){
  
  ntree =  set_number_of_trees(ntree)
  mtry = get_mtry(ncol(countdata))
  validate = check_input_samples(conditions, lower_limit)
  
  if(!validate$valid_comparison){
    log_warning(paste0("Cannot run since sample class is less than ",lower_limit))
    return(NULL)
  }else if(balanced){
    log_warning(paste("balanced sampling with number of samples in each class:", validate$min_sample_number))
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
    log_warning("Unbalanced sampling")
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

#' randomly split count data in train and test set. test = 1/fold; train = 1 - test
#' @param fold fold
#' @param countdata count data matrix
#' @param conditions labels of count data matrix
#' @param seed seed to make sampling reproducable 
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

#'Number of variables randomly sampled as candidates at each split.
#'@param nfeature number of maximal features/variables (genes) 
#'@export
get_mtry <- function(nfeature){
  mtry=floor(max(sqrt(nfeature),1))
  return(mtry)
}
  
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

write_to_json <- function(datasets, data, description, output_dirm, file_name){
  meta = list(dataset = datasets,  description = description)
  file_name = generate_filename(file_name, output_dir, '.json')
  write_json(list('metadata' = meta, 'data' = data),file_name, pretty = TRUE)
}

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

roc <- function(final_rf_model){
  rf_predictions <- predict(final_rf_model, type = 'prob')
  pred<-prediction(rf_predictions[,2], conditions)
  roc<-performance(pred, 'tpr', 'fpr')
  auc<-performance(pred, 'auc')
  x <- unlist(slot(roc, "x.values"))
  y <- unlist(slot(roc, "y.values"))
  auc_val = unlist(slot(auc, "y.values"))
  
  return(list("FP"=x,"TP"=y,"AUC"=auc_val))
}

oob <- function(final_rf_model){
  oob<-final_rf_model$err.rate
  oob_df = as.data.frame(oob)
  oob_df["ntree"] = rownames(oob_df)
  return(oob_df)
}

