
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
    print(paste0("Cannot run since sample class is less than ",lower_limit))
    return(0)
  }else if(balanced){
    print("Running stratified sampling")
    rfProfile <- rfe(countdata, conditions, sizes = subset_sizes, sampsize=c(validate$min_sample_number,validate$min_sample_number), strata=conditions, rfeControl = ctrl, metric = metric)
  }else{
    print("Running without stratified sampling")
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
#'WIP - Work in Progress
#'@param opt_data   Optimised data
#'@param agg_data - Aggregated data from rfprofile
#'@param rfProfile
#'@export
#'
feature_logscale <- function(opt_data, agg_data, rfProfile, output_dir){
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
  meta = list(dataset = datasets,  description = "Log Scale")
  file_name = generate_filename('log_scale', output_dir, '.json')
  
  write_json(list('metadata' = meta, 'data' = featNewJson),file_name, pretty = TRUE)
}


#'@param count_data
#'@param conditions
#'@param model_type eg: 'rf', 'ranger'
#'@param resampling
#'@param repeats
#'@param mtry
#'@param stratified
#'strata refers to which feature to do stratified sampling on.
#'sampsize refers to the size of the bootstrap samples to be taken from each class. These samples will be taken as input
# for each tree.
train_model <- function(count_data, conditions, model_type, resampling, repeats, mtry, stratified=FALSE) {
  
  check = stratified_sampling_check(conditions)
  control <- trainControl(method="repeatedcv", number = resampling, repeats = repeats)
  
  if(model_type == "rf"){
    tuneGrid = expand.grid(.mtry=mtry)
  }else{
    tuneGrid = expand.grid(mtry = mtry, min.node.size = 1, splitrule = 'gini')
  } 
  
  if(all(c(check$stratified, stratified), TRUE)){
    
    print('Running Stratified Sampling')
    rf.model <- train(x = count_data,
                      y = as.factor(conditions),
                      method = model_type, 
                      metric = "Accuracy",
                      trControl = control,
                      verbose = F,
                      tuneGrid = tuneGrid,
                      strata=as.factor(conditions),
                      sampsize = c(check$min_class,check$min_class))
  }else{
    
    print('Running Without Stratified Sampling')
    rf.model <- train(x = count_data,
                      y = as.factor(conditions),
                      method = model_type, 
                      metric = "Accuracy",
                      trControl = control,
                      verbose = F,
                      tuneGrid = tuneGrid)
  }
  
  return(rf.model)
}







