#' Read metadata file
#'
#' @param path path to meta data file
#' @export
#'
read.metadata.table <- function(path) {
  if(!file.exists(path)){
    stop(paste0("file not exists: ",path))
  }
  meta_data = read.table(
    meta_data_path,
    sep = "\t",
    header = T,
    quote = "",
    na.strings = "N/A"
  )
  return(meta_data)
}


#' Get sample ids
#'
#' @param dataset dataset number
#' @param meta_data meta data table
#' @export
#'
get.sample.ids <- function(dataset, meta_data) {
  return(as.vector(unique(meta_data[meta_data$dataset %in% dataset, 2])))
}


#' Get sample and run ids
#'
#' @param dataset dataset number
#' @param meta_data meta data table
#' @export
#'
get.run.ids <- function(dataset, meta_data) {
  run_id = meta_data[meta_data$dataset %in% dataset, 2:3]
  colnames(run_id) = c("sample_id", "run_id")
  return(run_id)
}


get_sample_group <- function(sub_sample){
  new_sample_group = as.character(sub_sample$new.sample.group[1])
  new_sample_group = gsub(" ", ".", new_sample_group)
  sample_group = strsplit(new_sample_group, "\\+")[[1]]
  return(sample_group)
}

#' combine all features from coloumn new.sample.group to condition
#'
#' @param dataset dataset number
#' @param meta_data meta data table
#' @return dataframe of metadata including requested dataset
#' @export
#'
get.metadata <- function(sub_sample) {
  sample_group = get_sample_group(sub_sample)
  
  # check if we have feature combinations
  if (length(sample_group) > 1) {
    new_condition = do.call(paste, c(sub_sample[, sample_group], list(sep =
                                                                        "&")))
  } else{
    new_condition = sub_sample[, sample_group]
  }
  
  # add feature combination as main condition to dataset
  sub_sample$condition = as.factor(new_condition)
  
  return(sub_sample)
}


#' Get transcript counts of sample list
#'
#' @param sample sample meta data
#' @export
#'
get.sample.counts <- function(sample) {
  gsm = sample["sample"]
  
  if (file.exists(sample["path"])) {
    
    feature_counts <-  read.table(sample["path"], header = TRUE)
    colnames(feature_counts)[7] = c("Counts")
    tmp = data.frame(feature_counts$Geneid,
                     as.integer(feature_counts$Counts))
    colnames(tmp) <- c("Geneid", gsm)
    return(tmp)
  }else{
      stop(paste0("file not exists: ",sample["path"]))
  }
}

#' Get path of transcript count files
#'
#' @param sample sample meta data
#' @export
#'
get.path <- function(sample) {
  path = paste(
    sample["path"],
    sample["dataset"],
    sample["sample"],
    "/summary/",
    sample["file"],
    sep = "/"
  )
  return(path)
}


#' Get transcript count matrix, conditions and feature length of requested dataset
#'
#' @param dataset dataset numbers
#' @param meta_data meta data table
#' @param working_dir base folder including all datasets
#' @param file name of transcript count file
#' @param transcript_selection transcript selection pattern to grep specific transcripts
#' @param formula formula to select conditions
#' @param filter filter data
#' @param filter_column column to grep by filter
#' @param pseudo_count add pseudo count to feature counts
#' @param min_counts minimum sum of counts per feature
#' @export
#'
get.data <-
  function(datasets,
           meta_data,
           working_dir,
           file = "counts.1.tsv",
           formula = as.formula("~condition"),
           transcript_selection = NA,
           filter = NA,
           filter_column = "disease",
           pseudo_count = 0,
           min_counts = 1) {

    sub_samples = meta_data[meta_data$dataset == datasets, ] 
    
    #filter sub_sample
    if (!is.na(filter)) {
      sub_samples = filter.meta.data(sub_samples, filter, filter_column)
    }
    
    #summary_json = read_summary_json(working_dir,datasets,sub_samples$sample[1])
    
    # extract only one meta data line per sample, if more than one run exists
    sub_samples = sub_samples[!duplicated(sub_samples$sample), ]
    
    sub_samples$path = working_dir
    sub_samples$file = file
    sub_samples$path = apply(sub_samples, 1, get.path)
    tmp = apply(sub_samples, 1, get.sample.counts)
    
    merge.all <- function(x, y) {
      merge(x, y, all = TRUE, by = "Geneid")
    }
    
    count_matrix = Reduce(merge.all, tmp)
    
    # filter by feature id
    if(!is.na(transcript_selection)){
      count_matrix = count_matrix[grep(transcript_selection, count_matrix$Geneid), ]
    }
   
    rownames(count_matrix) <- count_matrix$Geneid
    count_matrix[, 1] <- NULL
    
    #filter by gene count
    count_matrix = filter.count.matrix(count_matrix, min_counts)
    
    # add pseudo count
    count_matrix = count_matrix + pseudo_count
    
    # extract feature length
    tmp = read.table(sub_samples$path[1], header = T)
    feature.length = as.data.frame(tmp$Length)
    rownames(feature.length) = tmp$Geneid
    colnames(feature.length) = "length"
    feature.length = feature.length[rownames(count_matrix), ]
    
    # extract meta data
    dataset_metadata = get.metadata(datasets, sub_samples)
    conditions_out = get.conditons(dataset_metadata, formula)
    condition = conditions_out$condition
    formula = conditions_out$formula
      
    # only select conditions inlcuded in the count matrix after filtering
    conditions = condition[colnames(count_matrix), ,drop = FALSE]
    

    return(
      list(
        count.matrix = count_matrix,
        feature.length = feature.length,
        conditions = conditions,
        formula=as.formula(formula),
        sample_group = get_sample_group(sub_samples)
        
      )
    )
  }

filter.count.matrix <- function(count_matrix, min_counts) {
  filtered = count_matrix[rowSums(count_matrix) > min_counts, ]
  
  if (nrow(filtered) == 0) {
    stop(paste("No features left with counts >", min_counts))
  }
  
  return(filtered)
}


#' Filter dataset based on selected column and filter string.
#' The function returns all samples of a dataset if one sample
#' fulfils the filter criteria
#'
#' @param meta_data meta data
#' @param filter string to grep in column
#' @param column column to select
#' @return data frame of filtered samples
#' @export
#' 
filter.meta.data <- function(meta_data, filter, column) {
  filtered = meta_data[grep(filter, meta_data[, column]), ]
  ids = unique(filtered$dataset)
  filtered = meta_data[meta_data$dataset %in% ids,]
  
  if (nrow(filtered) == 0) {
    stop(paste("No samples left for filter criteria:", filter))
  }
  
  return(filtered)
}


#' Replace all <NA> values with new factor value of selected condition
#'
#' @param coldata dataframe with condition values
#' @param condition condition to select 
#' @param new_value replace <NA> with new_value 
#' @return condition dataframe 
#' @export
#' 
replace.na <- function(coldata, condition, new_value){
  if (condition %in% colnames(coldata)) {
    if(any(is.na(coldata[[condition]]))==T){
      levels(coldata[[condition]]) =  as.factor(c(levels(coldata[[condition]]),new_value))
      coldata[is.na(coldata[[condition]]), condition] <- new_value
      coldata[,condition] =  as.factor(coldata[,condition])
    }
  }
  return(coldata)
}

#' Extract all conditions out of meta data and return a dataframe
#' compatible with deSeq2. Check if 
#'
#' @param meta_data meta data table
#' @param formula formula
#' @export
#' 
get.conditons <- function(dataset_metadata, formula) {
  # extract only one meta data line per sample, if more than one run exists
  meta_data = dataset_metadata[!duplicated(dataset_metadata$sample), ]
  formula_list = all.vars(formula)
  coldata = data.frame(meta_data[, formula_list])
  rownames(coldata) = meta_data$sample
  colnames(coldata) = formula_list
  
  # replace N/A of condition row
  coldata = replace.na(coldata,"condition","N/A")

  # keep only condtions without N/A
  names = colnames(coldata[!apply(coldata,2,function (x){
    any(is.na(x) | x == "")
  })])
  coldata = coldata[,names,drop = FALSE]
 
  # find the quantiles and group the ages based on the quantile boarders 
  if ("age" %in% colnames(coldata)) {
    if(!any(is.na(coldata$age))){
      quant = quantile(coldata$age,na.rm=T)
      quant = unique(quant)
      if(length(quant)==1){
        quant = c(quant,quant)
      }
      quant[1] = quant[1]-1
      coldata$age = cut(coldata$age, breaks = quant)
    }
  }

  # keep only condtions with different values
  names = colnames(coldata[apply(coldata,2,function (x){
    length(unique(x))>1
  })])
  coldata = coldata[,names,drop = FALSE]
  
  # drop levels, if some age bins are not used
  coldata = droplevels(coldata)
  
  
  # build new formula
  formula = create.new.formula(names)
  
  # find collinearity of columns 
  names = find.collinearity(coldata, formula)
  coldata = coldata[,names,drop = FALSE]
  formula = create.new.formula(names)
  
  return(list(conditions=coldata,formula=formula))
}

#' @export
#' 
create.new.formula <-function(names){
  # create new formula
  new_formula = paste(names,collapse="+")
  new_formula = paste0("~",new_formula)
  formula = as.formula(new_formula)
  return(formula)
}


#' Builds a model matrix of the condition matrix and the formula. 
#' Tests if the Variance inflation factor of each factor is below threshold. Otherwise, 
#' the factor is collinear to other fator. Remove factor from condition matrix.   
#'
#' @param coldata condition matrix of dataset
#' @param formula formula
#' @param maxVif VIF threshold 
#' @export
#' 
find.collinearity <- function(coldata,formula, maxVif=10){
  
  formula_list = all.vars(formula)
  # can't calculate VIF with only one factor
  if(length(formula_list)==1){
    return(formula_list)
  }
  
  # calculate model matrix
  test_model = model.matrix(as.formula(formula),coldata)
  
  # exclude first intercept column 
  test_model_df = as.data.frame(test_model[,2:length(colnames(test_model))])
  
  # calculate VIF
  out = usdm::vif(test_model_df)
  
  # if VIF > maxVif, condition is collinear 
  collinar = out[out$VIF>maxVif,]$Variables
  
  # only check collinearity for cofactors 
  sub_formula = formula_list[1:(length(formula_list)-1)]
  tmp = unlist(lapply(sub_formula, function(x) pmatch(x,collinar)>0))
  tmp[is.na(tmp)]=FALSE
  if(length(sub_formula[tmp])>0){
    # print(out)
    warning(paste("The following factors are collinear:",paste(sub_formula[tmp],collapse = ",")))
  }
  
  return(c(sub_formula[!tmp],formula_list[length(formula_list)]))
}

read_summary_json <- function(working_dir,dataset,sample){
  path = file.path(paste(working_dir,dataset,sample,"summary/resources.json",sep="/"))
  summary_json = fromJSON(path)
  return(summary_json)
}
