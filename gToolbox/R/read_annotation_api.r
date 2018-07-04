
#' @export
#' @import httr
#'
ignore_ssl <-function(){
  set_config(httr::config(ssl_verifypeer = F))
}

#' @import httr
#' @export
#'
get_json <- function(group,id, max_request=10){
  status = 500
  num_request = 0
  response = NULL
  while(status!=200 & num_request<max_request){
    response <- GET(paste("https://tomcat.genevention.com/annotapi",group,id,sep="/"))
    status = response$status_code
    print(paste0("reconect: ", num_request))
    num_request = num_request+1
  }
  
  if(status!=200){
    response = NULL
  }
  
  return(httr::content(response))
}

#' @export
#'
is_annotated <- function(json){
  status = json$customProperties$annotStatus
  if(status==100 | status==10){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

#' @export
#'
get_sample_ids <- function(json){
  samples = unlist(json$samples)
  return(samples)
}

#' @export
#'
get_value <- function(sample,term){
  out = list(term=sample$annotations[[term]][[1]]$curatedTerm)
  names(out) <- c(gsub(" ",".",term))
  return(out)
}

#' @export
#'
get_metadata_sample <- function(sample_name,dataset){
 
  sample = get_json("samples",sample_name)
  header = c("disease stage","tissue phenotype","phenotype","gender","age","cell line","disease","latency type","time","new sample group","virus","tissue","treatment","disease details","genotype","tissue from cell line","cell type","experimental process","tissue part")
  
  line = unlist(lapply(header, get_value, sample=sample))
  line["sample"] = sample_name
  line["dataset"] = dataset
 
  for(name in header){
    name = gsub(" ",".",name)
    if(!(name %in% names(line))){
      line[name] = NA
    }
  }
  
  tmp = list(line)
  names(tmp) = sample_name
  return(data.frame(t(data.frame(tmp))))
}

#' @import jsonlite
#' @export
#'
get_metadata <- function(dataset){
  
  ignore_ssl()
  
  json = get_json("datasets",dataset)
  
  if(is_annotated(json)){
    samples = get_sample_ids(json)
    output = lapply(samples, get_metadata_sample, dataset = dataset)
    json_meta_data = dplyr::bind_rows(output)
    
    json_meta_data$age = as.numeric(json_meta_data$age)
    
    return(json_meta_data)
  }else{
    stop(paste("Dataset is not not completed:", dataset))
  }

}