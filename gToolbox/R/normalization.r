#' Normalaize the dataset with selected normalization method
#' 
#' @param dataset_data count data, meta data and feature length
#' @param method available methods: rpm, rpkm, deseq
#' @return normalized data 
#' @export
normalize_data <- function(dataset_data, method="rpm"){

  methods = list("rpm"=rpm, "rpkm"=rpkm, "deseq"=deseq_normalization)
  
  if(!method %in% names(methods)){
    stop(paste0("Method ",method, " doesn't exists. Possible methods are: ", paste0( names(methods), collapse = ", ")))
  }
  
  return(methods[[method]](dataset_data))
}

#' Normalaize the dataset with RPM or RPKM normalization method
#' 
#' @param dataset_data count data, meta data and feature length
#' @param rpkm if true, returns rpkm
#' @return rpkm 
#' @export
rpm <- function(dataset_data, rpkm=F){
  k = 1
  if(rpkm){
    k = 1000
  }
  
  l = dataset_data$feature.length
  cS <- colSums(dataset_data$count.matrix) #Total mapped reads per sample 
  rpm = t(dataset_data$count.matrix)/cS
  rpm = t(rpm)
  return(as.data.frame(10^6*k*(rpm/l)))
}

rpkm <- function(dataset_data){
  rpm(dataset_data,T)
}

#' Normalaize the dataset with DESeq normalization methods
#' @param dataset_data including count matrix, conditions and formula
#' @export
#'
deseq_normalization <- function(dataset_data){
  dds <-
    DESeqDataSetFromMatrix(
      countData = dataset_data$count.matrix,
      colData = dataset_data$conditions,
      design = dataset_data$formula
    )
  cds = estimateSizeFactors(dds)
  normalizedCounts <- counts(cds, normalized=TRUE)
  return(as.data.frame(normalizedCounts))
}