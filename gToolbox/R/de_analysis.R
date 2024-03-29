#' Create DDS object file
#' @param dataset_data including count matrix, conditions and formula
#' @import DESeq2
#' @export
#'
de_result <- function(dataset_data) {
  dds <-
    DESeqDataSetFromMatrix(
      countData = dataset_data$count.matrix,
      colData = dataset_data$conditions,
      design = dataset_data$formula
    )
  dds <- DESeq(dds,betaPrior = FALSE,parallel=FALSE)
  return(dds)
}

#' Count number of identical features in condition_1 and condition_2.
#' If one of the features is NA, count this feature as identical
#' 
#' @param condition_1 list of features of condition 1
#' @param condition_2 list of features of condition 2
#' @export
#'
intersection_count <- function(condition_1,condition_2){
  count = 0
  
  for(i in 1:length(condition_1)){
    if(condition_1[i]==condition_2[i]){
      count = count +1
    }else{
      if(condition_1[i]=="NA" | condition_2[i]=="NA"){
        count = count +1
      }
    }
  }
  
  # check if one of the feature list only consists of NA, set count to 0 
  if(length(condition_1[!condition_1 %in% "NA"])==0 | length(condition_2[!condition_2 %in% "NA"])==0){
    count = 0
  }
  
  return(count)
}

#' To always get the same "fold change direction" we have to set the healthy comparison as first condition
#' return position of comparison
#' @export
get_healthy_position <- function(comparison_with_healthy,condition_matrix,disease_pos){
  healthy_position = 1
  # set healthy sample as condition 1
  if(comparison_with_healthy){
    healthy_position = which(condition_matrix[,disease_pos] == "healthy")
  }else{
    row_len = length(condition_matrix[1,])
    healthy_position = grep("healthy",condition_matrix)
    healthy_position = ceiling(healthy_position / row_len)[1]
  }
  
  if(is.na(healthy_position)){
    healthy_position = 1
  }
  
  return(healthy_position)
}


#' Get Comparisons for all feasible combinations. Uses the concept of degree of freedoms and compares only
#' conditions in the sample group with the least differences. (n-1)
#' @import stringr
#' @export
get_comparisons <- function(conditions, sample_group){
  
  unique_conditions = unique(conditions$condition)
  cond_combn = combn(unique(unique_conditions), 2)
  
  conditions_1 = c()
  conditions_2 = c()
  
  sample_conditions = list()
  
  # check if feature disease is included in sample_group
  disease_pos = which(sample_group == "disease")
  
  for (i in 1:(length(cond_combn) / 2)) {
    # Get number of different sample groups. Add 1 since the character &" will be one less
    count = str_count(as.character(cond_combn[, 1][1]), "&") + 1
    
    # Get the intersection of strings to check similarity.
    splits = strsplit(as.character(cond_combn[, i]), "&")
    
    common = intersection_count(splits[[1]],splits[[2]])
    # common = intersect(splits[[1]],splits[[2]])
    condition_1 = toString(cond_combn[, i][1])
    condition_2 = toString(cond_combn[, i][2])
    
    # check if it is a comparison with healthy
    condition_matrix = do.call(rbind,splits)
    comparison_with_healthy = "healthy" %in% condition_matrix[,disease_pos]

    # identify condition with healthy label
    healthy_position = get_healthy_position(comparison_with_healthy,condition_matrix,disease_pos)
    # swap conditions if healthy is not the first condition
    if(healthy_position!=1){
      tmp = condition_2
      condition_2 = condition_1
      condition_1 = tmp
    }
    
    
    # degree of freedom should be <= 1
    df = count - common
    if(df<=1 | comparison_with_healthy){
      conditions_1 = c(conditions_1,condition_1)
      conditions_2 = c(conditions_2,condition_2)
      sample_conditions[[condition_1]] = rownames(conditions[conditions$condition==condition_1,,drop=F]) 
      sample_conditions[[condition_2]] = rownames(conditions[conditions$condition==condition_2,,drop=F]) 
    }
  }
  
  conditions_out = data.frame(condition_1=conditions_1, condition_2=conditions_2)

  return(list(conditions=conditions_out, samples=sample_conditions))
}


row_Means <- function(lvl,dds){
    sub_counts = counts(dds,normalized=TRUE)[,dds$condition == lvl, drop=F]
  if(length(colnames(sub_counts))>1){
    return(rowMeans(sub_counts))
  }else{
    return(sub_counts)
  }
  
}

#' DE Comparisons for all feasible combinations. Uses the concept of degree of freedoms and compares only
#' conditions in the sample group with the least differences. (n-1)
#' @param dds - dds object
#' @param conditions - conditions from sample group
#' @param datasets - datasets used for the analysis
#' @param outdir - output directory
#' @import DESeq2
#' @import stringr
#' @import SummarizedExperiment
#' @export
#'
de_comparisons <- function(dds, datasets, outdir, sample_group, formula, genome, overwrite=T, shrink=T) {
  # create an output directory
  conditions = colData(dds)
  valid_comparisons = get_comparisons(conditions,sample_group)
  condition_pairs = valid_comparisons$conditions
  condition_samples = valid_comparisons$samples
  
  if(is.null(condition_pairs$condition_1)){
    log_warning("Can't find valid comparisons")
  }else{
    for (i in 1:(length(condition_pairs$condition_1))) {
      tryCatch({
        condition_1 = as.character(condition_pairs$condition_1[i])
        condition_2 = as.character(condition_pairs$condition_2[i])
        
        contrast = c("condition", condition_1, condition_2)
        name = paste0(datasets,condition_1,condition_2)
        file_name = generate_filename('de_comparison',outdir, '.json', name=name ,hash=TRUE)
        
        # skip comparison if overwrite is false
        if(file.exists(file_name) & !overwrite){
          log_warning(paste0("File already exists: ", file_name))
          next
        }
        
        sample_con1 = rownames(conditions[as.character(conditions$condition)==condition_1,,drop=F])
        sample_con2 = rownames(conditions[as.character(conditions$condition)==condition_2,,drop=F])
          
        res = results(dds, contrast)
          
        if(shrink){
           res = lfcShrink(dds=dds, contrast=contrast, res=res, type = "normal")
        }
    
        # we add lfcSE, because if we use shrinkage with type normal, no lfcSE is included
        if(!("lfcSE" %in% names(res))){
          res["lfcSE"] = NA
        }
        
        # set flag for zero count genes
        res$zero_count = rowSums(counts(dds, normalized = F)[,c(sample_con1,sample_con2)])==0
    
        # get base mean for each condition
        baseMeanPerLvl <- sapply( levels(dds$condition),row_Means,dds=dds)
        baseMeanPerLvl_df = as.data.frame(baseMeanPerLvl)
        
        res$baseMean_A = baseMeanPerLvl_df[[condition_1]]
        res$baseMean_B = baseMeanPerLvl_df[[condition_2]]
        
        # sort by p-value
        resOrdered <- as.data.frame(res[order(res$pvalue), ])
        
        # TODO get genome id
        # genome = "gencode_28"
            
        # set meta da information: samples per condition, condition names and sample key values
        meta = generate_res_meta(datasets,
                                sample_con1,
                                sample_con2,
                                condition_1,
                                condition_2,
                                sample_group,
                                all.vars(formula),
                                genome)
        
        # output heatmap of top differential expressed genes
        dds_heatmap_differential_gene_counts(dds, res, 100, datasets, outdir, meta, name, c(sample_con1,sample_con2))
    
        write_json(list('metadata' = meta, 'data' = as.data.frame(resOrdered)),
                     file_name,
                     pretty = TRUE, na='string')
      
      }, error = function(e) {
        log_error(contrast)
        cat("ERROR: ", e$message, "\nin ")
        print(e$call)
        err <<- T
      })
    }
  }
}


#' Generates MA plot
#' @param dds - dds object
#' @export
#'
dds_MAPlot <- function(dds) {
  pdf("MaDiffExp.pdf")
  plotMA(dds)
  dev.off()
}


#' Generates heatmap for top {n} differential expressed genes and outputs to JSON file
#' based on normalized count data 
#' 
#' @param dds - dds object
#' @param dds_result - dds result under specific conditions
#' @param subset - number of genes required
#' @param datasets - datasets used for the analysis
#' @param outdir - output directory
#' @param meta - meta data information
#' @param name - file name 
#' @param sample_list - list of samples coresponding to conditions
#' @export
#'
dds_heatmap_differential_gene_counts <- function(dds, dds_result, subset, datasets, outdir, meta, name, sample_list) {
  select <- order(dds_result$pvalue)[1:subset]
  # create hash name
  name = substr(digest(name,algo="md5"),1,8)
  heatmapData = counts(dds, normalized = TRUE )[select,sample_list]
  
  # add meta information
  meta["row_number"] = length(rownames(heatmapData))
  meta["description"] = paste0("heatmap for top ",length(rownames(heatmapData))," differential expressed genes")
  meta["method"] = "normalized"
  heatmap(heatmapData, datasets, name, outdir, "normalized", description, meta)
}


#' Generates heatmap for top {n} expressed genes and outputs to JSON file
#' based on normalized count data 
#' 
#' @param dds - dds object
#' @param subset - number of genes required
#' @param datasets - datasets used for the analysis
#' @param outdir - output directory
#' @export
#'
dds_heatmap_gene_counts <- function(dds, subset, datasets, outdir) {
  mat = counts(dds, normalized = TRUE)
  heatmapData = mat[order(rowMeans(mat), decreasing = TRUE)[1:subset], ]
  description = paste0("heatmap for top ",length(rownames(heatmapData))," expressed genes")
  
  heatmap(heatmapData, datasets, "normalized", outdir, "normalized", description)
}

#' Generates heatmap for the top {n} genes from Regularized Log transform 
#' of the count data and outputs to JSON file
#' @param dds - dds object
#' @param subset - number of genes required
#' @param datasets - datasets used for the analysis
#' @param outdir - output directory
#' @param method - name of data transformation method
#' @export
#'
dds_heatmap_rld_vst <- function(dds, subset, datasets, outdir, method="raw") {
  mat = assay(dds)
  heatmapData = mat[order(rowMeans(mat), decreasing = TRUE)[1:subset], ]
  description = paste0("heatmap for top ",length(rownames(heatmapData))," expressed genes")

  heatmap(heatmapData, datasets, method, outdir, method, description)
}


#' Generates heatmap of the sample distances and outputs to JSON file
#' @param dds - dds object
#' @param datasets - datasets used for the analysis
#' @param outdir - output directory
#' @param method - name of data transformation method
#' @export
#'
dds_heatmap_distance <- function(dds, datasets, outdir, method="raw") {

  distsRL <- dist(t(assay(dds)))
  heatmapData <- as.matrix(distsRL)
  description = "heatmap of sample distances"

  heatmap(heatmapData, datasets, "sample_distance", outdir, method, description)
}


#' Generates PCA dataframe for the dds object, outputs PCA coordinates 
#' to JSON file.
#' @param dds - dds object
#' @param dataset - dataset names
#' @param outdir - output directory
#' @param sample_group - condition key names
#' @param transformation - name of transformation
#' @param plot - make scatter plot
#' @import DESeq2
#' @import SummarizedExperiment
#' @export
#'
dds_PCA <- function(dds, dataset_name, outdir, sample_group, transformation_name, plot=F) {
  conditions = SummarizedExperiment::colData(dds)
  mat = assay(dds)
  labels = conditions[colnames(mat),,drop = FALSE]$condition
  
  pca(mat, labels, dataset_name, outdir, sample_group, transformation_name , json=T, plot=plot)
  
}


#' Generates TSNE dataframe for the dds object for the first 3 dimensions
#' Explantion for the choice of parameters: https://cran.r-project.org/web/packages/Rtsne/Rtsne.pdf
#' @param dds - dds object
#' @param dataset - dataset names
#' @param outdir - output directory
#' @param sample_group - condition key names
#' @param perplexity - perplexity of TSNE analysis
#' @param transformation - name of transformation
#' @param plot - make scatter plot
#' @import DESeq2
#' @import SummarizedExperiment
#' @export
#'
dds_TSNE <- function(dds, dataset_name, outdir, sample_group, perplexity, transformation_name, plot=F) {
  
  conditions = SummarizedExperiment::colData(dds)
  mat = t(assay(dds))
  labels = conditions[rownames(mat),,drop=FALSE]$condition
  
  tsne(mat, labels, perplexity, dataset_name, outdir, sample_group, transformation_name , json=T, plot=plot)
}
