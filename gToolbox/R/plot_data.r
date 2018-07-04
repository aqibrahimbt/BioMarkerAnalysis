
#' @import ggplot2
#' @export
scatter.plot <-function(data, x, y, label, file_name){
  g1 = ggplot(data, aes_string(x=x, y=y, colour=label)) + geom_point() + theme(legend.position="none")
  ggsave(paste0(file_name,".png"),plot=g1,width = 12, height = 12, device = "png", units = "cm")
}

#' Generates TSNE dataframe, outputs TSNE coordinates 
#' to JSON file or plot.
#' 
#' @param mat - matrix
#' @param labels - labels for pca points
#' @param perplexity - TSNE perpelexity 
#' @param dataset - dataset name
#' @param outdir - output directory
#' @param sample_group - condition key names
#' @param name - name prefix for output file
#' @param json - create json output
#' @param plot - make scatter plot
#' @import ggplot2
#' @import Rtsne
#' @export
#'
tsne <-function(mat, labels, perplexity=5, dataset_name="dataset", outdir="./", sample_group, name="default", json=T, plot=F){
  dim = 3
  n = nrow(mat) - 1
  perplexity = adjust_perplexity(n,perplexity)
  
  # calculate TSNE
  data = Rtsne(
    mat,
    check_duplicates = FALSE,
    pca = TRUE,
    perplexity = perplexity,
    theta = 0.5,
    dims = dim
  )
  data = as.data.frame(data$Y)
  colnames(data) <- paste("PC", 1:dim, sep = "")
  rownames(data) <- rownames(mat)
  
  data$labels = labels
  
  meta = list(dataset = dataset_name, samples=rownames(data), description = paste0("TSNE data analysis"), sample_group = sample_group, method=name)
  file_name = generate_filename(paste0('tsne_',name), outdir, '.json')
  
  if(json){
    write_json(list('metadata' = meta, 'data' = data), file_name,pretty = TRUE, na='string')
  }
  
  if(plot){
    scatter.plot(data,"PC1","PC2",data$labels,file_name)
  }
  
  return(data)
}


#' Generates PCA dataframe, outputs PCA coordinates 
#' to JSON file or plot.
#' 
#' @param mat - matrix
#' @param labels - labels for pca points
#' @param dataset - dataset name
#' @param outdir - output directory
#' @param sample_group - condition key names
#' @param name - name prefix for output file
#' @param json - create json output
#' @param plot - make scatter plot
#' @import ggplot2
#' @export
#'
pca <- function(mat, labels, dataset_name="dataset", outdir="./", sample_group, name="default", json=T, plot=F){
  
  # calculate pca
  pc <- prcomp(mat)
  data <- as.data.frame(pc$rotation)
  
  # set labels
  data$labels <- labels
  
  meta = list(dataset = dataset_name, samples=rownames(data), description = paste0("Principal components data analysis"), sample_group = sample_group, method=name)
  file_name = generate_filename(paste0('pca_',name), outdir, '.json')
  
  
  if(json){
    write_json(list('metadata' = meta, 'data' = data), file_name,pretty = TRUE, na='string')
  }
  
  if(plot){
    scatter.plot(data,"PC1","PC2",data$labels,file_name)
  }
  
  return(data)
}

#' Write heatmap to json and plot heatmap
#' @import ggplot2
#' @export
heatmap <- function(mat, datasets, name, outdir, method="",description="", meta=NA, plot=F) {

  if(length(meta)==1){
    if(is.na(meta)){
      meta = list(dataset = datasets, description = description, row_number=length(rownames(mat)), samples=colnames(mat), method=method)
    }
  }
  
  file_name = generate_filename(paste0("heatmap_",name), outdir, '.json')
  write_json(list('metadata' = meta, 'data' = as.data.frame(mat)),
             file_name,
             pretty = TRUE)
  
  if(plot){
    pheatmap(mat[select,], cluster_rows=F, show_rownames=T,cluster_cols=T)
  }
}