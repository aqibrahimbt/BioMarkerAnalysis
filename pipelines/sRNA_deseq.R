#Initialize libraries and global variables

suppressWarnings({
  suppressMessages({
    library("gToolbox")
    packages <-c("ggplot2","jsonlite","stringr","Rtsne","digest","optparse")
    check.packages(packages)
    bio_packages <-c("DESeq2","genefilter")
    check.bioconductor.packages(bio_packages)
    Sys.setenv(TZ="Europe/Berlin")
  })
})


#options(error=recover, show.error.locations=TRUE, warn=2)

# read input parameter
option_list = list(
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output directory", metavar="character"),
  make_option(c("-m", "--meta"), type="character", default= NULL, 
              help="meta data table", metavar="character"),
  make_option(c("-i", "--input"), type="character", default= NULL, 
              help="working directory of input files", metavar="character"),
  make_option(c("-d", "--dataset"), type="character", default= NULL, 
              help="dataset id", metavar="character"),
  make_option(c("-w", "--overwrite"), type="logical", default= TRUE, 
              help="don't overwrite existing files")
  ); 

opt_parser = OptionParser(usage = "usage: %prog [options]", option_list=option_list,  add_help_option = TRUE);
opt = parse_args(opt_parser,print_help_and_exit = TRUE);

datasets = opt$dataset
output_dir = paste(opt$out,datasets,sep="/")
meta_data_path =  opt$meta
count_data_path =  opt$input
overwrite = opt$overwrite
dir.create(output_dir, showWarnings = FALSE)

datasets = "GSE69089"
output_dir = "./output/"
meta_data_path =  "./data/meta_data/meta_data_20180209.tsv"
count_data_path =  "./data/test_data/"
dir.create(output_dir, showWarnings = FALSE)

# redirect R output to file
log_file = paste(output_dir,paste0(datasets,".log"),sep="/")
log_level <<- "INFO"
log_running(datasets)
fp = file(log_file, open = "w")
sink(file = fp)
sink(file = fp, type = "message")
log_running(datasets)

err = F

tryCatch({
    
  ## Read meta data flat file ##
  log_debug("get_metadata")
  meta_data = get_metadata(datasets)
  #meta_data = read.metadata.table(meta_data_path)
  formula = as.formula("~age+gender+condition")
  
  #write list of packages
  file_name = paste(output_dir,"packages.json",sep="/")
 
  write_package_list(file_name, datasets)
  log_debug("get.data")
  dataset_data = get.data(datasets,meta_data, count_data_path, file = "featureCounts.tsv", formula=formula)
  
  sample_group = dataset_data$sample_group
  
  ## Calculate differential feature expression
  log_debug("dds")
  dds <- de_result(dataset_data)
  de_comparisons(dds, datasets, output_dir,sample_group,dataset_data$formula, overwrite=overwrite, shrink=T)
  
  
  ## Count data transformations 
  vst <- varianceStabilizingTransformation(dds)
  rld <- rlog(dds, blind=FALSE)

  dds_heatmap_gene_counts(dds, 100, datasets,output_dir)
  
  dds_heatmap_rld_vst(vst, 100, datasets,output_dir,"vst")
  dds_heatmap_rld_vst(rld, 100, datasets,output_dir,"rld")
  
  dds_heatmap_distance(vst, datasets,output_dir,"vst")
  
  pca = dds_PCA(vst, datasets, output_dir,sample_group,"vst",T)
  tsne = dds_TSNE(vst, datasets, output_dir,sample_group,5, "vst",T)
  
  pca = dds_PCA(rld, datasets, output_dir,sample_group,"rld",T)
  tsne = dds_TSNE(rld, datasets, output_dir,sample_group,5, "rld",T)

}, error = function(e) {
  e <<- e
  cat("ERROR: ", e$message, "\nin ")
  print(e$call)
  err <<- T
})

if(err){
  log_error(datasets)
  sink()
  sink(type = "message")
  log_error(datasets)
}else{
  log_done(datasets)
  sink()
  sink(type = "message")
  log_done(datasets)
}


