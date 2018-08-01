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
              help="don't overwrite existing files", action="store_false")
  ); 

opt_parser = OptionParser(usage = "usage: %prog [options]", option_list=option_list,  add_help_option = TRUE);
opt = parse_args(opt_parser,print_help_and_exit = TRUE);

BiocParallel::register(BiocParallel::SerialParam())

datasets = opt$dataset
output_dir = paste(opt$out,datasets,sep="/")
meta_data_path =  opt$meta
count_data_path =  opt$input
overwrite = opt$overwrite
dir.create(output_dir, showWarnings = FALSE)

datasets = "GSE65367"
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

classification <- function(dataset_data,samples_con1,samples_con2){
  
  norm_data = normalize_data(dataset_data, "deseq")
  
  sub_samples = c(samples_con1,samples_con2)
  countdata = as.data.frame(t(norm_data[,sub_samples]))
  conditions = dataset_data$conditions[rownames(countdata),]
  
  #Step 2: Recursive feature elimination
  
  profiles = recursive_feature_elimination(countdata, conditions, repeats = 2, number = 5)
  
  if(!length(profiles[[1]])){
    return(NULL)
  }
  
  #####################################################################################################################
  #Step 3: Aggregate data for results
  profiles[[1]]$resample
  aggregated_data = aggregate_profiles(profiles)
  
  profile_size = length(profiles)
  
  #####################################################################################################################
  #Step 4: Get dataset with Optimized variables and save logscale
  optimized_data <-as.data.frame(countdata[,profiles[[profile_size]]$optVariables])
  feature_logscale(optimized_data, aggregated_data, profiles[[profile_size]], output_dir)
  
  #####################################################################################################################
  #Step 5: Run Random Forest Classification
  mtry = c(1:5)
  rf.model <- train_model(optimized_data, conditions, "rf", mtry  ,repeats=2, number=5, balanced = F)
  
  rf.model
}

tryCatch({
    
  ## Read meta data flat file ##
  log_debug("get_metadata")
  meta_data = get_metadata(datasets)
  #meta_data = read.metadata.table(meta_data_path)
  
  if(is.null(meta_data)){
    log_error("Can't find valid meta data")
    stop()
  }
  
  formula = as.formula("~age+gender+condition")
  
  #write list of packages
  file_name = paste(output_dir,"packages.json",sep="/")
 
  write_package_list(file_name, datasets)
  log_debug("get.data")
  dataset_data = get.data(datasets,meta_data, count_data_path, file = "featureCounts.tsv", formula=formula)
  
  # Select condition of interrest (count data, labels)
  comparisons = get_comparisons(dataset_data$conditions, dataset_data$sample_group)
  
  for(i in 1:nrow(comparisons$conditions)){
    selected_comparision = comparisons$conditions[i,]
    samples_con1 = comparisons$samples[[as.character(selected_comparision[1,1])]]
    samples_con2 = comparisons$samples[[as.character(selected_comparision[1,2])]]
    classification(dataset_data,samples_con1,samples_con2)
  }
  


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


