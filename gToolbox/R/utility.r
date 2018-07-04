

#' Check to see if packages are installed. Install them if they are not, then load them into the R session.
#' @param pkg - relevant packages for the project
#' @export
#'
check.packages <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE, repos = "http://cran.us.r-project.org")
  sapply(pkg, require, character.only = TRUE)
}

#'Check for bioconductor packages
#'@param pkg - relevant packages
#'@export
check.bioconductor.packages<- function(pkgs){
  source("https://bioconductor.org/biocLite.R")
  for(pkg in pkgs){
    if(!require(pkg,character.only=TRUE)){
      biocLite(pkg)
    }
  }
}


#' Generates a random string
#' @param x - length of variable
#' @export
randomizeString <- function(x) {
  a <- sample(1:1000, 1, replace = TRUE)
  return(a)
}

#' Generates log when running
#' @param msg - length of variable
#' @export
#'
log_running <- function(msg){
  return(paste0(format(Sys.time(), "%F %T"), " [RUNNING] ", msg))
}


#' Generates log when completed
#' @param msg - length of variable
#' @export
#'
log_done <- function(msg){
  return(paste0(format(Sys.time(), "%F %T"), " [DONE] ", msg))
}

#' Generates error log
#' @param msg - length of variable
#' @export
#'
log_error <- function(msg){
  return(paste0(format(Sys.time(), "%F %T"), " [ERROR] ", msg))
}

#' Generates filename based on the type of data generated
#' @param type - type pf data (eg. heatmap,comparisons etc.)
#' @param outdir - output directory of file
#' @param ext - file extension (e.g .json)
#' @import digest
#'
generate_filename <- function(type, outdir, ext, name="", hash=F) {
  file_name = type
  if(hash){
    name <- substr(digest(name,algo="md5"),1,8)
    file_name = paste(file_name,
                      name,
                      sep="_")
  }
  
  #file_name = paste(file_name,
  #                  format(Sys.time(), "%F"),
  #                  sep="_")
  
  file_name = paste0(file_name,ext)
  file_name = paste(outdir, file_name, sep = "/")
  return(file_name)
}


#' Generate meta information to be stored in JSON file
#' @param dtset - results
#' @param smp - condtions 1
#' @param cond_1 - condtions 2
#' @param cond_2 - dataset
#' @export
generate_res_meta <- function(dtset,sample_con1, sample_con2, cond_1, cond_2, sample_group, formula, genome) {
  meta <- list(dataset = dtset,
               samples_group_A=sample_con1, 
               samples_group_B=sample_con2,
               group_A = cond_1,
               group_B = cond_2,
               sample_group = sample_group,
               covariates = formula[!(formula %in% c("~","condition"))],
               genome=genome)
  return(meta)
}

#'Adjust perplexity
#'@param mat
#'@export
#'
adjust_perplexity <- function(n, perplexity){
  if(n < 3 * perplexity){
    perplexity = n/3
  }
  return(perplexity)
}


write_package_list <- function(file_name, dataset, version=0.01){
  si = sessionInfo()
  pkgs = si$otherPkgs
  
  pkg_list_df = data.frame(t(sapply(pkgs,function(pkg) return(c(pkg$Package,pkg$Version)))))
  colnames(pkg_list_df) = c("name","version")
  
  write_json(list(dataset=dataset, version=version,packages=pkg_list_df),
             file_name,
             pretty = TRUE, na='string')
}
