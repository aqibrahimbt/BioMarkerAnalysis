#'@param package_name
#'@param script_source
#'
create_package <- function(package_name, script_source, packages, output_dir="./"){
  package_name = package_name
  script_source = script_source
  package_dir = paste(output_dir,package_name,sep="/")
  print(package_dir)

  if(file.exists(package_dir))
    unlink(package_dir, recursive = TRUE)
  
  devtools::create(package_dir)
  list.of.files = list.files (script_source, recursive = TRUE, full.name = TRUE)
  list.of.files <- list.of.files[ !grepl("pipelines", list.of.files) ]
  file.copy(list.of.files, paste0(package_dir, "/R"), overwrite = TRUE)
  setwd(package_dir)

  lapply(packages, function(x) devtools::use_package(x))
  devtools::document()
  devtools::build()

}
