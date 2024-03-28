

install_packages_if_not_installed <- function(package_names) {
  for (package in package_names) {
    if (requireNamespace(package, quietly = TRUE)) {
      cat(paste("Package '", package, "' is already installed.\n"))
    } else {
      cat(paste("Package '", package, "' is not installed.
       Install package? (yes|no).\n"))
      user_input <- scan(what = character(), n = 1, quiet = TRUE)
      
      if (user_input == "yes") {
        install.packages(package)
        library(package, character.only = TRUE)
      } else {
        cat(paste("Skipping installation of '", package, "'.\n"))
      }
    }
  }
}




vec <- c("BiocManager","ggplot2","data.table","dplyr","umap",
"stringr","FactoMineR","ape","optparse",
"htmlwidgets","igraph","networkD3","circlize",
"pheatmap","scales")

install_packages_if_not_installed(vec)

print("Installing Biocmanager only packages:")
BiocManager::install("Biostrings")
BiocManager::install("ComplexHeatmap")

