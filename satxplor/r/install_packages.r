

install_packages_if_not_installed <- function(package_names) {
  for (package in package_names) {
    if (!requireNamespace(package, quietly = TRUE)) {
      cat(paste("Installing package '", package, "'.\n"))
      install.packages(package, dependencies = TRUE)
      library(package, character.only = TRUE)
    } else {
      cat(paste("Package '", package, "' is already installed.\n"))
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

