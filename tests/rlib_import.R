vec <- c("BiocManager","ggplot2","data.table","dplyr","umap",
"stringr","FactoMineR","ape","optparse",
"htmlwidgets","igraph","networkD3","circlize",
"pheatmap","scales","Biostrings","ComplexHeatmap")


import_libraries <- function(library_names) {
  for (library_name in library_names) {
    tryCatch({
      library(library_name, character.only = TRUE, quietly = TRUE)
      cat(paste("Library '", library_name, "' imported successfully.\n"))
    }, error = function(e) {
      cat(paste("Error: Failed to import library '", library_name, "': ", conditionMessage(e), "\n", sep = ""))
    })
  }
}

import_libraries(vec)