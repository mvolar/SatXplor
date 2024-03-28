options(conflicts.policy = list(warn = FALSE))
library(umap,quietly = TRUE)
library(optparse,quietly = TRUE)
library(data.table,quietly = TRUE)
library(FactoMineR,quietly = TRUE)
library(ape,quietly = TRUE)
library(stringr,quietly = TRUE)
library(Biostrings,quietly = TRUE)
library(ggplot2,quietly = TRUE)


# Define command-line options
option_list <- list(
  make_option(c("-a", "--sequence_alignment"),
  type = "character", help = "Sequence alignment file"),
  make_option(c("-d", "--dim_red"),
  type = "character", help = "Dimensionality reduction option ('umap' or 'pca' or 'both')"),
  make_option(c("-o", "--output_path"),
  type = "character", help = "Output files path"),
  make_option(c("-m", "--save_matrix"),default = TRUE,
  type = "logical", help = "Path to save matrices, if empty matrices wont be saved.")
  
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

print(opt)

# Check for required arguments
if (is.null(opt$sequence_alignment) || is.null(opt$dim_red) || is.null(opt$output_path)) {
  cat("Usage: Rscript dimensionality_reduction_optparse.R --sequence_alignment <file> --dimensionality_reduction <umap_or_pca> --output_file <file>\n")
  quit(save = "no", status = 1)
}


msa <- readDNAMultipleAlignment(opt$sequence_alignment,
                                format="fasta")



filename = str_remove(rownames(msa),"_R_")[1] %>% str_split(.,"_") %>% unlist() %>% .[2]




x <- dist.dna(as.DNAbin(msa),model="F81",as.matrix=TRUE,pairwise.deletion=TRUE)

cat("Finished distance calculating: \n",
    opt$sequence_alignment)


matrix <- as.data.table(x,row.names="V1")

if (opt$save_matrix==TRUE)
{
  fwrite(matrix,paste0(opt$output_path,filename,"_aligned_matrix.csv.gz"))
}


if (opt$dim_red=="pca"){
  

  
  pca_res <- PCA(matrix,graph=FALSE)
  cat("Finished PCA calculating: \n",
      opt$sequence_alignment)
  
  dt <- pca_res$var$coord %>% as.data.table(keep.rownames = TRUE) 
  
  
  eigenvalues <- pca_res$eig[, 2]
  
  # Scree plot using ggplot2
  scree_data <- data.table(
    Principal_Component = seq_along(eigenvalues),
    Variance_Explained = eigenvalues
  )
  
  
  p1 <- ggplot(scree_data[1:10], aes(x = as.factor(Principal_Component), y = Variance_Explained)) +
    geom_bar(stat = "identity", fill = "skyblue", width = 0.7,color="black") +
    geom_text(aes(label = sprintf("%.2f", Variance_Explained)), vjust = -0.5) +
    labs(x = "Principal Component", y = "Variance Explained", title = "Scree Plot") +
    theme_minimal()
  
  
  
  ggsave(filename=paste0(opt$output_path,"/",filename,"_monomers_PCA_scree_plot.png"),
         plot=p1)
  dt[,ar_id:=str_remove(rownames(msa),"_R_")]
  
  #default behaviour is to color based on seqnames values
  dt[,origin:=str_split(dt$ar_id,"_",simplify = TRUE)[,1]]
  
  p2 <- ggplot(dt) + geom_point(aes(x=Dim.1,y=Dim.2,color=origin),alpha=0.4) + theme_bw() + theme(legend.position = "none")
  ggsave(plot=p2, filename = paste0(opt$output_path,"/",filename,"_monomers_PCA_plot.png"))
}


if (opt$dim_red=="umap"){
  
  pca_res <- umap(setnafill(matrix,fill=0))
  
  cat("Finished UMAP calculating: \n",
      opt$sequence_alignment)
  
  pca_r_dt <- pca_res$layout %>% as.data.table()
  
  pca_r_dt[,ar_id:=str_remove(rownames(msa),"_R_")]
  
  #default behaviour is to color based on seqnames values
  pca_r_dt[,origin:=str_split(dt$ar_id,"_",simplify = TRUE)[,1]]
  
  p1 <- ggplot(pca_r_dt) + geom_point(aes(x=V1,y=V2,color=origin),alpha=0.4) + theme_bw() + theme(legend.position = "none")
  ggsave(plot=p1, filename = paste0(opt$output_path,"/",filename,"_monomers_UMAP_plot.png"))
  
}

if (opt$dim_red=="both"){
  
  pca_res <- PCA(matrix,graph=FALSE)
  cat("Finished PCA calculating: \n",
      opt$sequence_alignment)
  
  dt <- pca_res$var$coord %>% as.data.table(keep.rownames = TRUE) 

  eigenvalues <- pca_res$eig[, 2]
  
  # Scree plot using ggplot2
  scree_data <- data.table(
    Principal_Component = seq_along(eigenvalues),
    Variance_Explained = eigenvalues
  )
  
  
  p1 <- ggplot(scree_data[1:10], aes(x = as.factor(Principal_Component), y = Variance_Explained)) +
    geom_bar(stat = "identity", fill = "skyblue", width = 0.7,color="black") +
    geom_text(aes(label = sprintf("%.2f", Variance_Explained)), vjust = -0.5) +
    labs(x = "Principal Component", y = "Variance Explained", title = "Scree Plot") +
    theme_minimal()
  
  ggsave(filename=paste0(opt$output_path,"/",filename,"_monomers_PCA_scree_plot.png"),
         plot=p1)
  dt[,ar_id:=str_remove(rownames(msa),"_R_")]
  
  #default behaviour is to color based on seqnames values
  dt[,origin:=str_split(dt$ar_id,"_",simplify = TRUE)[,1]]
  
  p2 <- ggplot(dt) + geom_point(aes(x=Dim.1,y=Dim.2,color=origin),alpha=0.4) + theme_bw() + theme(legend.position = "none")
  ggsave(plot=p2, filename = paste0(opt$output_path,"/",filename,"_monomers_PCA_plot.png"))
  
  
  pca_res <- umap(setnafill(matrix,fill=0))
  cat("Finished UMAP calculating: \n",
      opt$sequence_alignment)
  pca_r_dt <- pca_res$layout %>% as.data.table()
  
  pca_r_dt[,ar_id:=str_remove(rownames(msa),"_R_")]
  
  
  #default behaviour is to color based on seqnames values
  pca_r_dt[,origin:=str_split(dt$ar_id,"_",simplify = TRUE)[,1]]
  
  
  p1 <- ggplot(pca_r_dt) + geom_point(aes(x=V1,y=V2,color=origin),alpha=0.4) + theme_bw() + theme(legend.position = "none")
  
  ggsave(plot=p1, filename = paste0(opt$output_path,"/",filename,"_monomers_UMAP_plot.png"))
  
}





