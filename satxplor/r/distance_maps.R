suppressPackageStartupMessages({
  require(pheatmap)
  require("ComplexHeatmap")
  require(circlize)
  require(ape)
  require(data.table)
  require(Biostrings)
  require(scales)
  require(optparse)
  require(dplyr)
  require(stringr)
})

option_list <- list(
  make_option(c("-a", "--flank_alignment"), type = "character", help = "Sequence alignment file"),
  make_option(c("-o", "--output_path"), type = "character", help = "Output files path")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)



msa <- readDNAMultipleAlignment(opt$flank_alignment,format = "fasta")

filename = str_remove(rownames(msa),"_R_")[1] %>% str_split(.,"_") %>% unlist() %>% .[3]

x <- dist.dna(as.DNAbin(msa),model="F81",as.matrix=TRUE,pairwise.deletion=TRUE) 

x<-1-((max(x,na.rm=TRUE) - x)/max(x,na.rm=TRUE))

x<-rescale(x, to = c(0, 1))

dt <- setnafill(as.data.table(x),fill=min(x,na.rm=TRUE))

col_fun = colorRamp2(c(0,0.5,1), c("#AD3212","#F4ED7E","#1A5276"))

h1=Heatmap(as.matrix(1-dt),show_column_names = FALSE,col = col_fun,name=" ",
             heatmap_legend_param = list(
               title = "similarity", at = c(0,0.5,1)
             ))

png(paste0(opt$output_path,filename,"_flank_distances.png"))
draw(h1)
dev.off()



