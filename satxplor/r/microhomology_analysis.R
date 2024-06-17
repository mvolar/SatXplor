
suppressPackageStartupMessages({
  require(Biostrings)
  require(ggseqlogo)
  require(optparse)
  require(stringr)
  require(dplyr)
  require(ggplot2)
  
  
})

option_list <- list(
  make_option(c("-a", "--flank_alignment"), type = "character", help = "Sequence alignment file"),
  make_option(c("-o", "--output_path"), type = "character", help = "Output files path")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


msa  <- readDNAMultipleAlignment(opt$flank_alignment,format="fasta")

filename = str_remove(rownames(msa),"_R_")[1] %>% str_split(.,"_") %>% unlist() %>% .[3]


matrix <- consensusMatrix(msa,as.prob = T)[1:4,]


p <- ggplot() + geom_logo( matrix,
                      method = "probability" ) + theme_logo() +
  theme(legend.position = "none")

ggsave(plot=p,filename=paste0(opt$output_path,filename,"_seqlogo.png"),width=15,height=3)