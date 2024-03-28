suppressPackageStartupMessages({
  require(optparse)
  require(data.table)
  require(stringr)
  require(igraph)
  require(networkD3)
  require(htmlwidgets)
})



# Define command-line options
option_list <- list(
make_option(c("-m",  "--matrix"),
            type = "character",
            help = "Path to the matrix from the dimensionality reduction step"), 
make_option(c("-a",  "--arrays"),
            type = "character",
            help = "Path to the matrix from the dimensionality reduction step"), 
make_option(c("-k",  "--k_neigh"),
            type = "integer",
            help = "Number of nearest neighboursh"),
make_option(c("-o",  "--output_path"),
            type = "character",
            help = "Output files path")
)


opt_parser <- OptionParser(option_list = option_list)

opt <- parse_args(opt_parser)

matrix <- fread(opt$matrix)


sat_name = names(matrix)[1] %>% str_remove("_R_") %>% str_split("_") %>% .[[1]] %>% .[2]

names <- colnames(matrix)
matrix <- cbind(names, matrix)
dt <- melt(matrix)


filename <- sat_name

dt[, variable  :=  str_remove(variable,paste0("_",sat_name))]
dt[, names :=  str_remove(names, paste0("_",sat_name))]

#locations of monomers
dt[, var_id  :=  str_remove(variable,"_R_")]
dt[, name_id :=  str_remove(names, "_R_")]

#locatiosn of arrays

arrays <- fread(opt$arrays)[V3==sat_name]


#overlap arrays with a synthetic data table contained

dt_synth <- dt[,.(var_id)] %>% unique(.)
dt_synth[,seqnames:=str_remove(var_id,"_\\d+")]
dt_synth[,start:={
  str_extract(var_id,"_\\d+") %>%
    str_remove(.,"_") %>% 
    as.double()
  }
  ]
dt_synth[,end:=start+1]

setkey(dt_synth,seqnames,start,end)
setkey(arrays,V1,V4,V5)

ar_id_var <- foverlaps(dt_synth,arrays)  %>% na.omit() %>% .[,.(var_id,var_array = paste0(seqnames,"_",V4))]


dt_synth <- dt[,.(name_id)] %>% unique(.)
dt_synth[,seqnames:=str_remove(name_id,"_\\d+")]
dt_synth[,start:={
  str_extract(name_id,"_\\d+") %>%
    str_remove(.,"_") %>% 
    as.double()
}
]
dt_synth[,end:=start+1]

setkey(dt_synth,seqnames,start,end)

name_id_var <- foverlaps(dt_synth,arrays) %>% na.omit() %>% .[,.(name_id,name_array = paste0(seqnames,"_",V4))]



ld_ar_dt <- dt[ar_id_var,on="var_id"][name_id_var,on="name_id"][var_array!=name_array]

#calculate mean distances between arrays
ld_ar_dt[, mval := mean(value, na.rm = TRUE), by = .(var_array, name_array)]

#find the closest array for each array
ld_ar_dt[, mmval := min(mval, na.rm = TRUE), by = .(name_array)]

tmp <- unique(ld_ar_dt[order(mval)][,.(var_array, name_array, mval)])[,head(.SD,  opt$k_neigh),  by=.(var_array)]

g <- graph_from_data_frame(tmp, directed = F)

p <- igraph_to_networkD3(g)

p$nodes$group <- str_remove(p$nodes$name, paste0("_",sat_name)) %>% str_remove(.,"_\\d+") 

graph <- forceNetwork(Links = p$links,  Nodes = p$nodes,  Source = "source", 
                     Target = "target",  NodeID = "name",
                     Group = "group",  Value = "value", 
                     zoom = TRUE,  linkDistance = 30, 
                     linkWidth = 1, 
                     arrows = FALSE, 
                     charge=-50, 
                     legend = TRUE,  opacity = 1)


saveWidget(graph, file = paste0(opt$output_path, filename, "_network.html"))
