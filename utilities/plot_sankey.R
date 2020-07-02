# install.packages("networkD3")
library(networkD3)
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly=TRUE)
table_file_in <- args[1]
sankey_file_out <- args[2]


links_sorted <- read.table(table_file_in, header = T, sep = "\t")

get_node_LUT <- function(my_links_tbl){

  new_source <- sapply(my_links_tbl$source, FUN = function(x){
    if(startsWith(x, 'UK')){
      x <- paste0(x, "_previous", collapse = '')
    }
    return(x)
  })
  names(new_source) <- NULL

  new_target <- sapply(my_links_tbl$target, FUN = function(x){
    if(startsWith(x, 'UK')){
      x <- paste0(x, "_latest", collapse = '')
    }
    return(x)
  })
  names(new_target) <- NULL

  x <- data.frame(cbind('source' = new_source, 'target' = new_target,'value' = my_links_tbl$value))

  sources <- unique(unlist(x$source))
  targets <- unique(unlist(x$target))

  new_slice <- sources == 'new'
  if(any(new_slice)){
    nodes <- c('new', sources[!new_slice])
  } else {
    nodes <- sources
  }

  for(i in seq_along(targets))(
    if(!targets[i] %in% nodes){
      nodes <- c(nodes, targets[i])
    }
  )

  idx <- seq(0, length(nodes) - 1, 1)
  df <- idx
  names(df) <- nodes

  return(df)
}

get_links <- function(my_links_tbl, LUT){
  new_source <- sapply(my_links_tbl$source, FUN = function(x){
    if(startsWith(x, 'UK')){
      x <- paste0(x, "_previous", collapse = '')
    }
    return(LUT[x])
  })
  names(new_source) <- NULL

  new_target <- sapply(my_links_tbl$target, FUN = function(x){
    if(startsWith(x, 'UK')){
      x <- paste0(x, "_latest", collapse = '')
    }
    return(LUT[x])
  })
  names(new_target) <- NULL

  flat <- data.frame(cbind('source' = new_source, 'target' = new_target,'value' = my_links_tbl$value))
  flat$source <- as.integer(flat$source)
  flat$target <- as.integer(flat$target)
  return(flat)
}

# print(links_sorted)
node_LUT <- get_node_LUT(links_sorted)
# print(node_LUT)
# print(names(node_LUT))
links <- get_links(links_sorted, node_LUT)
# print(links)
links$value <- unlist(links$value) / max(unlist(links$value))

nodes <- as.data.frame(cbind('node' = seq(0, length(node_LUT) - 1, 1),
                             'name' = sapply(names(node_LUT), FUN = function(x){
                               strsplit(x, '_')[[1]][1]
                             })))
saveNetwork(
  sankeyNetwork(Links = links, Nodes = nodes,
                Source = "source",
                Target = "target",
                Value = "value",
                NodeID = "name",
                fontSize = 30,
                nodeWidth = 50,
                height = 2000,
                width = 1000,
                nodePadding = 20,
                sinksRight = F,
                iterations = 10,
                margin= list('top' = 0, 'bottom' = 0, 'left' = 0, 'right' = 0)),
  file = sankey_file_out,
  selfcontained = TRUE)
