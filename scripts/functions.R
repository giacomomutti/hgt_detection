# rescale weight to be from 0 to 1 but if 0 add 0.0001!!
range01 <- function(x) { 
  new <- (x-min(x))/(max(x)-min(x))
  new <- ifelse(new==0, 0.0001, new)
}


filter_df <- function(df, min_pident=80, min_cov=80, minevalue=1e-10, byevalue=FALSE){

  if (byevalue) {
    df <- filter(df, evalue<=minevalue, pident >= min_pident, max_cov >= min_cov, evalue<=minevalue) %>% 
      mutate(weight=evalue)
    
  } else {
    df <- filter(df, evalue<=minevalue, pident >= min_pident, max_cov >= min_cov) %>%
      mutate(weight=bit_self)
  }
  
}


get_network <- function(filtered_df, simplify=FALSE) {

  g <- graph_from_data_frame(filtered_df, directed = F)

  E(g)$weight <- filtered_df$weight
  E(g)$bit_self <- filtered_df$bit_self
  
  # color scale from wesanderson package
  c_scale <- colorRamp(as.character(wesanderson::wes_palette("Zissou1")))
  
  if (simplify) {
    # simplify the graph and take mean weight if multimap
    g <- simplify(g, edge.attr.comb = list(weight="mean", "ignore"))
  } 
  #Applying the color scale to edge weights.
  #rgb method is to convert colors to a character vector.
  E(g)$color = apply(c_scale(range01(E(g)$bit_self)), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
  
  return(g)
}

# alternatively you can plot an heatmap
# sometimes two hsp are for the same pair of query and target!!
# ident_mat <- df %>% group_by(qseqid, sseqid) %>% 
#   summarise(pident=mean(pident)) %>% 
#   ungroup() %>% 
#   pivot_wider(names_from = sseqid, values_from = pident, values_fill = 0) %>% 
#   column_to_rownames("qseqid") %>% 
#   as.matrix()  

# heatmap(ident_mat)
