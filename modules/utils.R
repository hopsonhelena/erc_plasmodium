# utils.R

reorder_mymatrix <- function(x) {
  dd <- dist(x, method = "euclidean")
  hc <- hclust(dd, method = "complete")
  return(x[hc$order, hc$order])  # Reorder matrix based on clustering
}

get_lower_tri <- function(x) {
  x[upper.tri(x)] <- NA
  return(x)
}

cleanList <- function(genes, mapping_table, erc_matrix) {
  mapping <- mapping_table %>%
    filter(gene %in% genes) %>%
    distinct(gene, cluster)
  
  clusters <- mapping$cluster
  names(clusters) <- mapping$gene  # Store gene names as names
  
  clusters <- clusters[clusters %in% colnames(erc_matrix)]
  return(clusters)
}

pairmat <- function(gene_list, erc_matrix, mapping_table, na.val = -2) {
  clusters <- cleanList(gene_list, mapping_table, erc_matrix)
  
  mat <- erc_matrix[clusters, clusters]
  rownames(mat) <- names(clusters)  
  colnames(mat) <- names(clusters) 
  
  diag(mat) <- NA
  mat[mat == na.val] <- NA
  return(mat)
}

pairlist <- function(gene_list, erc_matrix, mapping_table, na.val = -2) {
  cluster_list <- cleanList(gene_list, mapping_table, erc_matrix)
  mat <- erc_matrix[cluster_list, cluster_list]
  mat[mat == na.val] <- NA
  diag(mat) <- NA  
  return(mat[lower.tri(mat)])  
}

permTestMat <- function(gene_list, erc_matrix, mapping_table, perms = 1000) {
  # Convert gene list to cluster list
  cluster_list <- cleanList(gene_list, mapping_table, erc_matrix)
  
  # Compute observed mean ERC
  obs_values <- pairlist(gene_list, erc_matrix, mapping_table)
  obs_mean <- mean(obs_values, na.rm = TRUE)
  
  # If no valid ERC values, return NaN
  if (is.nan(obs_mean)) {
    return(list(obs = NaN, p = NaN, null = NaN))
  }
  
  # Get all clusters in the ERC matrix
  all_clusters <- colnames(erc_matrix)
  
  # Compute null distribution
  null_dist <- numeric(perms)
  for (i in 1:perms) {
    # Sample random clusters
    random_clusters <- sample(all_clusters, length(cluster_list))
    
    # Extract ERC values for the random clusters
    random_erc_subset <- erc_matrix[random_clusters, random_clusters]
    random_erc_subset[random_erc_subset == -2] <- NA
    diag(random_erc_subset) <- NA
    random_lower_tri <- random_erc_subset[lower.tri(random_erc_subset)]
    
    # Calculate the mean ERC for the random sample
    null_dist[i] <- mean(random_lower_tri, na.rm = TRUE)
  }
  
  # Compute empirical p-value
  p_value <- sum(obs_mean <= null_dist, na.rm = TRUE) / perms
  
  return(list(obs = obs_mean, p = p_value, null = null_dist))
}