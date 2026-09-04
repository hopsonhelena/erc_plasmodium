# result_matrix_module.R
resultMatrixModuleUI <- function(id) {
  ns <- NS(id)
  tagList(
    htmlOutput(ns("missing_genes_warning"), class = "erc-warning"),
    htmlOutput(ns("mean_erc"), class = "erc-text"),
    htmlOutput(ns("p_value"), class = "erc-text"),
    DTOutput(ns("result_matrix")),
    helpText("Note: 'NA' indicates a gene pair without an ERC value."),
    downloadButton(ns("download_matrix"), "Download Matrix as CSV")
  )
}

resultMatrixModule <- function(id, gene_list, gene_to_cluster_map, erc_matrix) {
  moduleServer(id, function(input, output, session) {
    # Generate the result matrix
    result_matrix <- reactive({
      req(gene_list())
      mat <- pairmat(gene_list(), erc_matrix, gene_to_cluster_map)
      
      # Round actual values to 3 decimals
      mat <- round(mat, 3)
      
      # Store real NAs before masking
      real_na_mask <- is.na(mat)
      
      # Mask upper triangle
      mat[upper.tri(mat)] <- Inf  # Temporary marker
      
      # Replace values for display:
      display_mat <- mat
      display_mat[real_na_mask] <- "NA"       # Show real missing values
      display_mat[display_mat == Inf] <- ""   # Hide upper triangle
      
      return(display_mat)
    })
    
    # Calculate mean ERC and p-value
    test_results <- reactive({
      req(gene_list())
      permTestMat(gene_list(), erc_matrix, gene_to_cluster_map)
    })
    
    # Display mean ERC
    output$mean_erc <- renderText({
      paste("Mean ERC =", round(test_results()$obs, 3))
    })
    
    # Display p-value
    output$p_value <- renderText({
      paste("Empirical p-value =", round(test_results()$p, 3))
    })
    
    # Display the result matrix (lower triangular only, with "NA" for missing values)
    output$result_matrix <- renderDT({
      mat <- result_matrix()
      display_mat <- mat
      display_mat[is.na(display_mat)] <- "NA"  # Show "NA" as text for missing values
      datatable(display_mat, options = list(pageLength = 10))
    })
    
    # Download matrix as CSV
    output$download_matrix <- downloadHandler(
      filename = function() {
        "result_matrix.csv"
      },
      content = function(file) {
        write.csv(result_matrix(), file)
      }
    )
    
    missing_genes <- reactive({
      req(gene_list())
      input_genes <- unique(trimws(gene_list()))
      all_erc_genes <- unique(gene_to_cluster_map$gene)
      setdiff(input_genes, all_erc_genes)
    })
    
    # Render warning message if any genes are missing
    output$missing_genes_warning <- renderUI({
      missing <- missing_genes()
      if (length(missing) > 0) {
        HTML(paste0(
          "<div style='color: red; font-weight: bold;'>",
          "The following gene(s) are not found in the ERC dataset: ",
          paste(missing, collapse = ", "),
          ". This may result in missing values or errors.",
          "</div>"
        ))
      } else {
        NULL
      }
    })
    
    
    
    return(list(
      test_results = test_results
    ))
  }) # Close moduleServer
} # Close resultMatrixModule