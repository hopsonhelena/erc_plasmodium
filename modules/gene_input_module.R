# gene_input_module.R

geneInputModuleUI <- function(id) {
  ns <- NS(id)
  tagList(
    fileInput(ns("file_input"), "Upload Gene List (TXT or XLSX)", accept = c(".txt", ".xlsx")),
    textAreaInput(ns("gene_list_input"), "Or Enter Gene List (space or comma-separated)", rows = 3),
    actionButton(ns("submit_genes"), "Submit Genes")
  )
}

geneInputModule <- function(id) {
  moduleServer(id, function(input, output, session) {
    genes <- reactiveVal(NULL)  # Store the gene list
    observeEvent(input$submit_genes, {
      gene_list <- NULL  # Initialize gene_list
      
      if (!is.null(input$file_input)) {
        file_path <- input$file_input$datapath
        ext <- tools::file_ext(file_path)
        
        if (ext == "txt") {
          gene_list <- readLines(file_path)
        } else if (ext == "xlsx") {
          gene_list <- readxl::read_excel(file_path)[[1]]
        } else {
          showNotification("Invalid file format. Please upload a .txt or .xlsx file.", type = "error")
          return(NULL)
        }
      } else if (!is.null(input$gene_list_input) && input$gene_list_input != "") {
        gene_list <- unlist(strsplit(input$gene_list_input, "[,\\s]+", perl = TRUE)) %>% trimws()
      } else {
        showNotification("Please enter genes manually or upload a file.", type = "error")
        return(NULL)
      }
      
      # Validate the gene list
      if (length(gene_list) == 0) {
        showNotification("The gene list is empty. Please check your input.", type = "error")
        return(NULL)
      }
      
      # Update the reactive value
      genes(trimws(gene_list))
    })
    
    return(genes)  # Return the reactive expression
  })
}