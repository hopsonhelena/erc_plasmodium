topGeneClustersModuleUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    textInput(ns("gene_input"), "Enter Gene:", value = "PF3D7_0929400"),
    actionButton(ns("run_analysis"), "Get Top Genes"),
    hr(),
    DT::dataTableOutput(ns("result_table")),
    br(),
    downloadButton(ns("download_csv"), "Download as .CSV"),
    downloadButton(ns("download_excel"), "Download as .XLSX")
  )
}

topGeneClustersModule <- function(id, results_folder = "data/per_gene_rds") {
  moduleServer(id, function(input, output, session) {
    
    # Reactive: Load and process the RDS file
    results <- eventReactive(input$run_analysis, {
      req(input$gene_input)
      gene <- trimws(input$gene_input)
      file_path <- file.path(results_folder, paste0(gene, ".rds"))
      
      # 1. Check if the file exists
      if (!file.exists(file_path)) {
        showNotification(paste("No results found for gene:", gene),
                         type = "error", duration = 5)
        
        # Return empty dataframe with correct renamed column
        return(data.frame(
          `Input Gene` = character(0),
          `Output Gene` = character(0),
          `ERC Value` = numeric(0),
          `Percentile` = character(0),
          `Name` = character(0),
          `Description` = character(0),
          check.names = FALSE
        ))
      }
      
      # 2. Load the data
      df <- readRDS(file_path)
      
      # 3. Rename "Quantile" to "Percentile" if it exists
      if ("Quantile" %in% names(df)) {
        names(df)[names(df) == "Quantile"] <- "Percentile"
      }
      
      # 4. Column sanity check (using the new name)
      expected_cols <- c("Input Gene", "Output Gene", "ERC Value", "Percentile", "Name", "Description")
      missing <- setdiff(expected_cols, names(df))
      
      if (length(missing) > 0) {
        showNotification(
          paste("Warning: missing columns:", paste(missing, collapse = ", ")),
          type = "warning"
        )
      }
      
      return(df)
    })
    
    # Render results table
    output$result_table <- DT::renderDataTable({
      df <- results()
      req(df)
      if (nrow(df) == 0) return(NULL)
      
      DT::datatable(
        df,
        escape = FALSE,
        rownames = FALSE,
        options = list(
          pageLength = 25,
          lengthMenu = c(25, 50, 100),
          columnDefs = list(
            list(
              targets = 1,  # Index 1 is "Output Gene"
              render = JS(
                "function(data, type, row, meta) {",
                "  return '<a href=\"https://plasmodb.org/plasmo/app/record/gene/' + data + '\" target=\"_blank\">' + data + '</a>';",
                "}"
              )
            )
          )
        )
      )
    })
    
    # Download Handlers
    output$download_csv <- downloadHandler(
      filename = function() {
        paste0(input$gene_input, "-results-", Sys.Date(), ".csv")
      },
      content = function(file) {
        write.csv(results(), file, row.names = FALSE)
      }
    )
    
    output$download_excel <- downloadHandler(
      filename = function() {
        paste0(input$gene_input, "-results-", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        # openxlsx uses write.xlsx
        openxlsx::write.xlsx(results(), file)
      }
    )
  })
}
