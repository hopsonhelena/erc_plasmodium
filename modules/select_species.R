species_choices <- c(
  "P. falciparum 3D7" = "falc",
  "P. berghei ANKA" = "berg",
  "P. vivax P01" = "viv",
  "P. reichenowi CDC" = "rei",
  "P. coatneyi Hackeri" = "coat",
  "P. malariae UG01" = "mal",
  "P. ovale curtisi GH01" = "oval",
  "P. knowlesi H" = "know",
  "P. cynomolgi M" = "cyno",
  "P. inui SanAntonio1" = "inu",
  "P. fragile Nilgiri" = "frag",
  "P. chabaudi chabaudi" = "chab",
  "P. yoelii yoelii 17X" = "yoel",
  "P. vinckei vinckei vinckei" = "vin",
  "P. gaboni G01" = "gab",
  "P. adleri G01" = "adl",
  "P. blacklocki G01" = "blac",
  "P. billcollinsi G01" = "bill",
  "P. praefalciparum G01" = "prae",
  "P. vivax-like Pvl01" = "vivl",
  "P. relictum SGS1" = "rel",
  "P. gallinaceum 8A" = "gal"
)





geneToSpeciesModuleUI <- function(id) {
  ns <- NS(id)
  
  
  tagList(
    textInput(ns("gene_input"), "Enter Gene(s) (comma or space-separated):", value = "PF3D7_0929400"),
    
    fileInput(ns("gene_file"), "Upload Excel File with Gene List", 
              accept = c(".xlsx", ".xls")),
    
    checkboxGroupInput(ns("species_select"), "Select Species:",
                       choices = species_choices,
                       selected = c("falc", "berg")),
    
    actionLink(ns("select_all"), "Select All"),
    actionLink(ns("deselect_all"), "Deselect All"),
    
    div(
      style = "margin-top: 10px;",
      actionButton(ns("clear_selection"), "Clear Selection")
    ),
    
    div(
      style = "margin-top: 10px;",
      actionButton(ns("find_genes"), "Find Orthologs in Selected Species")
    ),
    
    div(
      style = "margin-top: 5px;",
      conditionalPanel(
        condition = paste0("output['", ns("result_table"), "_shown'] === true"),
        DT::dataTableOutput(ns("result_table"))
      )
    ),
    
    div(
      style = "margin-top: 10px;",
      downloadButton(ns("download_csv"), "Download as .CSV")
    ),
    
    div(
      style = "margin-top: 10px;",
      downloadButton(ns("download_excel"), "Download as .XLSX")
    ),
    div(
      style = "margin-top: 15px; padding: 10px; color: red; font-weight: bold; border: 1px solid red;",
      textOutput(ns("missing_genes_box"))
    )
  )
}


geneToSpeciesModule <- function(id, gene_to_cluster_long) {
  moduleServer(id, function(input, output, session) {
    # Define species names and prefix-to-family mapping
    missing_genes_rv <- reactiveVal(character(0))
    priority_order <- c("falc", "berg", "viv", "rei", "coat", "mal", "oval", "know",
                        "cyno", "inu", "frag", "chab", "yoel", "vin", "gab", "adl",
                        "blac", "bill", "prae", "vivl", "rel", "gal")
    
    prefix_to_family <- c(
      PBANKA = "berg",
      PF3D7 = "falc",
      PCOAH = "coat",
      PcyM = "cyno",
      AK88 = "frag",
      PGAL8A = "gal",
      C922 = "inu",
      PKNH = "know",
      PmUG01 = "mal",
      PRELSG = "rel",
      PVP01 = "viv",
      PVL = "vivl",
      PY17X = "yoel",
      PADL01 = "adl",
      PBILCG01 = "bill",
      PBLACG01 = "blac",
      PCHAS = "chab",
      PGABG01 = "gab",
      PocGH01 = "oval",
      PPRFG01 = "prae",
      PRCDC = "rei",
      YYE = "vin"
    )
    
    # Short code to full species name mapping
    species_full_names <- c(
      falc = "P. falciparum 3D7",
      berg = "P. berghei ANKA",
      viv = "P. vivax P01",
      rei = "P. reichenowi CDC",
      coat = "P. coatneyi Hackeri",
      mal = "P. malariae UG01",
      oval = "P. ovale curtisi GH01",
      know = "P. knowlesi H",
      cyno = "P. cynomolgi M",
      inu = "P. inui SanAntonio1",
      frag = "P. fragile Nilgiri",
      chab = "P. chabaudi chabaudi",
      yoel = "P. yoelii yoelii 17X",
      vin = "P. vinckei vinckei vinckei",
      gab = "P. gaboni G01",
      adl = "P. adleri G01",
      blac = "P. blacklocki G01",
      bill = "P. billcollinsi G01",
      prae = "P. praefalciparum G01",
      vivl = "P. vivax-like Pvl01",
      rel = "P. relictum SGS1",
      gal = "P. gallinaceum 8A"
    )
    
    # Clear species selection
    observeEvent(input$clear_selection, {
      updateSelectInput(session, "species_select", selected = character(0))
    })
    
    
    observeEvent(input$select_all, {
      updateCheckboxGroupInput(session, "species_select",
                               selected = unname(species_choices))
    })
    
    # Deselect All
    observeEvent(input$deselect_all, {
      updateCheckboxGroupInput(session, "species_select",
                               selected = character(0))
    })
    
    
    
    # Reactive expression to read and clean uploaded Excel file
    uploaded_genes <- reactive({
      req(input$gene_file)
      file <- input$gene_file
      genes_df <- read_excel(file$datapath, col_names = TRUE)  # Read with headers
      col_name <- names(genes_df)[1]  # Assume the first column contains gene names
      genes <- as.character(genes_df[[col_name]])  # Extract and convert to character
      genes <- trimws(genes)  # Remove extra spaces
      genes <- genes[!is.na(genes) & genes != ""]  # Remove empty values
      return(genes)
    })
    
    ########## Reactive expression to process gene-to-species mapping


    results <- eventReactive(input$find_genes, {
      req(input$species_select)
      
      # Get genes either from file or text input, support comma/space-separated
      genes <- if (!is.null(input$gene_file)) {
        uploaded_genes()
      } else {
        unlist(strsplit(gsub("[,\\n]+", " ", input$gene_input), "\\s+"))
      }
      genes <- trimws(genes)
      genes <- genes[genes != ""]
      
      # Keep track of missing genes
      missing_genes <- c()
      
      result_list <- lapply(genes, function(gene) {
        gene_entry <- gene_to_cluster_long[gene_to_cluster_long$gene == gene, ]
        
        if (nrow(gene_entry) == 0) {
          missing_genes <<- c(missing_genes, gene)
          return(NULL)
        }
        
        cluster_id <- unique(gene_entry$cluster)
        
        cluster_genes <- gene_to_cluster_long[gene_to_cluster_long$cluster == cluster_id, ]
        
        result_row <- data.frame(Gene = gene, stringsAsFactors = FALSE)
        
        for (species in input$species_select) {
          orthologs <- cluster_genes$gene[cluster_genes$family == species]
          
          result_row[[species]] <- if (length(orthologs) > 0) {
            paste(orthologs, collapse = ", ")
          } else {
            "Not available in this species"
          }
        }
        
        return(result_row)
      })
      
      valid_results <- Filter(Negate(is.null), result_list)
      
      # Always update the missing genes reactive value
      missing_genes_rv(missing_genes)
      
      # Always show modal if missing genes exist
      if (length(missing_genes) > 0) {
        showModal(modalDialog(
          title = "Missing Genes",
          paste("The following gene(s) are not found in the ERC dataset:", paste(missing_genes, collapse = ", ")),
          easyClose = TRUE,
          footer = modalButton("OK")
        ))
      }
      
      # If no valid results, return empty df
      if (length(valid_results) == 0) {
        return(data.frame(Gene = character(0)))
      }
      
      final_results <- do.call(rbind, valid_results)
      
      # Rename species short codes to full names
      species_cols <- setdiff(names(final_results), "Gene")
      for (col in species_cols) {
        if (col %in% names(species_full_names)) {
          new_col_name <- species_full_names[[col]]
          names(final_results)[names(final_results) == col] <- new_col_name
        }
      }
      
      return(final_results)
    })
    
    
    
    
    
    # Render results in a table (fixed alignment)
    output$result_table <- DT::renderDataTable({
      df <- results()
      
      # If no results yet, create an empty placeholder
      if (nrow(df) == 0) {
        df <- data.frame(Gene = character(0), stringsAsFactors = FALSE)
      }
      
      DT::datatable(
        df,
        rownames = FALSE,
        extensions = 'Buttons',
        options = list(
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel'),
          pageLength = 10,
          scrollX = TRUE,
          autoWidth = FALSE,  # Prevent auto misalignment
          columnDefs = list(list(width = '150px', targets = "_all")),  # Fix column widths
          drawCallback = JS('function(settings) { this.api().columns.adjust().draw(); }')  # Fix alignment after render
        )
      )
    })
    
    output$result_table_shown <- reactive({
      nrow(results()) > 0
    })
    outputOptions(output, "result_table_shown", suspendWhenHidden = FALSE)
    
    output$missing_genes_text_trigger <- reactive({
      length(missing_genes_rv()) > 0
    })
    outputOptions(output, "missing_genes_text_trigger", suspendWhenHidden = FALSE)
    
    output$missing_genes_box <- renderText({
      mg <- missing_genes_rv()
      if (length(mg) > 0) {
        paste("The following gene(s) are not found in the ERC dataset:", paste(mg, collapse = ", "))
      } else {
        NULL
      }
    })
    
    
    # CSV Download
    output$download_csv <- downloadHandler(
      filename = function() { paste("gene-to-species-results-", Sys.Date(), ".csv", sep = "") },
      content = function(file) { write.csv(results(), file, row.names = FALSE) }
    )
    
    # Excel Download
    output$download_excel <- downloadHandler(
      filename = function() { paste("gene-to-species-results-", Sys.Date(), ".xlsx", sep = "") },
      content = function(file) { write_xlsx(results(), file) }
    )
  })
}