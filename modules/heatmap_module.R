
  # heatmap_module.R
  
  heatmapModuleUI <- function(id) {
    ns <- NS(id)
    tagList(
      htmlOutput(ns("mean_erc"), class = "erc-text"),
      htmlOutput(ns("p_value"), class = "erc-text"),
      
      # Warning for missing genes
      uiOutput(ns("missing_genes_warning_heatmap")),
      
      div(style = "margin-top: 10px;",
          colourpicker::colourInput(ns("low_color"), "Low Color", value = "#FFF5F5"),
          colourpicker::colourInput(ns("high_color"), "High Color", value = "red")
      ),
      
      div(style = "overflow: auto; max-height: 85vh;",
          plotOutput(ns("heatmap"), width = "100%", height = "80vh")
      ),
      
      div(
        style = "font-size: 0.9em; color: #555; margin-top: 10px;",
        "NA (gene pair without an ERC value) are grey and negative ERC values are white. Negative ERC values are not interpretable ",
        "The heatmap is clustered with ERC values above 0. Color scaling is relative to values included in matrix."
      ),
      
      div(style = "margin-top: 10px;",
          downloadButton(ns("download_png"), "Download PNG"),
          downloadButton(ns("download_pdf"), "Download PDF")
      )
    )
  }
  
  
  heatmapModule <- function(id, gene_list, gene_to_cluster_map, erc_matrix, test_results) {
    moduleServer(id, function(input, output, session) {
      
      # -----------------------------
      # Render mean ERC
      # -----------------------------
      output$mean_erc <- renderText({
        req(gene_list())
        input_genes <- trimws(gene_list())
        found_genes <- intersect(input_genes, unique(gene_to_cluster_map$gene))
        if (length(found_genes) < 1) return("")  
        paste("Mean ERC =", round(test_results()$obs, 3))
      })
      
      # -----------------------------
      # Render p-value
      # -----------------------------
      output$p_value <- renderText({
        req(gene_list())
        input_genes <- trimws(gene_list())
        found_genes <- intersect(input_genes, unique(gene_to_cluster_map$gene))
        if (length(found_genes) < 1) return("")  
        paste("Empirical p-value =", round(test_results()$p, 3))
      })
      
      # -----------------------------
      # Missing genes warning
      # -----------------------------
      output$missing_genes_warning_heatmap <- renderUI({
        req(gene_list())
        input_genes <- trimws(gene_list())
        missing_genes <- setdiff(input_genes, unique(gene_to_cluster_map$gene))
        if (length(missing_genes) > 0) {
          HTML(paste0(
            "<div style='color: orange; font-weight: bold;'>",
            "The following gene(s) were not found in the ERC dataset and are excluded: ",
            paste(missing_genes, collapse = ", "),
            "</div>"
          ))
        } else {
          NULL
        }
      })
      
      # -----------------------------
      # Heatmap reactive expression
      # -----------------------------
      heatmap_plot <- reactive({
        req(gene_list())
        input_genes <- trimws(gene_list())
        found_genes <- intersect(input_genes, unique(gene_to_cluster_map$gene))
        
        # -- No genes found --
        if (length(found_genes) == 0) {
          plot.new()
          text(0.5, 0.5, "No input genes were found in the ERC dataset.", cex = 1.2)
          return()
        }
        
        # -- Need at least 3 genes --
        if (length(found_genes) < 3) {
          plot.new()
          text(0.5, 0.5, "At least 3 genes are required to generate a heatmap.", cex = 1.2)
          return()
        }
        
        # Build matrix
        ercmat <- pairmat(found_genes, erc_matrix, gene_to_cluster_map)
        ercmat_round <- round(ercmat, 1)
        reorder_ercmat <- reorder_mymatrix(ercmat_round)
        reorder_ercmat[reorder_ercmat < 0] <- Inf
        
        gene_names <- names(cleanList(found_genes, gene_to_cluster_map, erc_matrix))
        num_genes <- length(gene_names)
        
        if (!is.null(gene_names) && num_genes == ncol(ercmat_round)) {
          rownames(ercmat_round) <- colnames(ercmat_round) <- gene_names
        }
        
        # ----------------------------------
        #  Adaptive scaling rules
        # ----------------------------------
        axis_text_size <- case_when(
          num_genes < 25 ~ 12,
          num_genes < 45 ~ 10,
          num_genes < 70 ~ 8,
          TRUE           ~ 6
        )
        
        # -------------------------------------------------------
        # ULTRA-AGGRESSIVE tile text scaling (very tiny numbers)
        # -------------------------------------------------------
        tile_text_size <- case_when(
          num_genes < 20  ~ 3.0,
          num_genes < 40  ~ 2.0,
          num_genes < 60  ~ 1.5,
          num_genes < 80  ~ 1.0,
          num_genes < 110 ~ 0.9,
          TRUE            ~ 0.7    # SUPER tiny for very large matrices
        )
        
        # Always show tile labels (no cutoff)
        show_tile_text <- TRUE
        
        # ----------------------------------
        # Melt and plot
        # ----------------------------------
        lower_tri <- get_lower_tri(reorder_ercmat)
        melted_data <- reshape2::melt(lower_tri, na.rm = FALSE)
        melted_data$Var1 <- factor(melted_data$Var1, levels = rownames(reorder_ercmat))
        melted_data$Var2 <- factor(melted_data$Var2, levels = colnames(reorder_ercmat))
        
        finite_vals <- melted_data$value[is.finite(melted_data$value)]
        upper_limit <- if(length(finite_vals) == 0) 1 else max(finite_vals, na.rm = TRUE)
        
        low_color <- input$low_color
        high_color <- input$high_color
        
        p <- ggplot(data = melted_data, aes(x = Var1, y = Var2, fill = value)) +
          geom_tile(color = 'white') +
          scale_fill_gradient(low = low_color, high = high_color, na.value = "#d9d9d9",
                              name = "ERC Value", limits = c(0, upper_limit)) +
          theme_minimal() +
          geom_tile(data = subset(melted_data, is.infinite(value)), fill = "white") +
          theme(
            axis.text.x = element_text(size = axis_text_size, angle = 45, hjust = 1, vjust = 1),
            axis.text.y = element_text(size = axis_text_size),
            panel.grid.major = element_blank(),
            axis.ticks = element_blank(),
            text = element_text(size = 15)
          ) +
          labs(x = "Genes", y = "Genes") +
          coord_fixed()
        
        if (show_tile_text) {
          p <- p + geom_text(
            aes(label = ifelse(is.finite(value), round(value, 4), "")),
            size = tile_text_size
          )
        }
        
        p
      })
      
      # -----------------------------
      # Adaptive heatmap height
      # -----------------------------
      output$heatmap <- renderPlot({
        heatmap_plot()
      }, height = function() {
        num_genes <- length(trimws(gene_list()))
        base_height <- 700
        extra <- max(0, (num_genes - 30) * 12)   # 12 px per gene after 30
        base_height + extra
      })
      
      # -----------------------------
      # Download PNG
      # -----------------------------
      output$download_png <- downloadHandler(
        filename = function() paste0("heatmap_", Sys.Date(), ".png"),
        content = function(file) {
          ggsave(file, plot = heatmap_plot(), device = "png", width = 12, height = 12, dpi = 300)
        }
      )
      
      # -----------------------------
      # Download PDF
      # -----------------------------
      output$download_pdf <- downloadHandler(
        filename = function() paste0("heatmap_", Sys.Date(), ".pdf"),
        content = function(file) {
          ggsave(file, plot = heatmap_plot(), device = "pdf", width = 12, height = 12)
        }
      )
      
    })
  }
  
