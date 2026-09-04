library(shiny)
library(DT)
library(bslib)
library(ggplot2)
library(writexl)
library(dplyr)
library(reshape2)
library(readxl)
library(colourpicker)
library(shinycssloaders)

source("modules/utils.R")
source("modules/heatmap_module.R")
source("modules/top_gene_clusters.R")
source("modules/select_species.R")
source("modules/result_matrix_module.R")
source("modules/gene_input_module.R")


# Load data
erc_matrix <- readRDS("data/cluster_modified_to_cluster.rds")
gene_to_cluster_map <- readRDS("data/gene_to_cluster_long.rds")

# Custom Theme 
my_theme <- bs_theme(
  bootswatch = "journal",
  primary = "royalblue"
)

ui <- fluidPage(
  theme = my_theme,
  
  tags$head(
    tags$title("ERC Plasmodium Analysis"),
    tags$style(HTML("
      b, strong {
        font-weight: 700 !important;
      }
    "))
  ),
  
  titlePanel(
    div(
      style = "text-align: center; margin-bottom: 20px;",
      h1(HTML("Evolutionary Rate Covariation (ERC) Analysis in <i>Plasmodium</i>"))
    )
  ),
  
  navlistPanel(
    widths = c(2, 10),
    id = "navlist",
    
    # ---------------- Home Tab ----------------
    tabPanel("Home", icon = icon("home"),
             fluidRow(
               column(12,
                      div(style = "text-align: left; margin-top: 30px; max-width: 900px; margin-left: auto; margin-right: auto;",
                          h2("Home Page"),
                          
                          # --- Section 1: What is an evolutionary rate ---
                          tags$b("What is a protein’s evolutionary rate and why does it vary between species?"),
                          p("An evolutionary rate quantifies the amino acid changes that occurred in a protein along a branch of a phylogeny. 
                            This rate reflects how constrained the protein is. More constrained proteins evolve at slower rates, 
                            while less constrained proteins evolve at faster rates. The rates can also vary between species, as some functional pathways 
                            become more or less constrained."),
                          
                          # --- Section 2: What is ERC ---
                          tags$b("What is evolutionary rate covariation (ERC)?"),
                          p("ERC measures the correlation between the evolutionary rates of two proteins across a phylogeny.
                            Proteins that show high ERC may be functionally related, for example acting in the same pathway or complex,
                            because shifts in evolutionary rate are often more similar between functionally related proteins. In this way, 
                            ERC uses a unique source of information about function that has been shaped by thousands of years of evolution
                            and is orthogonal to lab experiments."),
                          
                          
                          # --- Section 3: How ERC can be used ---
                          tags$b("How can ERC be used to learn about protein function in ", tags$i("Plasmodium"), "?"),
                          p("Many ", tags$i("Plasmodium"), " proteins are not well-studied or annotated. This ERC dataset includes 4,360 genes
                          across 22 ", tags$i("Plasmodium")," species. Gene names from any of these species can be input into the following tools
                          to screen for new functional connects between proteins:"),
                          
                          tags$ol(
                            tags$li(tags$b("Group of genes:"), " Input a set of genes. Outputs ERC values for all pairs as well as an empirical p-value indicating whether the average ERC is significantly higher than in random gene sets."),
                            tags$li(tags$b("Single gene:"), " Input a single gene. Outputs the genes with the highest ERC with this gene from across the genome in a ranked list."),
                            tags$li(tags$b("1:1 Orthologs:"), " Input a single gene. Outputs the names of 1:1 orthologs across species. These are the 1:1 orthologs that we inferred with OrthoMCL and used to estimate evolutionary rates. This is not related to ERC but we thought it was a useful tool.")
                          ),

                          p(tags$b("Importantly, ERC is a screening tool.")),

                          tags$ul(
                            tags$li(
                            "Not all functional connections have an ERC signal, and not all high ERC values can be identified with a particular functional connection."
                            ),
                            tags$li(
                            "It is recommended to first check if your pathway/protein set of interest shows significantly elevated ERC under the ",
                            tags$b("1. Group of genes"),
                            " tab. Then, proceed with ERC to screen for new functional connections."
                            ),
                            tags$li(
                            "ERC is a rank-based approach. There is no cutoff for a high ERC value, and functional connections tend to be enriched at higher ERC values. ",
                            "A practical cutoff is the ERC of your pathway/protein set of interest."
                            ),
                            tags$li(
                            "Paralogous genes are excluded from ERC."
                            )
                          ),
                          
                          br(),
                          
                          # --- Section 4: Citations + Links ---
                          tags$b("If you use this resource, please cite:"), 
                          p(
                            "Hopson HD, Omelianczyk RI, Ramirez A, Little JH, Clark N, Sigala PA, Leffler EM. ",
                            "Evolutionary rate covariation across malaria parasite species enables inference of protein interactions. ",
                            tags$i("Genome Biology and Evolution"),
                            "2026. ",
                            tags$a(
                              "DOI: 10.1093/gbe/evag203",
                              href = "https://doi.org/10.1093/gbe/evag203",
                              target = "_blank"
                            ),
                            " ",
                            tags$a(
                              "PMID: 42664386",
                              href = "https://pubmed.ncbi.nlm.nih.gov/42664386/",
                              target = "_blank"
                            )
                          ),
                          
                          p("The underlying R functions used in this Shiny app are adapted from: ",
                            tags$a("https://github.com/nclark-lab/erc", href = "https://github.com/nclark-lab/erc", target = "_blank")),
                          
                          p("This website is also inspired by the resources available for other species: ",
                            tags$a("https://csbweb.csb.pitt.edu/erc_analysis/", href = "https://csbweb.csb.pitt.edu/erc_analysis/", target = "_blank")),
                          
                          br(),
                          
                          tags$b("For more discussion of ERC and applications in other taxa, please consult:"),
                          tags$div(style = "font-size: 90%; margin-left: 15px;",
                                   p("Clark NL, Alani E, Aquadro CF. ",
                                     "Evolutionary rate covariation reveals shared functionality and co-expression of genes. ", 
                                     tags$i("Genome Research"),
                                     "2012. ",
                                     tags$a(
                                       "DOI: 10.1101/gr.132647.111", 
                                       href = "https://doi.org/10.1101/gr.132647.111",
                                       target= "_blank"
                                     ),
                                     " ",
                                     tags$a(
                                       "PMID: 22287101",
                                       href = "https://pubmed.ncbi.nlm.nih.gov/22287101",
                                       target= "_blank"
                                     )
                                   ),
                                   
                                   p("Little JH, Meyer GH, Grover A, Francette AM, Partha R, Arndt KM, Smith M, Clark N, Chikina M. ",
                                     "ERC2.0 evolutionary rate covariation update improves inference of functional interactions across large phylogenies. ",
                                     tags$i("Genome Research"),
                                     "2025. ",
                                     tags$a(
                                       "DOI: 10.1101/gr.280586.125", 
                                       href = "https://doi.org/10.1101/gr.280586.125",
                                       target= "_blank"
                                     ),
                                     " ",
                                     tags$a(
                                       "PMID: 40774815",
                                       href = "https://pubmed.ncbi.nlm.nih.gov/40774815",
                                       target= "_blank"
                                     )
                                   ),
                                   
                                   p("Findlay GD, Sitnik JL, Wang W, Aquadro CF, Clark NL, Wolfner MF. ",
                                     "Evolutionary Rate Covariation Identifies New Members of a Protein Network Required for Drosophila Female Post-Mating Responses. ",
                                     tags$i("PLOS Genetics"),
                                     "2014. ",
                                     tags$a(
                                       "DOI: 10.1371/journal.pgen.1004108", 
                                       href = "https://doi.org/10.1371/journal.pgen.1004108",
                                       target= "_blank"
                                     ),
                                     " ",
                                     tags$a(
                                       "PMID: 24453993",
                                       href = "https://pubmed.ncbi.nlm.nih.gov/24453993",
                                       target= "_blank"
                                     )
                                   ), 
                                   
                                   p("Little J, Chikina M, Clark NL. ",
                                     "Evolutionary rate covariation is a reliable predictor of co-functional interactions but not necessarily physical interactions. ",
                                     tags$i("eLife"),
                                     "2024. ",
                                     tags$a(
                                       "DOI: 10.7554/eLife.93333", 
                                       href = "https://doi.org/10.7554/eLife.93333",
                                       target= "_blank"
                                     ),
                                     " ",
                                     tags$a(
                                       "PMID: 38415754",
                                       href = "https://pubmed.ncbi.nlm.nih.gov/38415754",
                                       target= "_blank"
                                     )
                                   ),
                                   
                                   p("Kowalczyk A, Gbadamosi O, Kolor K, Sosa J, Andrzejczuk L, Gibson G, St Croix C, Chikina M, Aizenman E, Clark NL, Kiselyov K. ", 
                                     "Evolutionary rate covariation identifies SLC30A9 (ZnT9) as a mitochondrial zinc transporter. ",
                                     tags$i("Biochemical Journal"),
                                     "2021. ",
                                     tags$a(
                                       "DOI: 10.1042/BCJ20210342", 
                                       href = "https://doi.org/10.1042/BCJ20210342",
                                       target= "_blank"
                                     ),
                                     " ",
                                     tags$a(
                                       "PMID: 34397090",
                                       href = "https://pubmed.ncbi.nlm.nih.gov/34397090",
                                       target= "_blank"
                                     )
                                   ),
                                   
                                   p("Raza Q, Choi JY, Li Y, O’Dowd RM, Watkins SC, Chikina M, Hong Y, Clark NL, Kwiatkowski AV. ", 
                                     "Evolutionary rate covariation identifies the GTPase activating protein Raskol as a signaling component of the cadherin adhesion network in Drosophila. ", 
                                     tags$i("PLOS Genetics"),
                                     "2019. ",
                                     tags$a(
                                       "DOI: 10.1371/journal.pgen.1007720", 
                                       href = "https://doi.org/10.1371/journal.pgen.1007720",
                                       target= "_blank"
                                     ),
                                     " ",
                                     tags$a(
                                       "PMID: 30763317",
                                       href = "https://pubmed.ncbi.nlm.nih.gov/30763317/",
                                       target= "_blank"
                                     )
                                   ),
                                   
                                   p("Brunette GJ, Jamalruddin MA, Baldock RA, Clark NL, Bernstein KA. ",
                                     "Evolution-based screening enables genome-wide prioritization and discovery of DNA repair genes. ",
                                     tags$i("Proc Natl Acad Sci U S A"),
                                     "2019 ",
                                     tags$a(
                                       "DOI: 10.1073/pnas.1906559116", 
                                       href = "https://doi.org/10.1073/pnas.1906559116",
                                       target= "_blank"
                                     ),
                                     " ",
                                     tags$a(
                                       "PMID: 31501324",
                                       href = "https://pmc.ncbi.nlm.nih.gov/articles/31501324/",
                                       target= "_blank"
                                     )
                                   )
                          ),
                          
                          # --- Contact & Credits ---
                          tags$hr(),
                          div(style = "margin-top: 20px; font-size: 90%; text-align: left; color: #555;",
                              tags$b("Questions or problems with this website?"),
                              p(
                                "If you encounter any issues, have questions, or suggestions for the website, please contact us at ",
                                a("u60514449@umail.utah.edu", href = "mailto:u60514449@umail.utah.edu"),
                                "."
                              ),
                              br(),
                              tags$b("About this site:"),
                              p(
                                "This ERC Shiny application was created by the Leffler Lab (University of Utah) to support exploration of evolutionary rate covariation in ",
                                tags$i("Plasmodium"), " species."
                              ),
                              p(
                                "For more information about our work, visit ",
                                a("https://lefflerlab.org/", href = "https://lefflerlab.org/", target = "_blank"),
                                "."
                              )
                          ) # end contact div
                      ) # end main div
               ) # end column
             ) # end fluidRow
    ),
    
    # ---------------- Group of Genes Tab ----------------
    tabPanel("Group of genes", icon = icon("chart-line"),
             fluidPage(
               card(
                 full_screen = FALSE,
                 card_header("How to Use This Tool"),
                 p("Given a set of input genes, this tool returns the ERC values for all pairs 
                   and returns an empirical p-value, reflecting whether the mean ERC across all 
                   pairs is higher than random sets of the same size. Results can be viewed and 
                   downloaded as a matrix or clustered heatmap."),
                 p(tags$b("Use case:"), " Determine if the genes in a pathway show an ERC signature in ",
                   tags$i("Plasmodium"), " and identify which gene pairs within the set are driving it."),
                 p(tags$b("Instructions:")),
                 tags$ul(
                   tags$li(HTML("Enter a gene set using gene IDs from any <em>Plasmodium</em> species.")),
                   tags$li("Press ", tags$b("Submit Genes"), " to run the analysis.")
                 ),
                 p(tags$b("Example input gene lists:")),
                 tags$div(style = "margin-left: 15px;",
                          p("PF3D7_1133400, PF3D7_1452000, PF3D7_1116000"),
                          p("or"),
                          p("PBANKA_0915000, PBANKA_1315700, PBANKA_0932000")
                 )
               ),
               
               fluidRow(
                 column(3,
                        div(style = "position: sticky; top: 20px;",
                            card(
                              full_screen = FALSE,
                              card_header("Gene Input"),
                              geneInputModuleUI("gene_input_module_id")
                            )
                        )
                 ),
                 column(9,
                        card(
                          full_screen = TRUE,
                          card_header("Analysis Output"),
                          tabsetPanel(
                            type = "pills",
                            tabPanel("Result Matrix", icon = icon("table"),
                                     withSpinner(resultMatrixModuleUI("result_matrix_module_id"))
                            ),
                            tabPanel("Heatmap", icon = icon("fire"),
                                     withSpinner(heatmapModuleUI("heatmap_module_id"))
                            )
                          )
                        )
                 )
               )
             )
    ),
    
    # ---------------- Single Gene Tab ----------------
    tabPanel("Single gene", icon = icon("star"),
             fluidPage(
               card(
                 card_header("How to Use This Tool"),
                 p("This tool takes a single gene as input and identifies other genes with the highest evolutionary rate covariation (ERC) with it across the genome. 
                   ERC measures the correlation between the evolutionary rates of two proteins across a phylogeny. 
                   High ERC values suggest that the two genes may share a function, act in the same pathway, or participate in the same protein complex."),
                 p("This approach can help uncover potential functional partners for a gene of interest, including genes that have not been experimentally characterized."),
                 p(tags$b("Instructions:")),
                 tags$ul(
                   tags$li("Enter a valid gene ID from any of the ", tags$i("Plasmodium"), " species in the dataset (e.g. PF3D7_0929400.)"),
                   tags$li("Specify the number of top correlated genes you want to retrieve."),
                   tags$li("Click ", tags$strong("Get Top Genes"), " to run the analysis."),
                   tags$li("Download the results using the CSV or Excel buttons for further exploration.")
                 ),
                 p(tags$b("Note:"), " High ERC values suggest possible shared function. Genes listed may include orthologs from other ",
                   tags$i("Plasmodium"), " species."),
                 p(tags$b("Use case:"), " Suppose you are interested in a gene involved in a specific pathway. By entering that gene here, 
                   you can quickly identify potential partner genes that may work with it or belong to the same functional network.")
               ),
               card(
                 full_screen = TRUE,
                 card_header("Single Gene Analysis"),
                 fluidRow(
                   column(12,
                          topGeneClustersModuleUI("top_gene_clusters_module_id")
                   )
                 )
               )
             )
    ),
    
    # ---------------- 1:1 Orthologs Tab ----------------
    tabPanel("1:1 Orthologs", icon = icon("dna"),
             fluidPage(
               card(
                 full_screen = TRUE,
                 card_header("How to Use This Tool"),
                 p("This tool reports 1:1 orthologs for a gene across 22 ", tags$i("Plasmodium"), " species. 
                   Orthologs are genes in different species that evolved from a common ancestral gene and 
                   usually retain the same function. Identifying 1:1 orthologs is useful for comparing 
                   genes across species and understanding evolutionary conservation."),
                 p("Orthologs in this dataset were assigned using OrthoMCL on reference genomes from PlasmoDB (v5). Only 1:1 orthologs are included; paralogs are excluded."),
                 p(tags$b("Instructions:")),
                 tags$ul(
                   tags$li("Enter a gene ID from the dataset (e.g., PF3D7_0929400, PBANKA_0830200, PVP01_0727900)."),
                   tags$li("If the gene is not found in a species, the output will show: ", tags$b("Not found in species"), "."),
                   tags$li("If the gene is not in the dataset, the output will show: ", tags$b("Not available in dataset"), ".")
                 ),
                 p(tags$b("Tip:"), " Start with a well-conserved gene to see complete ortholog results across all species."),
                 p(tags$b("Use case:"), " This tool is useful if you want to explore the evolutionary conservation of a gene across multiple ", 
                   tags$i("Plasmodium"), " species, which can help infer function or identify candidate genes for comparative studies.")
               ),
               div(style = "margin-top: 0px;",
                   card(
                     full_screen = TRUE,
                     card_header("1:1 Orthologs Results"),
                     fluidRow(
                       column(12,
                              div(style = "margin-top: 0px; padding-top: 0px;",
                                  geneToSpeciesModuleUI("gene_to_species_module_id")
                              )
                       )
                     )
                   )
               )
             )
    )
  ), # <- closes navlistPanel

  # --- Floating Memory Panel ---
  tags$div(
    style = "position:fixed; bottom:10px; right:10px; 
             background:#f8f9fa; padding:5px 10px; 
             border:1px solid #ccc; border-radius:5px;
             font-size:12px; z-index:9999;",
    textOutput("memory_usage")
  )
) # <- closes fluidPage


server <- function(input, output, session) {
  gene_list <- geneInputModule("gene_input_module_id")
  
  # Capture test_results from resultMatrixModule
  result_outputs <- resultMatrixModule("result_matrix_module_id", gene_list, gene_to_cluster_map, erc_matrix)
  
  # Pass test_results into heatmapModule
  heatmapModule("heatmap_module_id", gene_list, gene_to_cluster_map, erc_matrix, test_results = result_outputs$test_results)
  
  # Other modules
  topGeneClustersModule("top_gene_clusters_module_id", results_folder = "data/per_gene_rds")
  geneToSpeciesModule("gene_to_species_module_id", gene_to_cluster_map)
  
  # Optional Memory tracker (Uncomment if needed)
  # mem_timer <- reactiveTimer(5000)
  # output$memory_usage <- renderText({
  #   mem_timer()
  #   proc <- ps::ps_handle()
  #   mem <- ps::ps_memory_info(proc)
  #   rss_mb <- as.numeric(mem["rss"]) / (1024^2)
  #   paste0("Memory: ", round(rss_mb, 1), " MB")
  # })
}

shinyApp(ui, server)
