# ====================================================================================
# Design and Implementation: Aditya Sriram, Dr. Soyeon Kim, Dr. Hyun-Jung Park    
# Description: Expression Quantitative Trait Methylation (eQTM) Atlas
# Date: 12/10/25
# R Studio 2025.05.1 Build 513
# R version 4.5.0 
# ====================================================================================
Sys.setlocale("LC_CTYPE", "en_US.UTF-8")   #OS’s UTF-8 locale

# Load necessary libraries for application
x <- c("shiny", "tidyverse", "qqman", "gridExtra", "patchwork", "readxl", "pheatmap", "RColorBrewer", "heatmaply", "plotly", "data.table", "shinythemes", "mailtoR", "locuszoomr", "EnsDb.Hsapiens.v75", "dplyr", "manhattanly", "shinyjs", "bslib")
lapply(x, require, character.only = TRUE)


#Load all the cohort data (I have included a folder for the .Rdata files)
#For uploading to shiny server via crc, please change file paths as needed. Thank you. 

setwd("~/UPDATED_final_eQTM_browser_050125/") #change as necessary for crc app loading
load("~/UPDATED_final_eQTM_browser_050125/12_10_2025_eQTM_final_database_data.Rdata") #change path as necessary for crc app loading
library(dplyr)
library(stringr)


#CSS - round borders, styles, etc.
rounded_border_css <- "
.style-rounded-border {
  border: 2px solid #ddd;
  border-radius: 10px;
  padding: 5px;
  cursor: pointer;
}
.maximized {
  position: fixed;
  top: 50px;
  left: 50px;
  right: 50px;
  bottom: 50px;
  z-index: 9999;
  background-color: white;
  overflow: auto;
}
"


#JS function for making plot/table maximized
maximize_js <- "
$(document).on('click', '.style-rounded-border', function() {
  $(this).toggleClass('maximized');
});
"
sidebar_css <- "
.custom-panel .custom-sidebar {
    width: 250px; /* Adjust the width as needed */
    margin: 10px;
    padding: 10px;
}
"

#UI: START
ui <- fluidPage(
  useShinyjs(),  #Include shinyjs for styling
  theme = bs_theme(
    version = 5,  # Use Bootstrap 5
    bootswatch = "flatly",  #Modern theme (Options: cosmo, minty, lumen, simplex, etc.)
    primary = "#2C3E50",  #Professional dark blue
    secondary = "#4BB543",  #Teal accenting
    success = "#4BB543",
    danger = "#E74C3C"
  ),
  tags$head(
    tags$link(rel = "shortcut icon", href = "favicon.ico"), 
    tags$style(HTML("
    .homepage-intro {
      font-size: 16px;
      line-height: 1.6;
      color: #444;
      max-width: 750px;
      margin: 0 auto;
      text-align: center;
    }

    /* Inline fading segments (NOT block elements) */
    .fade-seg {
      opacity: 0;
      display: inline;
      transition: opacity 0.8s ease-out;
    }

    .fade-seg.visible {
      opacity: 1;
    }
  ")),
    tags$style(HTML("
  .panel-title-align {
    margin-left: 16px;   /* matches panel padding */
    font-size: 24px;
    font-weight: 700;
    margin-bottom: 8px;
  }
")),
    tags$style(HTML("
    /* General card styling */
    .ewas-card,
    .ewas-card-big {
      background: #ffffff;
      border-radius: 10px;
      border: 1px solid #e0e0e0;
      padding: 6px;                    /* reduced padding */
      margin-bottom: 14px;             /* reduced space between cards */
      box-shadow: 0 1px 4px rgba(0,0,0,0.05);
      display: flex;
      justify-content: center;
      align-items: center;
      overflow: hidden;
    }

    /* Larger top plot */
    .ewas-card-big {
      height: 340px;                   /* increased height */
    }

    /* Two smaller plots */
    .ewas-card {
      height: 260px;                   /* increased height */
    }

    /* Make images as large as possible inside the card */
    .ewas-card img,
    .ewas-card-big img {
      width: 100%;
      height: 100%;
      object-fit: contain !important;  /* fits perfectly without overflow */
      margin: 0;                       /* remove default spacing */
      padding: 0;
    }

    /* Reduce extra top space on right panel title */
    .ewas-stats-panel h4 {
      margin-top: 0;
      margin-bottom: 12px;
    }
  ")),
    tags$style(HTML("
  /* Base styling (use whatever you already had and add the lines below) */
  .tissue-panel {
    background: #ffffff;
    border-radius: 12px;
    border: 1px solid #e0e0e0;
    padding: 16px;
    box-shadow: 0 2px 6px rgba(0,0,0,0.06);

    /* CRITICAL: disable any animation on size/box-shadow */
    transition: none !important;
  }

  /* If any global theme adds hover effects, neutralize them */
  .tissue-panel:hover {
    box-shadow: 0 2px 6px rgba(0,0,0,0.06);
    transform: none !important;
  }

  /* Keep title spacing tidy */
  .tissue-panel h4 {
    margin-top: 0;
    margin-bottom: 10px;
  }
    .tissue-panel-title {
    margin-left: 16px;   /* matches panel padding */
    font-size: 24px;
    font-weight: 700;
    margin-bottom: 8px;
  }

  .tissue-summary-row {
    display: flex;
    justify-content: center;
    gap: 10px;
    flex-wrap: wrap;
    margin-bottom: 10px;
  }

  .tissue-summary-card {
    background-color: #ffffff;
    border: 1px solid #dddddd;
    border-radius: 8px;
    padding: 8px 12px;
    min-width: 180px;
    box-shadow: 0 1px 4px rgba(0,0,0,0.08);
  }

  .tissue-summary-card h5 {
    margin: 0 0 3px 0;
    font-size: 14px;
    font-weight: 600;
  }

  .tissue-summary-card p {
    margin: 0;
    font-size: 13px;
    color: #555;
  }

  /* Extra guard: make Plotly's modebar overlay, not reflow */
  .tissue-panel .js-plotly-plot .plotly .modebar {
    position: absolute !important;
    right: 10px;
    top: 10px;
  }
")),
    tags$script(HTML("
    document.addEventListener('DOMContentLoaded', function() {
      var segs = document.querySelectorAll('#homepage_intro_text .fade-seg');
      var delay = 350;  // ms between segments

      segs.forEach(function(seg, idx) {
        setTimeout(function() {
          seg.classList.add('visible');
        }, 300 + idx * delay);
      });
    });
  "))
  ),
  tags$style(HTML("
  /* === EXISTING STYLES === */
  .selectize-dropdown {
    transition: all 0.3s ease-in-out;
    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
    border-radius: 8px;
  }
  input:focus, 
textarea:focus, 
.selectize-input:focus {
  outline: none !important;
  border: 2px solid rgba(58, 95, 180, 0.35) !important;
  box-shadow: 0 0 10px rgba(58, 95, 180, 0.2) !important;
  background-color: rgba(255, 255, 255, 0.3) !important;
  backdrop-filter: blur(10px);
  -webkit-backdrop-filter: blur(10px);
  border-radius: 10px !important;
  transition: all 0.3s ease-in-out;
}
input, 
textarea, 
.selectize-input, 
.select:focus, 
.selectize-control.single .selectize-input.focus {
  transition: border 0.1s ease, box-shadow 0.1s ease, background-color 0.1s ease;
}


.selectize-control.single .selectize-input.focus {
  outline: none !important;
  border: 2px solid rgba(58, 95, 180, 0.4) !important;
  box-shadow: 0 0 10px rgba(58, 95, 180, 0.2) !important;
  background-color: rgba(255, 255, 255, 0.3) !important;
  backdrop-filter: blur(8px);
  -webkit-backdrop-filter: blur(8px);
  border-radius: 10px;
  transition: all 0.1s ease;
}

select:focus {
  outline: none !important;
  border: 2px solid rgba(58, 95, 180, 0.4) !important;
  box-shadow: 0 0 10px rgba(58, 95, 180, 0.2) !important;
  border-radius: 10px;
  transition: all 0.1s ease;
}

button:focus, 
.btn:focus, 
.btn:active:focus, 
.btn-primary:focus, 
.btn-danger:focus, 
.btn-secondary:focus,
.action-button:focus {
  outline: none !important;
  box-shadow: none !important;
  border-color: rgba(58, 95, 180, 0.4) !important; /* Optional subtle blue border */
  background-color: rgba(255, 255, 255, 0.3) !important;
  backdrop-filter: blur(10px);
  -webkit-backdrop-filter: blur(10px);
  transition: all 0.15s ease-in-out;
}

  .selectize-dropdown-content .option {
    transition: background-color 0.2s ease-in-out, padding 0.2s ease-in-out;
    padding: 6px 10px;
  }

  .selectize-dropdown-content .option:hover {
    padding-left: 12px; /* slight shift on hover */
  }
  .custom-table-container {
    display: flex;
    justify-content: center;
    align-items: center;
    width: 100%;
    margin-top: 600px;
  }

  .nav-tabs > li > a {
    color: black !important;
    font-size: 18px !important;
  }

  .description-text {
    font-size: 16px !important;
  }
  
  .tissue-panel {
  text-align: center;
  padding: 10px;
  transition: transform 0.2s ease-in-out;
  position: relative;
  overflow: hidden;
  border-radius: 8px;
  background-color: #ffffff;
  box-shadow: 0px 2px 4px rgba(0, 0, 0, 0.1);
  margin: 6px;
  font-size: 14px;
}

.tissue-panel:hover {
  transform: scale(1.03);
  background-color: #f9f9f9;
}

.tissue-panel .panel-content {
  z-index: 1;
  position: relative;
}

.tissue-panel .description-overlay {
  position: absolute;
  top: 0;
  left: 0;
  width: 100%;
  height: 100%;
  background-color: rgba(255, 255, 255, 0.85);
  opacity: 0;
  transition: opacity 0.25s ease-in-out;
  z-index: 2;
  display: flex;
  justify-content: center;
  align-items: center;
  padding: 6px;
}

.tissue-panel:hover .description-overlay {
  opacity: 1;
}

.tissue-panel h4 {
  margin: 5px 0 2px 0;
  font-size: 16px;
  font-weight: bold;
}

.tissue-panel p {
  margin: 0;
  font-size: 13px;
}

  .left-align-text {
    text-align: left !important;
  }


  /* === IMAGE SIZING === */
.homepage-panel img {
    width: 135px !important;
    height: auto !important;
    display: block;
    margin: 0 auto 5px auto !important;  /* Reduce space below the image */
    background: transparent !important;
  }
  
  /* Reduce extra space under the image */
  .homepage-panel h3 {
    margin-top: 0px !important;
    padding-top: 0px !important;
  }

  /* === PANEL STYLING WITH SMALL SPACING === */
.homepage-panel {
  text-align: center;
  padding: 15px;
  transition: transform 0.3s ease-in-out;
  position: relative;
  overflow: hidden;
  border-radius: 10px;
  background-color: #ffffff; /* restore solid white */
  box-shadow: 0px 4px 6px rgba(0, 0, 0, 0.1);
  margin: 5px;
}

  .homepage-panel:hover {
    transform: scale(1.05);
  }

  /* === BLUR EFFECT ON HOVER === */
  .panel-content {
    transition: filter 0.3s ease-in-out;
    position: relative;
    z-index: 1;
  }

  .homepage-panel:hover .panel-content {
    filter: blur(8px);  /* Strong blur so text stands out */
  }

  /* === TEXT OVERLAY === */
  .description-overlay {
    position: absolute;
    top: 0;
    left: 0;
    width: 100%;
    height: 100%;
    background-color: rgba(255, 255, 255, 0.85);  /* Semi-transparent white for contrast */
    opacity: 0;
    transition: opacity 0.3s ease-in-out;
    z-index: 2;
    display: flex;
    justify-content: center;
    align-items: center;
    padding: 10px;
  }

  .homepage-panel:hover .description-overlay {
    opacity: 1;  /* Show overlay when hovering */
  }

  .description-text {
    font-size: 16px;
    font-weight: lighter;
    color: black;
    text-align: center;
    z-index: 3;
  }
.btn {
  background-color: transparent !important;
  color: #3a5f8a !important;
  border: 2px solid #3a5f8a !important;
  border-radius: 8px;
  font-weight: bold;
  font-family: 'Segoe UI', sans-serif;
  box-shadow: 0 2px 4px rgba(0,0,0,0.15);
  transition: all 0.3s ease;
}

.btn:hover {
  background-color: rgba(106, 90, 205, 0.1) !important;
  color: #3a5f8a !important;
}
  /* --- Smoother, switch-style tab design --- */
  .nav-tabs {
    border-bottom: none;
    margin-bottom: 0;
  }

  .nav-tabs > li > a {
    background-color: #f9f9f9;
    color: #3a5f8a;
    border: none;
    border-radius: 10px 10px 0 0;
    margin-right: 5px;
    padding: 10px 20px;
    font-weight: 500;
    box-shadow: inset 0 0 0 rgba(0,0,0,0);  /* Reset */
    transition: all 0.25s ease-in-out;
  }

  .nav-tabs > li > a:hover {
    background-color: #eaeaea;
    box-shadow: inset 0 2px 5px rgba(0, 0, 0, 0.1);
  }

  .nav-tabs > li.active > a,
  .nav-tabs > li.active > a:focus,
  .nav-tabs > li.active > a:hover {
    background-color: #ffffff;
    color: #2c3e50;
    box-shadow: inset 0 -3px 5px rgba(0, 0, 0, 0.1), 0 1px 2px rgba(0, 0, 0, 0.05);
    font-weight: bold;
    transform: translateY(1px);  /* Pressed look */
  }

  .tab-content {
    padding: 25px;
    background-color: #ffffff;
    border-radius: 0 0 10px 10px;
    box-shadow: 0 4px 8px rgba(0, 0, 0, 0.06);
    border-top: 1px solid #dee2e6;
  }
  html, body {
    background: linear-gradient(120deg, rgba(255, 255, 255, 0.08), rgba(230, 230, 255, 0.06));
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif;
    color: #1c1c1e;
    overflow-x: hidden;
    scroll-behavior: smooth;
    -webkit-font-smoothing: antialiased;
  }

  .container-fluid, .tab-content, .mainPanel, .sidebarPanel, .panel, .well, .custom-glass-panel {
    background: linear-gradient(135deg, rgba(255,255,255,0.22) 0%, rgba(245,245,255,0.18) 100%);
    border-radius: 24px;
    border: 1px solid rgba(255,255,255,0.25);
    box-shadow: 0 8px 30px rgba(0,0,0,0.06), inset 0 0 0.5px rgba(255,255,255,0.4);
    backdrop-filter: blur(30px) saturate(180%);
    -webkit-backdrop-filter: blur(30px) saturate(180%);
    transition: all 0.3s ease-in-out;
    padding: 20px;
    position: relative;
    z-index: 1;
  }

  .panel:hover {
    transform: scale(1.002);
    box-shadow: 0 12px 36px rgba(0,0,0,0.12);
  }

  .btn, .download-button, .action-button {
    background: rgba(255,255,255,0.25);
    border-radius: 12px;
    border: 1px solid rgba(255,255,255,0.25);
    color: #1c1c1e;
    font-weight: 500;
    box-shadow: inset 0 0 0.5px rgba(255,255,255,0.6);
    backdrop-filter: blur(10px);
    transition: all 0.3s ease;
  }

  .btn:hover, .action-button:hover {
    background: rgba(255,255,255,0.35);
    box-shadow: none;
  }

  .navbar, .navbar-default {
    background-color: rgba(255,255,255,0.2);
    backdrop-filter: blur(20px);
    border-bottom: 1px solid rgba(255,255,255,0.15);
    box-shadow: 0 2px 8px rgba(0,0,0,0.04);
  }

.homepage-panel:hover {
  transform: scale(1.01);
  box-shadow: 0 12px 36px rgba(0, 0, 0, 0.07);  /* Slight lift on hover */
  background: rgba(255, 255, 255, 0.94);  /* Brighten just slightly on hover */
}
.btn:not(#submit_homepage_search):not(.download-button) {
  margin-bottom: 5px !important;
}


")),
  tags$style(HTML("
  .selectize-dropdown-content .option.active {
    background-color: rgba(58, 95, 180, 0.25) !important;  /* richer blue, subtle opacity */
    color: black !important;
  }
")),
  div(
    class = "custom-panel",
    titlePanel("The eQTM Atlas"),
    tabsetPanel(
      id = "main_tabs",
      tabPanel(
        title = "Home Page",
        fluidRow(
          column(12, align = "left",
                 a(href = "https://docs.google.com/document/d/1-AgP6dSbxAn85UGntj1v8DzitnkCU_ao/edit?usp=sharing&ouid=113345651194344459277&rtpof=true&sd=true", target = "_blank", "Home Page Overview")
          )
        ),
        fluidRow(
          column(12, 
                 div(style = "text-align: center; margin-bottom: 20px;", 
                     h1("Welcome to the eQTM Atlas.")),
                 div(
                   style = "max-width: 750px; margin: 0 auto; text-align: center; margin-top: -5px; margin-bottom: 25px; font-size: 16px; color: #444;",
                   p(
                     id = "homepage_intro_text",
                     class = "homepage-intro",
                     
                     tags$span(
                       class = "fade-seg",
                       "The expression quantitative trait methylation (eQTM) Atlas provides valuable data, informative genome maps, and other relevant information about eQTM genes that are significantly associated with various CpGs, some of which may have regulatory effects on gene expression. "
                     ),
                     
                     tags$span(
                       class = "fade-seg",
                       "The goal of eQTM analysis is to identify CpG sites that are significantly associated with genes and gene expression, or identify genes that are significantly associated with CpG sites. "
                     ),
                     
                     tags$span(
                       class = "fade-seg",
                       "Our database is also linked to the EWAS Atlas, allowing for cross-referencing across two powerful resources for methylation data."
                     )
                   )
                 ),
                 div(
                   style = "width: 600px; margin: 0 auto;",
                   div(
                     class = "input-group input-group-lg",
                     
                     # Text input
                     tags$input(
                       id = "homepage_search",
                       type = "text",
                       class = "form-control",
                       placeholder = "Search a gene or CpG (e.g. PDCD10, cg05575921)..."
                     ),
                     
                     # Actual actionButton (visually styled as just an icon)
                     actionButton(
                       inputId = "submit_homepage_search",
                       label = NULL,
                       icon = icon("magnifying-glass"),
                       class = "btn",
                       style = "
    background-color: #fff;
    border: 1px solid #ced4da;
    border-left: none;
    border-radius: 0 8px 8px 0;
    padding: 10px 15px;
    color: #6c757d;"
                     )
                   )
                 ),
                 
                 # Examples row
                 div(
                   style = "text-align: center; margin-top: 8px;",
                   tags$span("Examples: "),
                   actionLink("example_gene1", "PDCD10"), tags$span(", "),
                   actionLink("example_gene2", "USP6"), tags$span(", "),
                   actionLink("example_cpg1", "cg04983687"), tags$span(", "),
                   actionLink("example_cpg2", "cg00119778")
                 ),
                 div(id = "homepage_search_loading", style = "margin-top: 10px; display: none;",
                     tags$span("Searching...", style = "color: #888; font-style: italic;"))
                 
          )
        ),
        uiOutput("home_content")  #This now loads the six blur panels instead of the single static homepage.
      ),
      # tabPanel("Data Overview",
      #          fluidRow(
      #            column(12, align = "left",
      #                   a(href = "https://drive.google.com/file/d/1uFrylAztmlY0q1bWsTVfXPjP8FcbBwc8/view?usp=sharing", target = "_blank", "Guide for Data Overview Tab")
      #            )
      #          ),
      #          fluidRow(
      #            column(12, p("This page allows you to explore, filter and search across all data. You can filter by cohort if needed.", class = "left-align-text"))
      #          ),
      #          selectInput("cohort_filter", "Select Cohort:", choices = c("All", unique(combined_df_atopic_asthma_final$cohort))),
      #          mainPanel(
      #            DT::dataTableOutput("data_table_home")
      #          )
      # ),
      tabPanel(
        "Gene ↔ CpG Search",
        fluidRow(
          column(12, align = "left",
                 a(href = "https://docs.google.com/document/d/1lHYjb6O2DxYm7cdfCaYqONPhU-CBesFD/edit?usp=sharing&ouid=113345651194344459277&rtpof=true&sd=true", target = "_blank", "Gene ↔ CpG Usage Guide")
          )
        ),
        fluidRow(
          column(
            6,
            textAreaInput("gene_cpg_list", "Enter a list of comma-separated Genes or CpGs:", 
                          placeholder = "e.g., PDCD10, USP6", rows = 10),
            radioButtons("input_type", "Input Type:", choices = c("Gene", "CpG"), 
                         selected = "Gene", inline = TRUE),
            actionButton("search_gene_cpg", "Search", class = "btn-primary"),
            actionButton("load_gene_example_button", "Load Genes Example", class = "btn-secondary"),
            actionButton("load_cpg_example_button", "Load CpGs Example", class = "btn-secondary")
          ),
          column(
            6,
            selectInput("filter_cohort", "Filter by Cohort:", choices = NULL, 
                        selected = "All"),
            selectInput("filter_tissue", "Filter by Tissue:", choices = NULL, 
                        selected = "All"),
            downloadButton("download_gene_cpg_results", "Download Results", class = "download-button")
          )
        ),
        fluidRow(
          DT::dataTableOutput("gene_cpg_table")
        )
      ),
      tabPanel("Genome Browser",
               fluidRow(
                 column(12, align = "left",
                        a(href = "https://docs.google.com/document/d/1neFQyI70kKpD3JPlgEs6873UJ_gjNS0E/edit?usp=sharing&ouid=113345651194344459277&rtpof=true&sd=true", target = "_blank", "Genome Browser Usage Guide")
                 )
               ),
               fluidRow(
                 column(12, p("This page lets you view eQTMs and CpGs in the context of their genomic locations. Search by gene name or CpG probe ID. You can filter by cohort or by disease.", class = "left-align-text"))
               ),
               sidebarLayout(
                 sidebarPanel(
                   class = "custom-sidebar",
                   textInput("gene_input", "Enter Gene:", "", placeholder = "ex. PDCD10"),
                   textInput("probeid_input", "Enter Probe ID:", "", placeholder = "ex. cg04983687"),
                   selectInput("cohort_input", "Select Cohort:", 
                               choices = c("All", unique(combined_df_atopic_asthma_final$cohort)), selected = "All"),
                   selectInput("disease_input", "Select Disease:", 
                               choices = c("All", unique(combined_df_atopic_asthma_final$disease)), selected = "All"),
                   actionButton("search_button", "Search"),
                   actionButton("soft_reset", "Reset Search", icon = icon("refresh"), class = "btn-danger"),
                   actionButton("load_example_button", "Load Gene Example"),  # Add this line
                   actionButton("load_cpg_example_button", "Load CpG Example"),  # Add this line
                   uiOutput("download_button_ui")
                 ),
                 mainPanel(
                   conditionalPanel(
                     condition = "output.manhattanPlotAvailable",
                     div(style = "margin-bottom: 40px;", 
                         plotlyOutput("manhattan_plot"), 
                         class = "style-rounded-border")
                   ),
                   conditionalPanel(
                     condition = "output.transManhattanPlotAvailable",
                     div(style = "margin-bottom: 40px;", 
                         plotlyOutput("trans_manhattan_plot"), 
                         class = "style-rounded-border")
                   ),
                   conditionalPanel(
                     condition = "output.cpgManhattanPlotAvailable",
                     div(style = "margin-bottom: 40px;", 
                         plotlyOutput("cpg_manhattan_plot"), 
                         class = "style-rounded-border")
                   ),
                   conditionalPanel(
                     condition = "output.transManhattanPlotByCpgAvailable",
                     div(style = "margin-bottom: 40px;", 
                         plotlyOutput("trans_manhattan_plot_by_cpg"), 
                         class = "style-rounded-border")
                   ),
                   div(style = "margin-top: 20px;", 
                       div(style = "overflow-x: auto;", DT::dataTableOutput("result_table"))
                   )
                 )
               )
      ),
      tabPanel("Heatmap",
               fluidRow(
                 column(12, align = "left",
                        a(href = "https://docs.google.com/document/d/163zsQ738XsMTRExX1nPLR8QlCy-z6EnQ/edit?usp=sharing&ouid=113345651194344459277&rtpof=true&sd=true", target = "_blank", "Heatmap Usage Guide")
                 )
               ),
               fluidRow(
                 column(12, p("This page lets you view a heatmap of p-values for data points across various tissues. Scroll to the right of the heatmap to view the legend and other plot controls.", class = "left-align-text"))
               ),
               sidebarLayout(
                 sidebarPanel(
                   class = "custom-sidebar",
                   textInput("heatmap_gene_input", "Enter Gene for Heatmap:", placeholder = "ex. PDCD10"),
                   textInput("heatmap_probeid_input", "Enter Probe ID for Heatmap:", placeholder = "ex. cg04983687"),
                   selectInput("heatmap_cohort_input", "Select Cohort for Heatmap:", 
                               choices = c("All", unique(combined_df_atopic_asthma_final$cohort)), selected = "All"),
                   actionButton("search_heatmap_button", "Search", class = "btn-primary"),
                   actionButton("reset_heatmap_button", "Reset Search", icon = icon("refresh"), class = "btn-danger"),
                   actionButton("load_heatmap_gene_example_button", "Load Gene Example", class = "btn-secondary"),
                   actionButton("load_heatmap_cpg_example_button", "Load CpG Example", class = "btn-secondary")
                   
                   
                 ),
                 fluidPage(
                   div(style = "overflow-x: scroll; width: 75%; height: 500px;", 
                       plotlyOutput("gene_heatmap", width = "100%", height = "100%"))
                 )
               )
      ),
      # tabPanel("EWAS-eQTM Co-visualization",
      #          fluidRow(
      #            column(12, align = "left",
      #                   a(href = "https://drive.google.com/file/d/18NkVV30wLI0Pwii_Chnct-U7VC7CJYcn/view?usp=sharing", target = "_blank", "Guide for EWAS-eQTM Tab")
      #            )
      #          ),
      #          fluidRow(
      #            column(12, p("This page lets you visualize EWAS-eQTM data. You can upload your own two-column, tab-delimited text formatted EWAS data.", class = "left-align-text"))
      #          ),
      #          fluidRow(
      #            column(4,
      #                   textInput("coloc_gene_input", "Gene:", value = ""),
      #                   selectInput("coloc_cohort_input", "Cohort:", choices = c("All", unique(combined_df_atopic_asthma_final$cohort))),
      #                   selectInput("coloc_tissue_input", "Tissue:", choices = c("All", unique(combined_df_atopic_asthma_final$tissue))),
      #                   actionButton("coloc_plot_button", "Plot", class = "btn-primary"),
      #                   actionButton("example_button", "Load Example", class = "btn-secondary"),
      #                   textAreaInput("ewas_data_text", "Paste your EWAS data (cpgID, pValue):", "", rows = 10, placeholder = "cpgID\tpValue\n...")
      #            ),
      #            column(8,
      #                   div(style = "display: flex; justify-content: center; align-items: center;",
      #                       plotlyOutput("colocalization_plot", width = "100%", height = "500px")
      #                   )
      #            )
      #          )
      # ),
      tabPanel("EWAS-eQTM Integration",
               fluidRow(
                 column(12, align = "left",
                        a(href = "https://docs.google.com/document/d/1-k6lmYw4gYthDbItmAa_hcnfZzquVPFz/edit?usp=sharing&ouid=113345651194344459277&rtpof=true&sd=true", target = "_blank", "EWAS-eQTM Usage Guide")
                 )
               ),
               fluidRow(
                 column(12, p("This page lets you view EWAS Atlas associations in tandem with our eQTM database cohort data.", class = "left-align-text"))
               ),
               sidebarLayout(
                 sidebarPanel(
                   textInput("ewas_probe_input", "Enter CpG ID:", placeholder = "e.g., cg05575921"),
                   # textInput("ewas_gene_input", "Enter Gene Symbol:", placeholder = "e.g., PDCD10"),
                   selectizeInput(
                     "ewas_trait_input", 
                     "Enter Trait Keyword:",
                     choices = c("", sort(unique(iconv(ewas_atlas$trait, from = "", to = "UTF-8", sub = "")))),
                     selected = NULL,
                     multiple = FALSE,
                     options = list(
                       placeholder = 'e.g., Asthma',
                       create = FALSE
                     )
                   ), # allows typing even if not in list
                   selectInput("ewas_cohort_input", "eQTM Cohort (Required):", 
                               choices = c("Select a cohort", unique(combined_df_atopic_asthma_final$cohort)), 
                               selected = "Select a cohort"),
                   
                   actionButton("search_ewas", "Search EWAS Atlas"),
                   actionButton("load_ewas_example_button", "Load Example (trait optional)", class = "btn-secondary"),
                   actionButton("reset_ewas_inputs", "Reset Search", icon = icon("refresh"), class = "btn-danger")
                 ),
                 mainPanel(
                   uiOutput("ewas_headers")  # Dynamically generated tabbed results
                 )
               )
      ),
      tabPanel("Downloads",
               uiOutput("cohortTableUI")
      ),
      tabPanel("Upload Your Data/Site Feedback", 
               mailtoR(email = "ads303@pitt.edu",
                       text = "Please click here to send any issues you are having with the site. 
                       We are still working on development. 
                       Additionally, if you want to upload your own eQTM data to our data catalog, please click here to contact us about uploading your data.", 
                       subject = "eQTM Atlas help/upload - (insert specific issue here)"),
               use_mailtoR())
    )
  )
)

#Server: START
server <- function(input, output, session) {
  reset_state <- reactiveVal(FALSE)
  ewas_reset_state <- reactiveVal(FALSE)
  heatmap_reset_state <- reactiveVal(FALSE)
  output$cohortTableUI <- renderUI({
    # Move all EVA-PR cohorts to the top in desired order
    cohorts <- c(
      "EVA-PR: Atopic Asthma Only",
      "EVA-PR: Controls",
      "EVA-PR: White Blood Cell",
      "EVA-PR: All Subjects",
      setdiff(unique(combined_df_atopic_asthma_final$cohort), 
              c("EVA-PR: Atopic Asthma Only", "EVA-PR: Controls", 
                "EVA-PR: White Blood Cell", "EVA-PR: All Subjects"))
    )
    #change name of EVA-PR to EVA-PR: Atopic asthma + healthy controls
    #Manually set FDR threshold values for each cohort; have to update every time there is new data added
    original_cohorts <- unique(combined_df_atopic_asthma_final$cohort)
    # Add author/publication links
    # Add author/publication links
    author_lookup <- setNames(c(
      '<a href="https://www.jacionline.org/article/S0091-6749(23)00713-3/fulltext" target="_blank">Kim et al. - JACI - 2023</a>',
      '<a href="https://www.jacionline.org/article/S0091-6749(23)00713-3/fulltext" target="_blank">Kim et al. - JACI - 2023</a>',
      '<a href="https://www.jacionline.org/article/S0091-6749(23)00713-3/fulltext" target="_blank">Kim et al. - JACI - 2023</a>',
      '<a href="https://www.sciencedirect.com/science/article/pii/S0012369220317165?via%3Dihub" target="_blank">Kim et al. - Chest - 2020</a>',
      '<a href="https://www.jacionline.org/article/S0091-6749(23)00713-3/fulltext" target="_blank">Kim et al. - JACI - 2023</a>',
      '<a href="https://www.pnas.org/doi/abs/10.1073/pnas.1814263116" target="_blank">Taylor et al. - PNAS - 2019</a>',
      '<a href="https://www.nature.com/articles/s41598-023-39936-3" target="_blank">Keshawarz et al. - Sci Rep - 2023</a>',
      '<a href="https://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-021-01148-9" target="_blank">Diez-Villanueva et al. - Clin Epigenetics - 2021</a>',
      '<a href="https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1601-3" target="_blank">Husquin et al. - Genome Biol - 2018</a>',
      '<a href="https://elifesciences.org/articles/65310" target="_blank">Ruiz Arenas et al. - eLife - 2022</a>',
      '<a href="https://www.nature.com/articles/s41588-022-01248-z" target="_blank">Oliva et al. - Nat Genet - 2023</a>',
      '<a href="https://www.nature.com/articles/s41588-022-01248-z" target="_blank">Oliva et al. - Nat Genet - 2023</a>',
      '<a href="https://www.nature.com/articles/s41588-022-01248-z" target="_blank">Oliva et al. - Nat Genet - 2023</a>',
      '<a href="https://www.nature.com/articles/s41588-022-01248-z" target="_blank">Oliva et al. - Nat Genet - 2023</a>',
      '<a href="https://www.nature.com/articles/s41588-022-01248-z" target="_blank">Oliva et al. - Nat Genet - 2023</a>',
      '<a href="https://www.nature.com/articles/s41588-022-01248-z" target="_blank">Oliva et al. - Nat Genet - 2023</a>',
      '<a href="https://www.nature.com/articles/s41588-022-01248-z" target="_blank">Oliva et al. - Nat Genet - 2023</a>',
      '<a href="https://www.nature.com/articles/s41588-022-01248-z" target="_blank">Oliva et al. - Nat Genet - 2023</a>'
    ), cohorts)
    
    
    
    threshold_lookup <- setNames(c(
      "FDR < 0.05", 
      "FDR < 0.01", 
      "cis: p < 1E-7, trans p < 1E-14",
      "p < 1.7e−08", 
      "p < 4.99e-05",
      "FDR < 0.05", 
      "FDR < 0.01", 
      "FDR < 0.05",
      "p < 0.05", 
      "Bonferroni-adjusted p < 0.05, FDR < 0.05",
      "Bonferroni-adjusted p < 0.05, FDR < 0.05",
      "Bonferroni-adjusted p < 0.05, FDR < 0.05",
      "Bonferroni-adjusted p < 0.05, FDR < 0.05",
      "Bonferroni-adjusted p < 0.05, FDR < 0.05",
      "Bonferroni-adjusted p < 0.05, FDR < 0.05",
      "Bonferroni-adjusted p < 0.05, FDR < 0.05",
      "Bonferroni-adjusted p < 0.05, FDR < 0.05",
      "FDR < 0.05"
    ), original_cohorts)
    # Reordered thresholds
    threshold_values <- threshold_lookup[cohorts]
    author_values <- author_lookup[cohorts]
    
    
    
    table_data <- data.frame(
      Cohort = cohorts,
      Author = author_values,
      Disease = sapply(cohorts, function(cohort) {
        combined_df_atopic_asthma_final$disease[which(combined_df_atopic_asthma_final$cohort == cohort)[1]]
      }),
      NumCases = sapply(cohorts, function(cohort) {
        combined_df_atopic_asthma_final$num_cases[which(combined_df_atopic_asthma_final$cohort == cohort)[1]]
      }),
      NumControls = sapply(cohorts, function(cohort) {
        combined_df_atopic_asthma_final$num_controls[which(combined_df_atopic_asthma_final$cohort == cohort)[1]]
      }),
      Tissue = sapply(cohorts, function(cohort) {
        combined_df_atopic_asthma_final$tissue[which(combined_df_atopic_asthma_final$cohort == cohort)[1]]
      }),
      Threshold = threshold_values,
      stringsAsFactors = FALSE
    )
    
    
    #Create safe IDs by replacing spaces and special characters with underscores
    safe_ids <- gsub("[^a-zA-Z0-9]", "_", cohorts)
    
    #Set up download handlers for each cohort
    lapply(seq_along(cohorts), function(i) {
      cohort <- cohorts[i]
      cohort_data_id <- paste0("download_", safe_ids[i])
      output[[cohort_data_id]] <- downloadHandler(
        filename = function() { paste(cohort, "data.csv", sep = "_") },
        content = function(file) {
          write.csv(combined_df_atopic_asthma_final[combined_df_atopic_asthma_final$cohort == cohort, ], file, row.names = FALSE)
        }
      )
    })
    
    #Create UI elements for each cohort
    cohort_ui <- lapply(seq_len(nrow(table_data)), function(i) {
      cohort <- paste(table_data$Cohort[i])
      cohort_data_id <- paste0("download_", safe_ids[i])
      
      tags$tr(
        tags$td(cohort, style = "padding-right: 20px; border-bottom: 1px solid #ddd;"),
        tags$td(HTML(table_data$Author[i]), style = "padding-right: 20px; border-bottom: 1px solid #ddd;"),
        tags$td(table_data$Disease[i], style = "padding-right: 20px; border-bottom: 1px solid #ddd;"),
        tags$td(table_data$NumCases[i], style = "padding-right: 20px; border-bottom: 1px solid #ddd;"),
        tags$td(table_data$NumControls[i], style = "padding-right: 20px; border-bottom: 1px solid #ddd;"),
        tags$td(table_data$Tissue[i], style = "padding-right: 20px; border-bottom: 1px solid #ddd;"),
        tags$td(table_data$Threshold[i], style = "padding-right: 20px; border-bottom: 1px solid #ddd;"), # Add Threshold column
        tags$td(downloadButton(cohort_data_id, "Download", class = "download-button"), style = "padding-right: 20px; border-bottom: 1px solid #ddd;")
      )
    })
    
    #Create the table UI with a border and enough spacing
    tags$table(
      style = "border-collapse: collapse; width: 100%;",
      tags$thead(
        tags$tr(
          tags$th("Cohort", style = "border-bottom: 2px solid #ddd; padding-bottom: 5px;"),
          tags$th("Author", style = "border-bottom: 2px solid #ddd; padding-bottom: 5px;"),
          tags$th("Disease", style = "border-bottom: 2px solid #ddd; padding-bottom: 5px;"),
          tags$th("Num Cases", style = "border-bottom: 2px solid #ddd; padding-bottom: 5px;"),
          tags$th("Num Controls", style = "border-bottom: 2px solid #ddd; padding-bottom: 5px;"),
          tags$th("Tissue", style = "border-bottom: 2px solid #ddd; padding-bottom: 5px;"),
          tags$th("Threshold", style = "border-bottom: 2px solid #ddd; padding-bottom: 5px;"), # Add Threshold header
          tags$th("Download", style = "border-bottom: 2px solid #ddd; padding-bottom: 5px;")
        )
      ),
      tags$tbody(
        do.call(tagList, cohort_ui)
      )
    )
  })
  
  
  observeEvent(input$submit_homepage_search, {
    shinyjs::show("homepage_search_loading")
    
    updateTabsetPanel(session, inputId = "main_tabs", selected = "Gene ↔ CpG Search")
    
    shinyjs::delay(250, {
      query <- input$homepage_search
      
      if (!is.null(query) && nzchar(query)) {
        # Detect CpG and switch radio button using Shiny, not JS
        new_input_type <- if (grepl("^cg\\d{5,8}$", query, ignore.case = TRUE)) "CpG" else "Gene"
        updateRadioButtons(session, "input_type", selected = new_input_type)
        
        shinyjs::delay(250, {
          updateTextAreaInput(session, "gene_cpg_list", value = query)
          click("search_gene_cpg")
          shinyjs::delay(500, shinyjs::hide("homepage_search_loading"))
        })
      } else {
        shinyjs::hide("homepage_search_loading")
        showNotification("Please enter a search term.", type = "error")
      }
    })
  })
  
  
  
  
  
  observeEvent(input$example_gene1, {
    updateTextInput(session, "homepage_search", value = "PDCD10")
  })
  observeEvent(input$example_gene2, {
    updateTextInput(session, "homepage_search", value = "USP6")
  })
  observeEvent(input$example_cpg1, {
    updateTextInput(session, "homepage_search", value = "cg04983687")
  })
  observeEvent(input$example_cpg2, {
    updateTextInput(session, "homepage_search", value = "cg00119778")
  })
  
  observeEvent(input$homepage_search, {
    # optionally debounce or monitor typing
  })
  
  
  
  observeEvent(input$soft_reset, {
    updateTextInput(session, "gene_input", value = "")
    updateTextInput(session, "probeid_input", value = "")
    updateSelectInput(session, "cohort_input", selected = "All")
    updateSelectInput(session, "disease_input", selected = "All")
    
    reset_state(TRUE)  # Set flag to TRUE to indicate reset
  })
  
  observeEvent(input$reset_heatmap_button, {
    updateTextInput(session, "heatmap_gene_input", value = "")
    updateTextInput(session, "heatmap_probeid_input", value = "")
    updateSelectInput(session, "heatmap_cohort_input", selected = "All")
    
    heatmap_reset_state(TRUE)
    
    # Optional: show notification
    showNotification("Heatmap input query reset. Please try a new search.", type = "message")
  })
  
  observeEvent(input$load_heatmap_gene_example_button, {
    updateTextInput(session, "heatmap_gene_input", value = "PAX8")
    updateTextInput(session, "heatmap_probeid_input", value = "")
    updateSelectInput(session, "heatmap_cohort_input", selected = "All")
  })
  
  observeEvent(input$load_heatmap_cpg_example_button, {
    updateTextInput(session, "heatmap_gene_input", value = "")
    updateTextInput(session, "heatmap_probeid_input", value = "cg11747594")
    updateSelectInput(session, "heatmap_cohort_input", selected = "All")
  })
  
  
  
  #EWAS-eQTM integration example button
  observeEvent(input$load_ewas_example_button, {
    updateTextInput(session, "ewas_probe_input", value = "cg00117599")
    updateSelectizeInput(session, "ewas_trait_input", selected = "atopy")
    updateSelectInput(session, "ewas_cohort_input", selected = "EVA-PR: Atopic asthma + controls")
  })
  
  observeEvent(input$load_gene_example_button, {
    updateTextAreaInput(session, "gene_cpg_list", value = "PDCD10, USP6")
    updateRadioButtons(session, "input_type", selected = "Gene")
  })
  
  observeEvent(input$load_cpg_example_button, {
    updateTextAreaInput(session, "gene_cpg_list", value = "cg05575921, cg00119778")
    updateRadioButtons(session, "input_type", selected = "CpG")
  })
  
  
  #Load Example button
  observeEvent(input$load_example_button, {
    updateTextInput(session, "gene_input", value = "PDCD10")
    updateSelectInput(session, "cohort_input", selected = "EVA-PR: Atopic asthma + controls")
  })
  
  #Load CpG Example button
  observeEvent(input$load_cpg_example_button, {
    updateTextInput(session, "probeid_input", value = "cg04983687")  # Set your example CpG probe ID here
    updateSelectInput(session, "cohort_input", selected = "EVA-PR: Atopic asthma + controls")
  })
  
  observeEvent(input$cohort_input, {
    if (input$cohort_input != "All") {
      updateSelectInput(session, "disease_input", selected = "All")
      shinyjs::disable("disease_input")
    } else {
      shinyjs::enable("disease_input")
    }
  })
  
  observeEvent(input$disease_input, {
    if (input$disease_input != "All") {
      updateSelectInput(session, "cohort_input", selected = "All")
      shinyjs::disable("cohort_input")
    } else {
      shinyjs::enable("cohort_input")
    }
  })
  
  result <- reactive({
    reset_state(FALSE)
    req(input$search_button)
    isolate({
      withProgress(message = "Searching Ensembl annotation data to render plot...", value = 0, {
        
        # Step 1: Load base data
        data <- combined_df_atopic_asthma_final
        incProgress(0.15, detail = "Loading eQTM dataset...")
        
        # Step 2: Apply cohort and disease filters
        if (input$cohort_input != "All" & input$disease_input != "All") {
          data <- data[data$cohort == input$cohort_input & data$disease == input$disease_input, ]
        } else if (input$cohort_input != "All") {
          data <- data[data$cohort == input$cohort_input, ]
        } else if (input$disease_input != "All") {
          data <- data[data$disease == input$disease_input, ]
        }
        incProgress(0.35, detail = "Filtering by cohort or disease...")
        
        # Step 3: Apply gene/probe filters
        gene_filter <- input$gene_input != ""
        probe_filter <- input$probeid_input != ""
        
        data <- data[
          (data$gene == input$gene_input & gene_filter) |
            (data$probeid == input$probeid_input & probe_filter),
        ]
        incProgress(0.25, detail = "Matching gene or probe ID...")
        
        # Step 4: Final processing step
        Sys.sleep(0.5)  # (optional) simulate final processing lag
        incProgress(0.25, detail = "Preparing plot output...")
        
        return(data)
      })
    })
  })
  
  
  
  observeEvent(input$example_button, {
    updateTextInput(session, "coloc_gene_input", value = "PDCD10")
    updateSelectInput(session, "coloc_cohort_input", selected = "EVA-PR: Atopic asthma + controls")
    updateSelectInput(session, "coloc_tissue_input", selected = "All")
    
    #Prepare example data for the textarea input
    example_data <- merged_df %>%
      dplyr::filter(gene == "PDCD10", cohort == "EVA-PR: Atopic asthma + controls") %>%
      dplyr::select(ProbeID, pvalue) %>%
      head(10)  # Get the first 10 rows, for example
    
    #Create the text to update the textAreaInput
    example_text <- paste(example_data$ProbeID, example_data$pvalue, sep = "\t", collapse = "\n")
    updateTextAreaInput(session, "ewas_data_text", value = example_text)
  })
  
  #Process the pasted EWAS data
  processed_ewas_data <- reactive({
    req(input$ewas_data_text)  # Ensure there's input
    #Split input text into lines
    lines <- strsplit(input$ewas_data_text, "\n")[[1]]
    #Process each line assuming the first value is cpgID and the second is pValue
    data <- do.call(rbind, lapply(lines, function(line) {
      strsplit(trimws(line), "\\s+")[[1]][1:2]  # Adjust based on expected format
    }))
    #Convert to data frame without column names in the input
    ewas_data <- as.data.frame(data, stringsAsFactors = FALSE)
    colnames(ewas_data) <- c("cpgID", "pValue")  # Explicitly set column names
    ewas_data$pValue <- as.numeric(ewas_data$pValue)  # Convert pValue to numeric
    return(ewas_data)
  })
  
  output$individuals_img <- renderImage({
    list(src = "www/people_icon.png", width = 80, height = 80)
  }, deleteFile = FALSE)
  
  output$eqtms_img <- renderImage({
    list(src = "www/eqtm_icon.png", width = 80, height = 80)
  }, deleteFile = FALSE)
  
  output$genes_img <- renderImage({
    list(src = "www/genes_icon.png", width = 80, height = 80)
  }, deleteFile = FALSE)
  
  output$tissues_img <- renderImage({
    list(src = "www/tissues_icon.png", width = 80, height = 80)
  }, deleteFile = FALSE)
  
  output$diseases_img <- renderImage({
    list(src = "www/diseases_icon.png", width = 80, height = 80)
  }, deleteFile = FALSE)
  
  output$upload_img <- renderImage({
    list(src = "www/upload_icon.png", width = 80, height = 80)
  }, deleteFile = FALSE)
  
  output$ewas_numberofpubs_img <- renderImage({
    list(src = "www/ewas_numberofpubs.png",
         width = "100%")
  }, deleteFile = FALSE)
  
  output$ewas_top10tissues_img <- renderImage({
    list(src = "www/ewas_top10tissues.png",
         width = "100%")
  }, deleteFile = FALSE)
  
  output$ewas_top10traits_img <- renderImage({
    list(src = "www/ewas_top10traits.png",
         width = "100%")
  }, deleteFile = FALSE)
  
  
  output$tissue_panels <- renderUI({
    library(ggplot2)
    
    # Get fully unique rows
    tissue_counts <- combined_df_atopic_asthma_final %>%
      dplyr::filter(!is.na(tissue)) %>%
      dplyr::distinct() %>%
      dplyr::count(tissue, name = "entry_count") %>%
      arrange(desc(entry_count))
    
    # Split Nasal Epithelium and others
    nasal_panel <- tissue_counts %>%
      dplyr::filter(tissue == "Nasal Epithelium")
    white_blood_panel <- tissue_counts %>%
      dplyr::filter(tissue == "White Blood Cell")
    
    
    other_panels <- tissue_counts %>%
      dplyr::filter(!tissue %in% c("Nasal Epithelium", "White Blood Cell"))
    
    
    pie_panel <- column(
      6,
      
      # ---- TITLE OUTSIDE ----
      div(
        h3("eQTM Pairs per Tissue", class = "tissue-panel-title"),
      ),
      
      # ---- PANEL BOX ----
      div(
        class = "tissue-panel",
        
        # Summary cards (Nasal + WBC)
        div(
          class = "tissue-summary-row",
          div(
            class = "tissue-summary-card",
            h5("Nasal Epithelium"),
            p(paste0(format(nasal_panel$entry_count, big.mark = ","), " pairs"))
          ),
          div(
            class = "tissue-summary-card",
            h5("White Blood Cell"),
            p(paste0(format(white_blood_panel$entry_count, big.mark = ","), " pairs"))
          )
        ),
        
        # Pie chart
        div(
          style = "display:flex; justify-content:center; align-items:center; width:100%;",
          plotlyOutput("tissue_pie_static", height = "540px")
        )
      )
    )
    
    
    
    ewas_graphs_panel <- column(
      6,
      div(
        h3("EWAS Atlas: Database Statistics Overview", class = "panel-title-align"),

        
        # Row 1: big publications card (full width)
        fluidRow(
          column(
            12,
            div(
              class = "ewas-card-big",
              imageOutput("ewas_numberofpubs_img", height = "100%")
            )
          )
        ),
        
        # Row 2: two cards side-by-side
        fluidRow(
          column(
            6,
            div(
              class = "ewas-card",
              imageOutput("ewas_top10tissues_img", height = "100%")
            )
          ),
          column(
            6,
            div(
              class = "ewas-card",
              imageOutput("ewas_top10traits_img", height = "100%")
            )
          )
        ),
        
        tags$p(
          "Summary visualizations adapted from the EWAS Atlas statistics page. https://ngdc.cncb.ac.cn/ewas/statistics",
          style = "font-size: 12px; color: #777; margin-top: 8px;"
        )
      )
    )
    
    
    
    
    
    
    # # Dedicated nasal epithelium panel
    # nasal_ui <- column(
    #   4,
    #   div(class = "tissue-panel",
    #       div(class = "panel-content",
    #           h4("Nasal Epithelium"),
    #           p(paste0(format(nasal_panel$entry_count, big.mark = ","), " pairs"))
    #       ),
    #       div(class = "description-overlay",
    #           p("Unique eQTM pairs - EVA-PR cohort."))
    #   )
    # )
    
    output$tissue_pie_static <- renderPlotly({
      tissue_counts <- combined_df_atopic_asthma_final %>%
        dplyr::filter(!is.na(tissue), !tissue %in% c("Nasal Epithelium", "White Blood Cell")) %>%
        dplyr::distinct() %>%
        dplyr::count(tissue, name = "entry_count") %>%
        dplyr::arrange(dplyr::desc(entry_count))
      
      # Safe color handling for >12 slices
      cols <- if (nrow(tissue_counts) <= 12) {
        RColorBrewer::brewer.pal(max(nrow(tissue_counts), 3), "Pastel1")
      } else {
        grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Pastel1"))(nrow(tissue_counts))
      }
      
      plotly::plot_ly(
        data = tissue_counts,
        labels = ~tissue,
        values = ~entry_count,
        type = "pie",
        # show counts on-slice
        text = ~format(entry_count, big.mark = ","),
        textinfo = "label+text",
        textposition = "outside",
        # show percent on hover
        hovertemplate = "<b>%{label}</b><br>%{percent:.1%} of pairs<extra></extra>",
        marker = list(
          line = list(color = "#F0F0F0", width = 2),
          colors = cols
        ),
        pull = ~ifelse(entry_count / sum(entry_count) < 0.05, 0.1, 0.03),
        domain = list(x = c(0, 1), y = c(0.2, 1))
      ) %>%
        plotly::layout(
          showlegend = FALSE,
          margin = list(t = 20, b = 80, l = 40, r = 40),
          uniformtext = list(minsize = 10, mode = "hide"),
          automargin = TRUE
        )
    })
    
    
    
    
    
    fluidRow(
      pie_panel,
      ewas_graphs_panel
    )
    
    
  })
  
  
  
  
  output$home_content <- renderUI({
    tagList(
      # --- Main summary panels ---
      fluidRow(
        column(4, div(class = "homepage-panel",
                      div(class = "panel-content", imageOutput("individuals_img", height = "120px")),
                      div(class = "description-overlay", p("Our database contains information from 5,266 individuals across 8 independent studies: Epigenetic Variation and Childhood Asthma in Puerto Ricans (EVA-PR, 2 unique studies), GTEx, FUSION, Framingham Heart Study, Colonomics, HELIX, and EvoImmunoPop.", class = "description-text")),
                      h3("5,266 individuals, 8 studies")
        )),
        column(4, div(class = "homepage-panel",
                      div(class = "panel-content", imageOutput("eqtms_img", height = "120px")),
                      div(class = "description-overlay", p("We provide detailed eQTM data, covering over 11 million associations between gene expression and DNA methylation.", class = "description-text")),
                      h3("11,609,019 eQTMs, 173,886 CpG probes")
        )),
        column(4, div(class = "homepage-panel",
                      div(class = "panel-content", imageOutput("genes_img", height = "120px")),
                      div(class = "description-overlay", p("Our repository includes data on over 20,000 unique genes involved in various biological processes.", class = "description-text")),
                      h3("20,231 unique genes")
        )),
        column(4, div(class = "homepage-panel",
                      div(class = "panel-content", imageOutput("tissues_img", height = "100px")),
                      div(class = "description-overlay", p("We analyze eQTMs across 14 collected tissue sources, helping provide insights into tissue-specific regulatory mechanisms.", class = "description-text")),
                      h3("11 unique tissues cataloged")
        )),
        column(4, div(class = "homepage-panel",
                      div(class = "panel-content", imageOutput("diseases_img", height = "100px")),
                      div(class = "description-overlay", p("In addition to healthy patients, our platform currently houses eQTM associations relevant to four major diseases: asthma, type II diabetes and related traits, cardiometabolic disease, and colon cancer.", class = "description-text")),
                      h3("4 diseases of interest")
        )),
        column(4, div(class = "homepage-panel",
                      div(class = "panel-content", imageOutput("upload_img", height = "100px")),
                      div(class = "description-overlay", p("Users can upload their own datasets for comparisons and analyses with existing eQTM data.", class = "description-text")),
                      h3("Upload data for comparison")
        ))
      ),
      # --- Tissue Entry Count Section ---
      fluidRow(
        column(12,
               div(style = "margin-top: 40px; text-align: center;"
                   #h3("eQTM Pairs per Tissue")
               )
        )
      ),
      uiOutput("tissue_panels")  # <- dynamically rendered small tissue panels
    )
  })
  
  # output$data_table_home <- DT::renderDataTable({
  #   selected_cohort <- input$cohort_filter
  #   data_home <- combined_df_atopic_asthma_final %>%
  #     dplyr::select(gene, probeid, cohort, pvalue, FDR, tissue) %>%
  #     distinct() %>%
  #     {if (selected_cohort != "All") dplyr::filter(., cohort == selected_cohort) else .} # Filter by selected cohort if not "All"
  #   
  #   DT::datatable(data_home, options = list(
  #     pageLength = 10,
  #     autoWidth = TRUE,
  #     searchHighlight = TRUE,
  #     order = list(list(4, 'asc')) #Default sort by pvalue column (index 4) in ascending order
  #   ))
  # })
  
  output$result_table <- DT::renderDataTable({
    if (reset_state()) return(NULL)
    
    res <- result()
    if (nrow(res) == 0) {
      return(data.frame("Message" = "No match found."))
    }
    
    # Drop all-NA rows but keep types as-is
    res <- res %>% dplyr::filter(rowSums(is.na(.)) != ncol(.))
    
    dt <- DT::datatable(
      res,
      options = list(
        scrollX   = TRUE,
        pageLength = 10,
        autoWidth = TRUE, 
        order = list(list(which(colnames(res) == "pvalue") - 1, "asc"))
      ),
      rownames = FALSE
    ) %>%
      # scientific format for pvalue/FDR/beta, but still sortable numerically
      DT::formatSignif(c("pvalue", "FDR", "beta"), digits = 3) %>%
      # integer formatting for positions/distances
      DT::formatRound(c("TSS", "pos", "dist"), digits = 0)
    
    dt
  })
  
  #reactive expressions to make manhattan plots available
  
  output$manhattanPlotAvailable <- reactive({
    !is.null(result()) && nrow(result()) > 0 && input$gene_input != ""
  })
  outputOptions(output, "manhattanPlotAvailable", suspendWhenHidden = FALSE)
  
  output$cpgManhattanPlotAvailable <- reactive({
    !is.null(result()) && nrow(result()) > 0 && input$probeid_input != ""
  })
  outputOptions(output, "cpgManhattanPlotAvailable", suspendWhenHidden = FALSE)
  
  output$transManhattanPlotAvailable <- reactive({
    !is.null(result()) && nrow(result()) > 0 && input$gene_input != ""
  })
  outputOptions(output, "transManhattanPlotAvailable", suspendWhenHidden = FALSE)
  
  output$transManhattanPlotByCpgAvailable <- reactive({
    !is.null(result()) && nrow(result()) > 0 && input$probeid_input != ""
  })
  outputOptions(output, "transManhattanPlotByCpgAvailable", suspendWhenHidden = FALSE)
  
  output$manhattan_plot <- renderPlotly({
    if (reset_state()) return(NULL)
    if (input$gene_input == "" && input$probeid_input != "") return(NULL)
    
    res <- result()
    res_cis <- res %>% dplyr::filter(eQTM_label == "cis")
    
    if (nrow(res_cis) == 0) {
      showModal(modalDialog(
        title = "No Results Found",
        "No cis-eQTM data points match your search query.",
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
      return(NULL)
    }
    
    # Safely wrap plotting block
    plot <- tryCatch({
      
      # Construct locus plot
      test <- locus(
        data   = res_cis,
        gene   = res_cis$gene,
        seqname = "gene.chr",
        ens_db  = "EnsDb.Hsapiens.v75",
        chrom   = "gene.chr",
        pos     = "pos",     # CpG site (bp)
        p       = "pvalue",
        labs    = "probeid",
        flank   = 1e6
      )
      
      # --- FORCE ALL POINTS to pch 21 and slateblue ---
      test$data$pch <- 21                   # pch 21 = filled circle
      test$data$col <- "slateblue"          # outline color
      test$data$bg  <- "slateblue"          # fill color
      
      
      plot <- locus_plotly(test)
      
      # Recolor points
      plot$x$data <- lapply(plot$x$data, function(trace) {
        if (trace$mode == "markers") trace$marker$color <- "slateblue"
        trace
      })
      
      # ---- NEW: place the annotation at the gene’s TSS ----
      # res_cis$TSS is already in bp
      tss_bp <- as.numeric(res_cis$TSS[1])
      tss_mb <- tss_bp / 1e6   # convert to Mb for locus_plotly x-axis
      
      plot <- plot %>% layout(
        showlegend = FALSE,
        title = paste0("cis-CpG and gene genomic position plot: ", res_cis$gene[1]),
        annotations = list(
          list(
            text = paste0("<b>", res_cis$gene[1], "</b>"),
            
            # Position annotation at the TSS
            x    = tss_mb,
            xref = "x",
            
            # Keep visible near bottom of plot
            y    = 0.42,
            yref = "paper",
            
            showarrow = FALSE,
            arrowcolor = "rgba(100,100,100,0.3)",
            arrowwidth = 1.5,
            ax = 0,
            ay = -40,
            
            font = list(size = 8, color = "black"),
            bgcolor = "lightyellow",#"rgba(255,255,180,0.5)",
            bordercolor = "rgba(0,0,0,0.2)",
            borderwidth = 0.8
          )
        )
      )
      
      plot
      
    }, error = function(e) {
      showModal(modalDialog(
        title = "Plot Generation Error",
        "An issue occurred while generating the plot.",
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
      NULL
    })
    
    return(plot)
  })
  
  
  #Populate filter dropdowns on app load
  observe({
    updateSelectInput(
      session, "filter_cohort", 
      choices = c("All", unique(combined_df_atopic_asthma_final$cohort)), 
      selected = "All"
    )
    updateSelectInput(
      session, "filter_tissue", 
      choices = c("All", unique(combined_df_atopic_asthma_final$tissue)), 
      selected = "All"
    )
  })
  
  #Reactive function to process the input and apply filters
  gene_cpg_results <- reactive({
    req(input$search_gene_cpg)
    
    #Get input type and list
    input_type <- input$input_type
    search_list <- strsplit(input$gene_cpg_list, "[,\\n\\r]+")[[1]]
    search_list <- trimws(search_list)
    
    #Filter data based on input type
    results <- if (input_type == "Gene") {
      combined_df_atopic_asthma_final %>%
        dplyr::filter(gene %in% search_list) %>%
        dplyr::select(gene, probeid, cohort, tissue, pvalue, beta, FDR) %>%
        distinct()
    } else {
      combined_df_atopic_asthma_final %>%
        dplyr::filter(probeid %in% search_list) %>%
        dplyr::select(probeid, gene, cohort, tissue, pvalue, beta, FDR) %>%
        distinct()
    }
    
    #Apply cohort filter
    if (input$filter_cohort != "All") {
      results <- results %>% 
        dplyr::filter(cohort == input$filter_cohort)
    }
    
    #Apply tissue filter
    if (input$filter_tissue != "All") {
      results <- results %>% 
        dplyr::filter(tissue == input$filter_tissue)
    }
    
    return(results)
  })
  
  output$gene_cpg_table <- DT::renderDataTable({
    res <- gene_cpg_results()
    if (is.null(res) || nrow(res) == 0) return(NULL)
    
    # Make sure there is a 'pvalue' column
    if ("pvalue" %in% names(res)) {
      # Ensure pvalue is numeric for proper sorting
      if (!is.numeric(res$pvalue)) {
        suppressWarnings({
          res$pvalue <- as.numeric(res$pvalue)
        })
      }
      
      # Build options with default sort by pvalue (0-based index)
      pval_col_idx <- which(names(res) == "pvalue") - 1
      opts <- list(
        pageLength = 10,
        autoWidth  = TRUE,
        order      = list(list(pval_col_idx, "asc"))
      )
    } else {
      # Fallback: no pvalue column, no default ordering
      opts <- list(
        pageLength = 10,
        autoWidth  = TRUE
      )
    }
    
    DT::datatable(
      res,
      options  = opts,
      rownames = FALSE
    )
  })
  
  
  #Download filtered results
  output$download_gene_cpg_results <- downloadHandler(
    filename = function() { paste0("gene_cpg_results_", Sys.Date(), ".csv") },
    content = function(file) {
      write.csv(gene_cpg_results(), file, row.names = FALSE)
    }
  )
  
  
  output$cpg_manhattan_plot <- renderPlotly({
    if (reset_state()) return(NULL)  # Do not render if reset flag is active
    if (input$probeid_input == "" && input$gene_input != "") {
      return(NULL)  # Return NULL if only a gene is entered
    }
    
    res <- result()
    cpg_res <- res %>%
      dplyr::filter(probeid == input$probeid_input & eQTM_label == "cis")
    
    # If no cis CpG data is found, show a modal
    if (nrow(cpg_res) == 0) {
      showModal(modalDialog(
        title = "No Results Found",
        "No cis CpG data points match your search query. Please refresh and try another CpG site.",
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
      return(NULL)
    }
    
    # Clean chr prefix for compatibility with EnsDb.Hsapiens.v75
    cpg_res$meth.chr <- gsub("^chr", "", cpg_res$meth.chr)
    
    # Check if gene(s) have valid start/end coordinates
    gene_info <- AnnotationDbi::select(
      EnsDb.Hsapiens.v75,
      keys    = cpg_res$gene,
      keytype = "GENENAME",
      columns = c("GENEID", "SEQNAME", "GENESEQSTART", "GENESEQEND")
    ) %>% dplyr::distinct()
    
    gene_info$SEQNAME <- gsub("^chr", "", gene_info$SEQNAME)
    
    if (nrow(gene_info) == 0 ||
        all(is.na(gene_info$GENESEQSTART)) ||
        all(is.na(gene_info$GENESEQEND))) {
      showModal(modalDialog(
        title = "Gene Cannot Be Plotted",
        paste0(
          "None of the genes associated with '", input$probeid_input,
          "' have sufficient structure (start/end coordinates) to render a plot."
        ),
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
      return(NULL)
    }
    
    # Try to create the plot
    plot <- tryCatch({
      # 1) POINT POSITIONS: use gene TSS as x
      test <- locus(
        data    = cpg_res,
        gene    = cpg_res$gene,
        seqname = "meth.chr",
        ens_db  = "EnsDb.Hsapiens.v75",
        chrom   = "meth.chr",
        pos     = "TSS",      # <- genes plotted at TSS
        p       = "pvalue",
        labs    = "gene",
        flank   = 2e6
      )
      
      test$data$pch <- 21                   # pch 21 = filled circle
      test$data$col <- "slateblue"          # outline color
      test$data$bg  <- "slateblue"          # fill color
      
      
      plot <- locus_plotly(test)
      
      # Color all marker points uniformly
      plot$x$data <- lapply(plot$x$data, function(trace) {
        if (trace$mode == "markers") {
          trace$marker$color <- "slateblue"
        }
        trace
      })
      
      # 2) CENTERING + LABEL: center window on CpG (pos), not TSS
      label_row <- cpg_res %>%
        dplyr::filter(probeid == input$probeid_input) %>%
        dplyr::slice(1)
      
      # CpG position in bp
      cpg_bp <- as.numeric(label_row$pos[1])
      window <- 1.5e6   # +/- 1 Mb around CpG, to match flank
      
      cpg_mb   <- cpg_bp / 1e6
      x_min_mb <- (cpg_bp - window) / 1e6
      x_max_mb <- (cpg_bp + window) / 1e6
      
      # --- Add ONE CpG point (red) on the x-axis at x = cpg_mb ---
      
      marker_idx <- which(
        vapply(
          plot$x$data,
          function(tr) !is.null(tr$mode) && tr$mode == "markers",
          logical(1)
        )
      )[1]
      
      if (!is.na(marker_idx)) {
        tr <- plot$x$data[[marker_idx]]
        
        # Make sure x/y/text are vectors
        tr$x <- as.numeric(tr$x)
        tr$y <- as.numeric(tr$y)
        n_old <- length(tr$x)
        
        if (is.null(tr$text)) {
          tr$text <- rep("", n_old)
        }
        
        # Append CpG point at y = 0 (x-axis)
        tr$x <- c(tr$x, cpg_mb)
        tr$y <- c(tr$y, 0)
        tr$text <- c(
          tr$text,
          paste0(
            "CpG: ", input$probeid_input,
            "<br>Chr: ", label_row$meth.chr,
            "<br>Pos: ", format(cpg_bp, big.mark = ",")
          )
        )
        
        # Marker list
        if (is.null(tr$marker)) tr$marker <- list()
        
        ## Sizes
        if (is.null(tr$marker$size)) {
          tr$marker$size <- rep(5, n_old)
        } else if (length(tr$marker$size) == 1L) {
          tr$marker$size <- rep(tr$marker$size, n_old)
        } else {
          tr$marker$size <- rep(tr$marker$size, length.out = n_old)
        }
        tr$marker$size <- c(tr$marker$size, 8)  # CpG point a bit larger
        
        ## Colors: keep existing, last one red
        if (is.null(tr$marker$color)) {
          tr$marker$color <- rep("slateblue", n_old)
        } else if (length(tr$marker$color) == 1L) {
          tr$marker$color <- rep(tr$marker$color, n_old)
        } else {
          tr$marker$color <- rep(tr$marker$color, length.out = n_old)
        }
        tr$marker$color <- c(tr$marker$color, "red")
        
        plot$x$data[[marker_idx]] <- tr
      }
      
      
      plot <- plot %>% layout(
        xaxis = list(range = c(x_min_mb, x_max_mb)),  # <- CpG-centered axis
        showlegend = FALSE,
        title = paste0("cis-CpG and gene genomic position plot: ", input$probeid_input),
        annotations = list(
          list(
            text = paste0(
              "<b>",
              paste(unique(cpg_res$probeid), collapse = ", "),
              "</b>"
            ),
            x    = cpg_mb,   # label at CpG position
            xref = "x",
            y    = 0.44,
            yref = "paper",
            showarrow = FALSE,
            arrowcolor = "rgba(100, 100, 100, 0.3)",
            arrowwidth = 1.5,
            ax = 0,
            ay = -40,
            font = list(size = 8, color = "black"),
            bgcolor = "lightyellow",
            bordercolor = "rgba(0, 0, 0, 0.2)",
            borderwidth = 1
          )
        )
      )
      
      plot
    }, error = function(e) {
      showModal(modalDialog(
        title = "Plot Generation Error",
        "An issue occurred while generating the cis CpG plot. Please try another CpG site.",
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
      NULL
    })
    
    plot
  })
  
  
  
  
  observe({
    if (input$heatmap_gene_input != "" || input$heatmap_probeid_input != "") {
      heatmap_reset_state(FALSE)
    }
  })
  
  heatmap_data_processed <- eventReactive(input$search_heatmap_button, {
    withProgress(message = "Processing data for heatmap...", value = 0.1, {
      
      heatmap_data <- combined_df_atopic_asthma_final %>%
        dplyr::filter((gene == input$heatmap_gene_input & input$heatmap_gene_input != "") | 
                        (probeid == input$heatmap_probeid_input & input$heatmap_probeid_input != ""))
      
      if (input$heatmap_cohort_input != "All") {
        heatmap_data <- heatmap_data %>% dplyr::filter(cohort == input$heatmap_cohort_input)
      }
      
      incProgress(0.5)
      
      heatmap_data <- heatmap_data %>%
        mutate(transformed_pvalue = -log10(pvalue + 1e-10))
      
      incProgress(1)
      
      return(heatmap_data)
    })
  })
  
  
  # Normalize once
  normalized_eqtm <- reactive({
    combined_df_atopic_asthma_final %>%
      mutate(
        gene   = trimws(gene),
        probeid= trimws(probeid),
        tissue = trimws(tissue),
        pvalue = suppressWarnings(as.numeric(pvalue))
      )
  })
  
  heatmap_data_processed <- eventReactive(input$search_heatmap_button, {
    req(nzchar(input$heatmap_gene_input) || nzchar(input$heatmap_probeid_input))
    
    dat <- normalized_eqtm()
    
    if (input$heatmap_cohort_input != "All") {
      dat <- dat %>% dplyr::filter(cohort == input$heatmap_cohort_input)
    }
    
    g_in <- trimws(input$heatmap_gene_input)
    p_in <- trimws(input$heatmap_probeid_input)
    
    if (nzchar(g_in) && !nzchar(p_in)) {
      dat <- dat %>% dplyr::filter(toupper(gene) == toupper(g_in))
    } else if (!nzchar(g_in) && nzchar(p_in)) {
      dat <- dat %>% dplyr::filter(toupper(probeid) == toupper(p_in))
    } else if (nzchar(g_in) && nzchar(p_in)) {
      dat <- dat %>% dplyr::filter(toupper(gene) == toupper(g_in) | toupper(probeid) == toupper(p_in))
    }
    
    validate(need(nrow(dat) > 0, "No rows after filtering for heatmap."))
    
    dat <- dat %>%
      mutate(transformed_pvalue = -log10(pvalue + 1e-10))
    
    # --- Deduplicate before pivot ---
    if (nzchar(g_in) && !nzchar(p_in)) {
      # tissue x probeid
      dat2 <- dat %>%
        dplyr::filter(!is.na(tissue), !is.na(probeid)) %>%
        dplyr::group_by(tissue, probeid) %>%
        dplyr::summarise(transformed_pvalue = max(transformed_pvalue, na.rm = TRUE), .groups = "drop")
      
      mat <- dat2 %>%
        tidyr::pivot_wider(
          names_from  = probeid,
          values_from = transformed_pvalue,
          values_fill = 0,
          values_fn   = list(transformed_pvalue = max)
        )
    } else {
      # tissue x gene
      dat2 <- dat %>%
        dplyr::filter(!is.na(tissue), !is.na(gene)) %>%
        dplyr::group_by(tissue, gene) %>%
        dplyr::summarise(transformed_pvalue = max(transformed_pvalue, na.rm = TRUE), .groups = "drop")
      
      mat <- dat2 %>%
        tidyr::pivot_wider(
          names_from  = gene,
          values_from = transformed_pvalue,
          values_fill = 0,
          values_fn   = list(transformed_pvalue = max)
        )
    }
    
    # Long form for Plotly heatmap
    mat_long <- mat %>%
      tidyr::pivot_longer(-tissue, names_to = "Identifier", values_to = "LogPValue")
    
    list(mat = mat, mat_long = mat_long)
  })
  
  output$gene_heatmap <- renderPlotly({
    req(!heatmap_reset_state())
    hdat <- heatmap_data_processed()
    validate(need(nrow(hdat$mat_long) > 0, "No data to draw heatmap."))
    
    plot_ly(
      data = hdat$mat_long,
      x = ~Identifier,
      y = ~tissue,
      z = ~LogPValue,
      type = "heatmap",
      colors = colorRampPalette(RColorBrewer::brewer.pal(3, "Reds"))(256),
      hovertemplate =
        paste(
          "<b>Gene/Probe:</b> %{x}<br>",
          "<b>Tissue:</b> %{y}<br>",
          "<b>-log10(p-value):</b> %{z}<extra></extra>"
        )
    ) %>%
      layout(
        yaxis = list(title = "Tissue", autorange = "reversed"),
        xaxis = list(title = "Gene/Probe ID"),
        coloraxis = list(colorbar = list(title = "-log10(p-value)")),
        autosize = TRUE, width = 2000, height = 350
      )
  })
  
  # Re-enable after typing
  observe({
    if (nzchar(input$heatmap_gene_input) || nzchar(input$heatmap_probeid_input)) {
      heatmap_reset_state(FALSE)
    }
  })
  
  
  output$download_button_ui <- renderUI({
    if (reset_state()) return(NULL)  # Do not render if reset flag is active
    res <- result()
    if (nrow(res) > 0) {
      downloadButton("download_data", "Download Data", class = "download-button")
    }
  })
  
  output$trans_manhattan_plot <- renderPlotly({
    if (reset_state()) return(NULL)  # Do not render if reset flag is active
    if (input$gene_input == "" && input$probeid_input != "") {
      return(NULL)  # Return NULL to not render this plot
    }
    
    res <- result()
    res_trans <- res %>% dplyr::filter(eQTM_label == "trans")
    
    # If no trans-eQTM data is found, show a modal pop-up instead of a notification
    if (nrow(res_trans) == 0) {
      showModal(modalDialog(
        title = "No Results Found",
        "No trans-eQTM data points match your search query. Please refresh the page and try another gene or CpG site.",
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
      return(NULL)
    }
    
    # Try to create the plot and catch any errors
    plot <- tryCatch({
      # Ensure that the necessary columns are in the correct format
      res_trans$SNP <- res_trans$probeid
      res_trans$CHR <- as.numeric(gsub("chr", "", res_trans$meth.chr))
      res_trans$BP <- as.numeric(res_trans$pos)
      res_trans$P <- as.numeric(res_trans$pvalue)
      
      # Prepare data for plotly
      res_trans <- res_trans %>% dplyr::select(SNP, CHR, BP, P, cohort)
      
      chromosome_lengths <- res_trans %>%
        group_by(CHR) %>%
        summarize(max_bp = max(BP)) %>%
        arrange(CHR)
      
      chromosome_offsets <- c(0, cumsum(chromosome_lengths$max_bp))[1:nrow(chromosome_lengths)]
      names(chromosome_offsets) <- chromosome_lengths$CHR
      
      res_trans <- res_trans %>%
        mutate(cumulative_bp = BP + chromosome_offsets[as.character(CHR)])
      
      # Range for x-axis so there’s no extra white space inside the box
      x_min <- min(res_trans$cumulative_bp, na.rm = TRUE)
      x_max <- max(res_trans$cumulative_bp, na.rm = TRUE)
      x_pad <- (x_max - x_min) * 0.01   # 1% padding on each side
      
      
      # Add chromosome labels at the center of each chromosome
      chr_labels <- res_trans %>%
        group_by(CHR) %>%
        summarize(CenterBP = mean(cumulative_bp, na.rm = TRUE))
      
      # Alternate colors for chromosomes
      res_trans <- res_trans %>%
        mutate(color = ifelse(CHR %% 2 == 0, 'firebrick', 'dodgerblue'))
      
      # Create the interactive Manhattan plot with plotly
      p <- plot_ly(
        data = res_trans,
        x = ~cumulative_bp,
        y = ~-log10(P),
        type = 'scatter',
        mode = 'markers',
        marker = list(
          color = ~color,
          line = list(
            color = "black",
            width = 0.3   # try 0.3–0.6 for perfect thickness
          )
        ),
        text = ~paste('SNP:', SNP, '<br>Chr:', CHR, '<br>Pos:', BP, '<br>P:', P),
        hoverinfo = 'text'
      ) %>%
        layout(
          title = 'trans-CpG and gene genomic position plot',
          
          # ---- X AXIS ----
          xaxis = list(
            title = 'Chromosome Position',
            tickvals = chr_labels$CenterBP,
            ticktext = chr_labels$CHR,
            range     = c(x_min - x_pad, x_max + x_pad),  # <-- controlled whitespace
            tickangle = 0,                # <- never auto-tilt the labels
            automargin = TRUE,
            tickfont = list(
              family = "Arial",   # or "Helvetica", "Sans-Serif"
              size = 11,
              color = "black"
            ),
            showgrid = FALSE,      # ❌ remove gridlines
            zeroline = FALSE,
            showline = TRUE,       # ✔ draw main axis line
            linecolor = "black",
            linewidth = 1.5,
            mirror = TRUE          # ✔ draws right-hand axis (box effect)
          ),
          
          # ---- Y AXIS ----
          yaxis = list(
            title = '-log10(P)',
            
            showgrid = FALSE,      # ❌ remove gridlines
            zeroline = FALSE,
            showline = TRUE,       # ✔ draw main axis line
            linecolor = "black",
            linewidth = 1.5,
            mirror = TRUE          # ✔ draws top axis (box effect)
          ),
          
          legend = list(title = list(text = 'Cohort')),
          
          shapes = list(
            list(
              type = 'line',
              x0 = 0,
              x1 = 1,
              xref = 'paper',
              y0 = -log10(5e-8),
              y1 = -log10(5e-8),
              line = list(dash = 'dash', color = 'red')
            )
          ),
          
          # Optional: ensure no weird background shading
          paper_bgcolor = "white",
          plot_bgcolor  = "white"
        )
      
      
      return(p)
    }, error = function(e) {
      showModal(modalDialog(
        title = "Plot Generation Error",
        "There was an issue generating the trans Manhattan plot. Please try subsetting by cohort, or try a different gene or CpG site.",
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
      return(NULL)
    })
    
    return(plot)
  })
  
  
  output$trans_manhattan_plot_by_cpg <- renderPlotly({
    if (reset_state()) return(NULL)  # Do not render if reset flag is active
    req(input$probeid_input)  # Ensure a CpG site is provided
    res <- result()
    
    # Filter data for trans Manhattan plot when searched by CpG site
    res_trans <- res %>% dplyr::filter(probeid == input$probeid_input & eQTM_label == "trans")
    
    # If no trans CpG data is found, show a modal pop-up
    if (nrow(res_trans) == 0) {
      showModal(modalDialog(
        title = "No Results Found",
        "No trans CpG data points match your search query. Please refresh the page and try another CpG site.",
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
      return(NULL)
    }
    
    # Try to create the plot and catch any errors
    plot <- tryCatch({
      res_trans$SNP <- res_trans$probeid
      res_trans$CHR <- as.numeric(gsub("chr", "", res_trans$gene.chr))
      res_trans$BP <- as.numeric(res_trans$TSS)
      res_trans$P <- as.numeric(res_trans$pvalue)
      
      # Prepare data for plotly
      res_trans <- res_trans %>% dplyr::select(SNP, CHR, BP, P, gene, gene.chr, TSS, cohort)
      
      chromosome_lengths <- res_trans %>%
        group_by(CHR) %>%
        summarize(max_bp = max(BP)) %>%
        arrange(CHR)
      
      chromosome_offsets <- c(0, cumsum(chromosome_lengths$max_bp))[1:nrow(chromosome_lengths)]
      names(chromosome_offsets) <- chromosome_lengths$CHR
      
      res_trans <- res_trans %>%
        mutate(cumulative_bp = BP + chromosome_offsets[as.character(CHR)])
      
      # Range for x-axis so there’s no extra white space inside the box
      x_min <- min(res_trans$cumulative_bp, na.rm = TRUE)
      x_max <- max(res_trans$cumulative_bp, na.rm = TRUE)
      x_pad <- (x_max - x_min) * 0.01   # 1% padding on each side
      
      # Add chromosome labels at the center of each chromosome
      chr_labels <- res_trans %>%
        group_by(CHR) %>%
        summarize(CenterBP = mean(cumulative_bp, na.rm = TRUE))
      
      # Alternate colors for chromosomes
      res_trans <- res_trans %>%
        mutate(color = ifelse(CHR %% 2 == 0, 'firebrick', 'dodgerblue'))
      
      p <- plot_ly(
        data = res_trans,
        x = ~cumulative_bp,
        y = ~-log10(P),
        type = 'scatter',
        mode = 'markers',
        marker = list(
          color = ~color,
          line = list(
            color = "black",
            width = 0.3   # try 0.3–0.6 for perfect thickness
          )
        ),
        text = ~paste('Gene:', gene, '<br>Chr:', gene.chr, '<br>TSS:', TSS, '<br>P:', P),
        hoverinfo = 'text'
      ) %>%
        layout(
          title = 'trans-CpG and gene genomic position plot',
          
          # ---------------- X AXIS ----------------
          xaxis = list(
            title = 'Chromosome Position',
            tickvals = chr_labels$CenterBP,
            ticktext = chr_labels$CHR,
            range     = c(x_min - x_pad, x_max + x_pad),  # <-- controlled whitespace
            tickangle = 0,                # <- never auto-tilt the labels
            automargin = TRUE,
            tickfont = list(
              family = "Arial",   # or "Helvetica", "Sans-Serif"
              size = 11,
              color = "black"
            ),
            showgrid   = FALSE,     # remove vertical gridlines
            zeroline   = FALSE,
            showline   = TRUE,      # draw solid x-axis
            linecolor  = "black",
            linewidth  = 1.5,
            mirror     = TRUE       # draw top border → box effect
          ),
          
          # ---------------- Y AXIS ----------------
          yaxis = list(
            title = '-log10(P)',
            
            showgrid   = FALSE,     # remove horizontal gridlines
            zeroline   = FALSE,
            showline   = TRUE,      # draw solid y-axis
            linecolor  = "black",
            linewidth  = 1.5,
            mirror     = TRUE       # draw right border → box effect
          ),
          
          legend = list(title = list(text = 'Cohort')),
          
          shapes = list(
            list(
              type  = 'line',
              x0    = 0,
              x1    = 1,
              xref  = 'paper',
              y0    = -log10(5e-8),
              y1    = -log10(5e-8),
              line  = list(dash = 'dash', color = 'red')
            )
          ),
          
          # Make background clean white like the GWAS figure
          paper_bgcolor = "white",
          plot_bgcolor  = "white"
        )
      
      
      return(p)
    }, error = function(e) {
      showModal(modalDialog(
        title = "Plot Generation Error",
        "There was an issue generating the trans CpG Manhattan plot. Please try subsetting by cohort, or try a different CpG site.",
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
      return(NULL)
    })
    
    return(plot)
  })
  
  
  output$download_data <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      data_to_download <- result() %>%
        dplyr::filter(rowSums(is.na(.)) != ncol(.))  # Filter out rows where all columns are NA
      write.csv(data_to_download, file, row.names = FALSE)
    }
  )
  
  observeEvent(input$reset_ewas_inputs, {
    updateTextInput(session, "ewas_probe_input", value = "")
    updateSelectizeInput(session, "ewas_trait_input", selected = "")
    updateSelectInput(session, "ewas_cohort_input", selected = "Select a cohort")
    
    ewas_reset_state(TRUE)  # Reset flag for EWAS tab
  })
  
  
  # Ensure you have the function defined correctly:
  get_local_ewas_data <- function(probe_id = NULL, gene_symbol = NULL) {
    if (is.null(probe_id) && is.null(gene_symbol)) {
      return(data.frame(Message = "Please provide a CpG ID or a Gene Symbol."))
    }
    
    if (!is.null(gene_symbol) && gene_symbol != "") {
      res <- ewas_atlas %>%
        dplyr::filter(gene %in% gene_symbol)
    } else if (!is.null(probe_id) && probe_id != "") {
      res <- ewas_atlas %>%
        dplyr::filter(probe_ID %in% probe_id)
    } else {
      return(data.frame(Message = "Invalid input provided."))
    }
    
    if (nrow(res) == 0) {
      return(data.frame(Message = "No data found in EWAS Atlas."))
    } else {
      return(res)
    }
  }
  
  observeEvent(input$search_ewas, {
    ewas_reset_state(FALSE)  # <--- reset the flag so results render
    if (input$ewas_cohort_input == "Select a cohort") {
      showModal(modalDialog(
        title = "Cohort Selection Required",
        "Please select an eQTM cohort to proceed with the search.",
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
      return()
    }
    
    if (input$ewas_probe_input == "" && input$ewas_gene_input == "" && input$ewas_trait_input == "") {
      showModal(modalDialog(
        title = "Invalid Input",
        "Please enter either a CpG ID, a Gene Symbol, or a Trait keyword to search.",
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
      return()
    }
    
    ewas_tables <- list()
    
    # ----------------- 1. PROBE SEARCH -----------------
    if (nzchar(input$ewas_probe_input)) {
      associations_df <- get_local_ewas_data(probe_id = input$ewas_probe_input) %>%
        dplyr::filter(!is.na(p_value) & p_value != "") %>%
        dplyr::select(probe_ID, trait, p_value, study_ID) %>%
        dplyr::rename(EWAS_pvalue = p_value) %>%
        distinct()
      
      matched_eqtm_genes <- combined_df_atopic_asthma_final %>%
        dplyr::filter(probeid == input$ewas_probe_input,
                      cohort == input$ewas_cohort_input) %>%
        dplyr::select(probe_ID = probeid, eQTM_gene = gene, eQTM_pvalue = pvalue, eQTM_FDR = FDR) %>%
        distinct()
      
      associations_expanded <- associations_df %>%
        left_join(matched_eqtm_genes, by = "probe_ID")
      # Sort by eQTM_pvalue if available
      if ("eQTM_pvalue" %in% colnames(associations_expanded)) {
        associations_expanded <- associations_expanded %>%
          arrange(eQTM_pvalue)
      }
      
      
      if (nrow(associations_expanded) == 0) {
        associations_expanded <- data.frame(Message = "No matching EWAS or eQTM data found.")
      }
      
      ewas_tables[["EWAS Associations"]] <- associations_expanded
      
      if (nrow(matched_eqtm_genes) == 0) {
        matched_eqtm_genes <- data.frame(Message = "No matching eQTM genes found for this probe.")
      }
      
      # ewas_tables[["Matched eQTM Genes"]] <- matched_eqtm_genes
    }
    
    
    # ----------------- 2. TRAIT SEARCH -----------------
    if (nzchar(input$ewas_trait_input)) {
      library(stringi)
      
      matched_trait_df <- ewas_atlas %>%
        dplyr::filter(
          stri_detect_regex(
            trait,
            input$ewas_trait_input,
            opts_regex = list(case_insensitive = TRUE)
          )
        ) %>%
        distinct(probe_ID, trait, p_value, study_ID) %>%
        dplyr::rename(EWAS_pvalue = p_value) %>%
        dplyr::filter(!is.na(EWAS_pvalue), EWAS_pvalue != "")
      
      
      #f probe is specified, filter to just that one
      if (nzchar(input$ewas_probe_input)) {
        matched_trait_df <- matched_trait_df %>%
          dplyr::filter(probe_ID == input$ewas_probe_input)
      }
      
      #Join to eQTM dataset using probe ID
      matched_eqtm_df <- matched_trait_df %>%
        dplyr::inner_join(combined_df_atopic_asthma_final, by = c("probe_ID" = "probeid")) %>%
        dplyr::filter(cohort == input$ewas_cohort_input) %>%
        dplyr::select(probe_ID, trait, EWAS_pvalue, study_ID,
                      eQTM_gene = gene, cohort, tissue, eQTM_pvalue = pvalue, eQTM_FDR = FDR) %>%
        distinct()
      
      # Sort by eQTM_pvalue when present
      if ("eQTM_pvalue" %in% colnames(matched_eqtm_df)) {
        matched_eqtm_df <- matched_eqtm_df %>% arrange(eQTM_pvalue)
      }
      
      
      if (nrow(matched_eqtm_df) == 0) {
        matched_eqtm_df <- data.frame(Message = "No matching probes in eQTM for this trait.")
      }
      
      ewas_tables[["EWAS Associations"]] <- matched_eqtm_df
    }
    
    
    # UI rendering logic (unchanged)
    output$ewas_headers <- renderUI({
      if (ewas_reset_state()) return(NULL)
      bslib::navs_tab(
        !!!lapply(names(ewas_tables), function(table_name) {
          table_id <- paste0("ewas_table_", gsub("\\s+", "_", table_name))
          download_id <- paste0("download_", gsub("\\s+", "_", table_name))
          
          tabPanel(
            title = table_name,
            downloadButton(download_id, paste0("Download ", table_name), class = "download-button"),
            div(style = "overflow-x: auto;", DT::dataTableOutput(table_id))
          )
        })
      )
    })
    
    lapply(names(ewas_tables), function(table_name) {
      table_id <- paste0("ewas_table_", gsub("\\s+", "_", table_name))
      download_id <- paste0("download_", gsub("\\s+", "_", table_name))
      
      output[[table_id]] <- DT::renderDataTable({
        if (ewas_reset_state()) return(NULL)
        DT::datatable(ewas_tables[[table_name]],
                      options = list(pageLength = 20, autoWidth = TRUE),
                      rownames = FALSE)
      })
      
      output[[download_id]] <- downloadHandler(
        filename = function() paste0(gsub("\\s+", "_", table_name), "_", Sys.Date(), ".csv"),
        content = function(file) write.csv(ewas_tables[[table_name]], file, row.names = FALSE)
      )
    })
  })
}

# Run the application 
shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))
