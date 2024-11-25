# Load required libraries
library(shiny)
library(ggplot2)
library(dplyr)
library(readr)
library(httr)
library(jsonlite)
library(utils)

# Define a mapping of column names to prettier terms
pretty_names <- c(
  "ambig_basecount" = "Ambiguous bases within barcode",
  "ambig_full_basecount" = "Ambiguous bases in the full sequence",
  "nuc_basecount" = "Sequence length within barcode",
  "nuc_full_basecount" = "Overall sequence length",
  "stop_codons" = "Number of stop codons in the barcode translation"
)

# Function to get TSV files from GitHub
get_github_tsv_files <- function(organization) {
  url <- paste0("https://api.github.com/repos/naturalis/barcode_validator/contents/data?ref=", organization)
  response <- GET(url)
  content <- content(response, "text")
  files <- fromJSON(content)
  tsv_files <- files$name[grepl("\\.tsv$", files$name)]
  return(tsv_files)
}

# UI
ui <- fluidPage(
  titlePanel("BGE Genome Skimming Results"),

  sidebarLayout(
    sidebarPanel(
      selectInput("organization", "Select Organization",
                  choices = c("naturalis", "nhm", "unifi"),
                  selected = "naturalis"),
      selectizeInput("files", "Select TSV File(s)",
                     choices = NULL,
                     multiple = TRUE),
      selectInput("field", "Select Field for Histogram",
                  choices = setNames(names(pretty_names), pretty_names)),
      actionButton("analyze", "Analyze Selected Files")
    ),

    mainPanel(
      plotOutput("histogram"),

      h3("Summary Statistics"),
      htmlOutput("summary_stats")
    )
  )
)

# Server
server <- function(input, output, session) {

  # Update file list when organization changes
  observe({
    tsv_files <- get_github_tsv_files(input$organization)
    updateSelectizeInput(session, "files", choices = tsv_files)
  })

  data <- eventReactive(input$analyze, {
    req(input$files)

    all_data <- lapply(input$files, function(file) {
      encoded_file <- URLencode(file, reserved = TRUE)
      url <- paste0("https://raw.githubusercontent.com/naturalis/barcode_validator/",
                   input$organization, "/data/", encoded_file)

      tryCatch({
        df <- read_tsv(url, na = c("", "NA", "None"))

        # Convert relevant columns to numeric, replacing 'None' with NA
        numeric_cols <- names(pretty_names)
        df <- df %>%
          mutate(across(all_of(numeric_cols), ~as.numeric(as.character(.))))

        return(df)
      }, error = function(e) {
        warning(paste("Error reading file:", file, "-", e$message))
        return(NULL)
      })
    })

    # Remove NULL entries from failed reads
    all_data <- all_data[!sapply(all_data, is.null)]

    if (length(all_data) == 0) {
      return(NULL)
    }

    bind_rows(all_data)
  })

  output$histogram <- renderPlot({
    req(data())
    if (nrow(data()) == 0) return(NULL)

    ggplot(data(), aes_string(x = input$field)) +
      geom_histogram(binwidth = function(x) diff(range(x, na.rm = TRUE))/30) +
      theme_minimal() +
      labs(title = paste("Histogram of", pretty_names[input$field]),
           x = pretty_names[input$field],
           y = "Count") +
      scale_x_continuous(labels = scales::comma)
  })

  output$summary_stats <- renderUI({
    req(data())
    if (nrow(data()) == 0) {
      return(HTML("<p>No valid data available for analysis.</p>"))
    }

    tryCatch({
      # Calculate BIN compliance
      nuc_ambig_stat <- data() %>%
        summarise(percentage = mean(nuc_basecount >= 500 & ambig_basecount <= 6, na.rm = TRUE) * 100) %>%
        pull(percentage)

      # Calculate taxonomy match percentage with safer handling
      id_obs_stat <- data() %>%
        rowwise() %>%
        mutate(
          match = tryCatch({
            if (is.na(identification) || is.na(obs_taxon)) {
              FALSE
            } else {
              obs_taxa <- strsplit(obs_taxon, ",")[[1]]
              identification %in% obs_taxa
            }
          }, error = function(e) FALSE)
        ) %>%
        ungroup() %>%
        summarise(percentage = mean(match, na.rm = TRUE) * 100) %>%
        pull(percentage)

      # Format output
      HTML(paste(
        "<ul>",
        "<li><strong>BIN compliant records (barcode length &ge; 500 and ambiguous bases &le; 6):</strong> ",
        sprintf("%.2f%%", nuc_ambig_stat), "</li>",
        "<li><strong>Records where expected taxon is recovered through reverse taxonomy:</strong> ",
        sprintf("%.2f%%", id_obs_stat), "</li>",
        "</ul>"
      ))
    }, error = function(e) {
      HTML(paste("<p>Error calculating statistics:", htmlEscape(as.character(e)), "</p>"))
    })
  })
}

# Run the app
shinyApp(ui = ui, server = server)