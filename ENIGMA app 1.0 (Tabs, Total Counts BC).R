# Load packages ----------------------------------------------------------------
# app.R for testing purposes
# Interrupt is shift+ctrl+F10
# Jennifer Moreno 6/4/25 CCI DOE
# must call these every time
# library calls the packages after manual install (install.packages("packagename"))
# These packages are for general data base(db) and app-building
library(duckdb)
library(shiny)
# These packages are for app appearance
library(shinydashboard)
library(shinyWidgets)
library(bslib)
library(shinythemes)
# These are for statistical analysis
library(RSQLite)
library(tidyverse) # This package has tidyr, dplyr, and ggplot2 (among others)
library(phyloseq)
#___________________________________________
#user interface (front end)
ui <- navbarPage(
  #change the title in quotes here
  'ENIGMA Enrichments',
  theme=shinytheme("journal"),
  tabPanel("Metadata Overview"
    # Will eventually have all main graphs for metadata (from the Google drive)
    #Possible integrate the ability to select specific graphs (no data manip)
  ),
  tabPanel("Total Counts",
    # CURRENTLY WORKING ON: having the ability to select the meta data tabs specifically and see the samples and their TOTAL counts(reads)
    sidebarLayout(
      sidebarPanel(
        uiOutput("metadata_selection_ui"),
        uiOutput("submetadata_selection_ui")
      ),
      mainPanel(
        # Will eventually contain the outputted plot for Total Counts tabs
        plotOutput("metadata_barchart"),
        dataTableOutput("filtered_metadata"),
        dataTableOutput("filtered_totalcounts")
      )
    )
  ),
)
#____________________________________________________________
#back-end
server <- function (input, output, session) {
    # This opens the db given the path (for persistent memory)
    con <- dbConnect(duckdb::duckdb(), dbdir = "C:/Users/JMorenoRamirez/OneDrive/new/ENIGMA web-app directory/Data Output/exampleduckbyjen.db")
    #disconnects from the data base when app is exited 
    on.exit({(dbDisconnect(con))})
#_____________________________________________________________
##TESTING CODE HERE for meta data drop down menu based on the asv model
# Any/All manipulation for meta data mass selection -> TOTAL ASV counts / Sample ID bar chart
  #________________________________________________________
  # THIS WORKS
  # Getting the sum table from the counts table
    countsdata <- dbGetQuery(con,"SELECT * FROM counts_data")
    # Ensures the correct row names
    rownames(countsdata) <- countsdata[, 1]
    countsdata <- countsdata[, -1]
    # Ensures all values in the table are number
    countsdatanumeric<- countsdata[, sapply(countsdata, is.numeric)]
    # takes the sums of the columns
    countsums <- colSums(countsdatanumeric, na.rm=FALSE)
    # makes a data frame of the table for manip
    # naming the columns in parens
    countsums_df <- data.frame(SampleID = names(countsdata), Sum=countsums)
    countsums_df$SampleID <- gsub("X", "", countsums_df$SampleID)
  #______________________________________________________
  # THIS WORKS
  # Getting user input for meta data tab
    # Query the db to get the data from the meta data table
    metadata <-dbGetQuery(con, "SELECT * FROM meta_data")
    
    #1 Get the column(s) name input
      output$metadata_selection_ui <- renderUI({
        selectInput("metadata_selection", "Select meta-data type(s) of interest:",
          choices = colnames(metadata),
          multiple = TRUE)
     })
    #2 Get the selected column(s) value(s) input
     output$submetadata_selection_ui <- renderUI({
       # requires a selection from the meta data headers
       req(input$metadata_selection)
       # Combines the elements for sub selection into one UI element
       tagList(
         # creating the function to apply to all sub meta selections (using headers as vectors)
         lapply(input$metadata_selection, function(x) {
           # subsets the metadata values by header+ removes duplicate options
           values <- unique(metadata[[x]])
           # creates elements for each metadata type
           # paste0 makes unique element IDs for each element variant (diff header choices)
           selectInput(paste0("submetadata_selection_", x),
             paste0("Select Value(s) from ", x, ":"),
             choices = values,
             multiple = TRUE)
         })
       )
     })
  #________________________________________________________
  ## TESTING
  # Getting the filtered sample IDs into a data subset
  metadatasubset <- reactive ({
     # If no selections are made:
    if (is.null(input$metadata_selection) == TRUE || length(input$metadata_selection) == 0) {
      # returns metadata (non maniped values)
      return(metadata)
    } else {
      # initializes to non manip metadata
      subset <- metadata
      # for the column(s) selected:
      for (meta in input$metadata_selection) {
        # grabs the value(s) from selected column
        submeta_input <- input[[paste0("submetadata_selection_", meta)]]
        # if it exists and is valid in the table, grab the rows that they belong to
        if (!is.null(submeta_input) && length(submeta_input) > 0) {
          subset <- subset[subset[[meta]] %in% submeta_input, ]
        }
      }
      # if the filters don't have matching rows, let the user know
      if (nrow(subset) == 0) {
        return(data.frame(ncol=ncol(metadata), nrow = 0))
      } else {
        return(subset)
      }
    }
  })

  
  #____________________________________________________
  ## TESTING
  # Making a bar chart for the filtered data
  output$metadata_barchart <- renderPlot({
    # require some data to make the bar plot
    req(metadatasubset())
    subset <- metadatasubset()
      # If the subset is a string (no data matches the filters)
      if (nrow(subset)==0) {
        # This doesn't actually output anything
        validate(need(nrow(subset)!=0, "There are no matches in the dataset."))
      } else {
        # Generate the plot if there is data
        # find the matching columns for SampleID in the countsums data frame as is present in the subset
        sample_sums <- countsums_df[countsums_df$SampleID %in% subset$SampleID, ]
        ggplot(sample_sums, aes(x = SampleID, y = Sum)) +
          geom_bar(stat = "identity") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          labs(title = "Total Counts by Sample ID", x = "Sample ID", y = "Total Counts")
      }
  })
  # Render the filtered metadata table
  output$filtered_metadata <- renderDT({
    req(metadatasubset())
    datatable(metadatasubset(),
      # title
      caption = "Filtered Metadata Matches",
      # these options allow for the scroll bar and fitting to screen
      options = list(
        scrollX = TRUE,
        autoWidth = TRUE
        ))
  })
  
  # Render the filtered total counts table
  output$filtered_totalcounts <- renderDT({
    req(metadatasubset())
    subset <- metadatasubset()
    # If subset has no matches:
    if (nrow(subset) == 0) {
      # returns an empty data frame
      return(data.frame(matrix(ncol = 2, nrow = 0, dimnames = list(NULL, c("SampleID", "Sum")))))
    } else {
      # otherwise find the matching vals and output
      sample_sums <- countsums_df[countsums_df$SampleID %in% subset$SampleID, ]
      datatable(sample_sums, caption = "Filtered Total Count Matches", rownames = FALSE)
    }
  })
}
#___________________________________________________________________
# Compiling all the components
shinyApp(ui, server)
# AS of 6/13/25 This "works" but require Jen Testing and interaction