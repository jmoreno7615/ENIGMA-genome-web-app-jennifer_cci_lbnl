# This script is for making the actual data base, tables, and viewing db contents
# Jennifer Moreno 6/5/25
# Load in needed packages
library(duckdb)
library(dplyr)
library(DT)
library(RSQLite)
library(tidyverse) # This package has tidyr, dplyr, and ggplot2 (among others)
library(phyloseq)
library(DT)

#_________________________________________________
#CREATING AND CONNECTING TO DB
  # Define the directory path and database file name
  # REPLACE with your own directory path
  db_dir <- "C:/Users/JMorenoRamirez/OneDrive/new/ENIGMA web-app directory/Data Output"
  # REPLACE with your own name.db for the data base
  db_file <- "exampleduckbyjenTEST.db"

  # Create the directory if it doesn't exist
  dir.create(db_dir, recursive = TRUE)

  # Connect to the database - no need to replace anything here
  con <- dbConnect(duckdb(file.path(db_dir, db_file)))

#_____________________________
#CREATING DB TABLES
  # Getting the data from the csv (raw data)
  # Headers set to true if the first row has the headers
  #METADATA
    # REPLACE with your own filename
    metadata <-read.csv("C:/Users/JMorenoRamirez/Downloads/Copy of Meta_data_test - Sheet4(1).csv")
    # Writing the data into the database
    # The words in the Quotes are what this table of data will be named in the database (for future referencing)
    dbWriteTable(con, "meta_data", metadata, overwrite = TRUE) #include append = TRUE argument if adding to the existing table
  
  #COUNTSDATA
    # REPLACE with your own filename
    countsdata <-read.csv("C:/Users/JMorenoRamirez/Downloads/counts_data_2 (1).csv")
    colnames(countsdata) <- gsub("X", "", colnames(countsdata))
    dbWriteTable(con, "counts_data", countsdata, overwrite = TRUE) #include append = TRUE argument if adding to the existing table
  
  #TAXADATA
    # REPLACE with your own filename
    taxadata <- "C:/Users/JMorenoRamirez/Downloads/tax_data - tax_data.csv"
    
    #filling in missing values
    taxa_df <- read.csv(taxadata, row.names = 1)
    taxa_df <- rownames_to_column(taxa_df, "ASV")
    # Function to find the most specific taxonomic designation
    fill_missing_taxa <- function(row) {
      # Find the indices of filled cells in each row
      filled_indices <- which(row != "")
      # If there are no filled cells, return the row as totally empty
      if (length(filled_indices) == 0) {
        return(row)
      #if there are filled cells
      } else {
        # Get the most specific taxonomic designation
        most_specific_taxa <- row[filled_indices[length(filled_indices)]]
        # Fill in the missing values with the most specific taxonomic designation available
        row[row == ""] <- most_specific_taxa
        return(row)
      }
    }
    
    filled_taxadata <- apply(taxa_df, 1, fill_missing_taxa)
    # Convert the result back to a data frame, ensuring rows and columns are correct
    filled_taxadata <- as.data.frame((filled_taxadata), stringsAsFactors = FALSE)
    # fixing the formatting of columns and headers
    final_taxadata <- data.frame(t(filled_taxadata[-1]))
    row.names(final_taxadata) <- NULL
    dbWriteTable(con, "taxa_data", final_taxadata, overwrite = TRUE) #include append = TRUE argument if adding to the existing table
   
  #TOTAL COUNTS  
    # Get the data
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
    
    #For debugging/verification before appending
    #datatable(countsums_df)
    dbWriteTable(con, "countsums_data", countsums_df, overwrite = TRUE)
    
    
  #RELATIVE ABUNDANCE
    # Get data from Counts 
    countsdata <- dbGetQuery(con,"SELECT * FROM counts_data")
    
    # Getting correct row names
    rownames(countsdata) <- countsdata[, 1]
    countsdata <- countsdata[, -1]
    countsdatanumeric<- countsdata[, sapply(countsdata, is.numeric)]
    countsums <- colSums(countsdatanumeric, na.rm=TRUE)
    countsums_df <- data.frame(SampleID = names(countsdata), Sum=countsums)
    countsums_df$SampleID <- gsub("X", "", countsums_df$SampleID)
    colnames(countsdata) <- gsub("X", "", colnames(countsdata))
    
    # Assuming countsums_df has columns 'SampleID' and 'Sum'
    # Assuming countsdata has ASV rownames and SampleIDs as column headers
    # Ensure SampleID is a character vector in both data frames
    countsums_df$SampleID <- as.character(countsums_df$SampleID)
    colnames(countsdata) <- as.character(colnames(countsdata))
    
    # Create a named vector of sums for quick lookup
    sums_lookup <- setNames(countsums_df$Sum, countsums_df$SampleID)
    # Initialize relative_abundance with the same dimensions as countsdata
    relative_abundance <- matrix(0, nrow=nrow(countsdata), ncol=ncol(countsdata))
    
    # Get column and row names
    sample_ids <- colnames(countsdata)
    asvs <- rownames(countsdata)
    
    # Extract sum values for all sample IDs at once
    sum_values <- sapply(sample_ids, function(id) sums_lookup[[id]])
    for (i in seq_along(sample_ids)) {
      sample_id <- sample_ids[i]
      sum_value <- sum_values[i]
      
      # Check if the sum is not NA and not zero
      if (!is.na(sum_value) && sum_value != 0) {
        # Vectorized operation for all ASVs
        counts <- countsdata[, sample_id]
        non_zero_counts <- counts > 0
        relative_abundance[non_zero_counts, i] <- counts[non_zero_counts] / sum_value
      }else {
        # If sum is zero or NA, set relative abundance to 0 for all ASVs
        relative_abundance[, i] <- 0
      }
    }
    relative_abundance_df <- as.data.frame(relative_abundance2.0)
    colnames(relative_abundance_df) <- sample_ids
    rownames(relative_abundance_df) <- asvs
    relative_abundance_df <- rownames_to_column(relative_abundance_df, "ASV")
    dbWriteTable(con, "relabund_data", relative_abundance_df, overwrite = TRUE)
    
#_____________________________
#VIEWING DB CONTENTS
  #TABLE NAMES
  tables <- dbListTables(con)
  for (table in tables) {
    print(table)
  }
  #TABLE CONTENTS
  #change the FROM table to the table name desired
  data <- dbGetQuery(con, "SELECT * FROM taxa_data")
  datatable(head(data))
  
  #ALTERATIONS:
  #if you need to change a column name
  #dbExecute(con, "ALTER TABLE relabund_data RENAME COLUMN row_names TO ASV;")
#___________________
dbDisconnect(con)