# This script is for making the actual data base 
# Jennifer Moreno 6/5/25
# Load in needed packages
library(duckdb)
library(dplyr)
library(DT)

#_________________________________________________
# Define the directory path and database file name
db_dir <- "C:/Users/JMorenoRamirez/OneDrive/new/ENIGMA web-app directory/Data Output"
db_file <- "exampleduckbyjen.db"

# Create the directory if it doesn't exist
dir.create(db_dir, recursive = TRUE)

# Connect to the database
con <- dbConnect(duckdb(file.path(db_dir, db_file)))

# Getting the data from the csv (raw data)
# Headers set to true if the first row has the headers
data <-read.csv("C:/Users/JMorenoRamirez/Downloads/Copy of Meta_data_test - Sheet4.csv")



# Writing the data into the database
# The words in the Quotes are what this table of data will be named in the database (for future referencing)
dbWriteTable(con, "meta_data", data)

# Here you can see the available tables in the db by name
tables <- dbListTables(con)
for (table in tables) {
  print(table)
}

# Here you can see the data within a table by referencing it by name
data <- dbGetQuery(con, "SELECT * FROM meta_data")
datatable(data)

#______________________________________________________
# # This was for testing purposes (un-comment if needed)
# # Read the CSV file
# data <- read.csv("C:/Users/JMorenoRamirez/Downloads/counts_data.tab")
# 
# # Print the data
# print(data)
# 
# # Run a SELECT query on the table
# result <- dbGetQuery(con, "SELECT * FROM counts_data")
# 
# # Print the result
# print(result)
#_______________________________________________________
# Close the database connection after everything is done
dbDisconnect(con)
