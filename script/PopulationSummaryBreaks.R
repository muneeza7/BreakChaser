#Set working directory
setwd("./bams/tmp_bams/bp/")

#####################################################################################################################################################

## load library
library(pacman)

p_load("tidyverse", 
       "Biostrings", 
       "BSgenome.Hsapiens.UCSC.hg19", 
       "dplyr", 
       "stringr",
       "randomForest",
       "DT", "ggplot2", "data.table")

pacman::p_load(readxl, tidyverse, xlsx, data.table, dplyr, purrr, matrixStats, ggplot2, openxlsx, magrittr, gapminder)
options(scipen = 999)

### Functions
new_df <- function(x){
  
  cols_to_retain <- c("chr", "donor_pos", "acceptor_pos", "geneName_5p", "geneName_3p", "DelSize", "avg_insert", "counts_wrt_5p_3p", 
                      "total_5p_counts", "total_3p_counts", "strand_5p", "strand_3p", "location_5p", "locNum_5p", 
                      "location_3p", "locNum_3p", "lastExDist_5p", "nextExDist_5p", "lastExDist_3p", "nextExDist_3p")
  x <- x %>% dplyr::select(all_of(cols_to_retain)) %>%
    mutate(chr__donor_pos__acceptor_pos__geneName_5p__geneName_3p__DelSize__strand_5p__strand_3p__location_5p__locNum_5p__location_3p__locNum_3p__lastExDist_5p__nextExDist_5p__lastExDist_3p__nextExDist_3p = paste(
      chr, donor_pos, acceptor_pos, geneName_5p, geneName_3p, DelSize, 
      strand_5p, strand_3p, location_5p, locNum_5p, location_3p, locNum_3p,
      lastExDist_5p, nextExDist_5p, lastExDist_3p, nextExDist_3p,
                        sep = "__"),
           .before = chr)
  x
}
#### region of interest
roi <- function(x){
  
  x <- x %>% 
    dplyr::filter(chr %in% c("chr7", "chr9", "chr12", "chr13")) %>%
    dplyr::mutate(flag = ifelse(((acceptor_pos > 50168000 & acceptor_pos < 50519000) | (acceptor_pos > 48950000 & acceptor_pos < 49175000) | 
                            (acceptor_pos > 21800000 & acceptor_pos < 22300000) | (acceptor_pos > 36802000 & acceptor_pos < 37053500) | 
                            (acceptor_pos > 92146000 & acceptor_pos < 92539400) | is.na(acceptor_pos)) & 
                           ((donor_pos > 50168000 & donor_pos < 50518000) | (donor_pos > 48950000 & donor_pos < 49175000) | 
                              (donor_pos > 21800000 & donor_pos < 22300000) | (donor_pos > 36802000 & donor_pos < 37053500) | 
                              (donor_pos > 92146000 & donor_pos < 92539400) | is.na(donor_pos)), paste("keep"), paste("fp") )) %>%
    dplyr::filter(flag == "keep") %>%
    dplyr::select(-flag)
  x
}

#Defining function to calculate stats (11) ## For comparing multiple datasets
calc_stat <- function(x){
  x1 <- x %>% remove_rownames() %>% 
    column_to_rownames(var = "chr__donor_pos__acceptor_pos__geneName_5p__geneName_3p__DelSize__strand_5p__strand_3p__location_5p__locNum_5p__location_3p__locNum_3p__lastExDist_5p__nextExDist_5p__lastExDist_3p__nextExDist_3p") %>%
    dplyr::select(starts_with("count")) %>%
    type_convert(.) %>%
    transmute(`#PatientsWithVarCall` = rowSums(!is.na(.)))
  
  x2 <- x %>% remove_rownames() %>% 
    column_to_rownames(var = "chr__donor_pos__acceptor_pos__geneName_5p__geneName_3p__DelSize__strand_5p__strand_3p__location_5p__locNum_5p__location_3p__locNum_3p__lastExDist_5p__nextExDist_5p__lastExDist_3p__nextExDist_3p") %>%
    dplyr::select(starts_with("total_5p")) %>%
    type_convert(.) %>% rowwise() %>%
    transmute(
              max_count_5p = max(c_across(everything(starts_with("total_5p"))), na.rm = T))
  
  x3 <- x %>% remove_rownames() %>% 
    column_to_rownames(var = "chr__donor_pos__acceptor_pos__geneName_5p__geneName_3p__DelSize__strand_5p__strand_3p__location_5p__locNum_5p__location_3p__locNum_3p__lastExDist_5p__nextExDist_5p__lastExDist_3p__nextExDist_3p") %>%
    dplyr::select(starts_with("total_3p")) %>%
    type_convert(.) %>% rowwise() %>%
    transmute(
      max_count_3p = max(c_across(everything(starts_with("total_3p"))), na.rm = T))
  
  
  
  x5 <- x %>% remove_rownames() %>%
    column_to_rownames(var = "chr__donor_pos__acceptor_pos__geneName_5p__geneName_3p__DelSize__strand_5p__strand_3p__location_5p__locNum_5p__location_3p__locNum_3p__lastExDist_5p__nextExDist_5p__lastExDist_3p__nextExDist_3p")
  x <- cbind(x5, x2, x3, x1)
}
# Genes of Interest
goi <- c( "ASXL1", "ATM","BCOR","BCORL1","BTG1","CALR","CDKN2A","CDKN2B","CUX1","DNMT3A","ETV6",
          "EZH2","GATA2","IKZF1","JAK2","KDM6A","KIT","KMT2D","PAX5","PHF6","PPM1D","RB1",
          "RUNX1","SETD1B","SETD2","SF3B1","SMC3","SRSF2","STAG2","TET2","TP53","UBE2A","WT1", "MTAP", "CDKN2A-AS1","CDKN2B-AS1", 
          "EP300", "NUP98", "CREBBP", "RAD21", "RCBTB2", "CBFB", "MYH11")


## Reading files:
all_samples <- list.files(pattern = "*_breaks.csv")
list_all_samples <- sapply(all_samples, read.csv, stringsAsFactors = F, simplify = F)

#Creating a df_list with equal number of columns:
shorterDF_list <- sapply(list_all_samples, new_df, simplify = F)

## filtering on region of interest
shorterDF_list1 <- sapply(shorterDF_list, roi, simplify = F)


#merging the list of dfs into a single df:
mergedDF <- bind_rows(shorterDF_list1, .id = "id")

mergedDF1 <- mergedDF %>%
  filter(geneName_5p %in% c(goi, "intergenic", NA) | geneName_3p %in% c(goi, "intergenic", NA)) %>%
  filter(!(geneName_5p %in% "intergenic" &  geneName_3p %in% "intergenic"))

#pivot_wider according to the IDs of dataframes
mergedDF1 %>% 
  pivot_wider(id_cols = chr__donor_pos__acceptor_pos__geneName_5p__geneName_3p__DelSize__strand_5p__strand_3p__location_5p__locNum_5p__location_3p__locNum_3p__lastExDist_5p__nextExDist_5p__lastExDist_3p__nextExDist_3p,
              names_from = id,
              values_from = c(counts_wrt_5p_3p,
                              total_5p_counts,
                              total_3p_counts,
                              avg_insert),
              names_vary = "slowest") -> df_wider
#Calculating some stats of the variants data and how many times are these variants been called (# of Patients & runs)
df_wider %>% calc_stat() %>% arrange(desc(`#PatientsWithVarCall`)) -> final_stat

final_stat1 <- final_stat %>% dplyr::select(`#PatientsWithVarCall`,  "max_count_5p", "max_count_3p", everything())
final_stat1 %>% rownames_to_column(var = "chr__donor_pos__acceptor_pos__geneName_5p__geneName_3p__DelSize__strand_5p__strand_3p__location_5p__locNum_5p__location_3p__locNum_3p__lastExDist_5p__nextExDist_5p__lastExDist_3p__nextExDist_3p") -> final_stat2

finale <- final_stat2 %>% separate(chr__donor_pos__acceptor_pos__geneName_5p__geneName_3p__DelSize__strand_5p__strand_3p__location_5p__locNum_5p__location_3p__locNum_3p__lastExDist_5p__nextExDist_5p__lastExDist_3p__nextExDist_3p, 
                                   c("chr", "donor_pos","acceptor_pos", "geneName_5p", "geneName_3p", "DelSize", "strand_5p", 
                                     "strand_3p", "location_5p", "locNum_5p", "location_3p", "locNum_3p", "lastExDist_5p", 
                                     "nextExDist_5p", "lastExDist_3p", "nextExDist_3p"),
                                   sep = "__") %>% suppressWarnings(type_convert(.)) %>%
  dplyr::select(-lastExDist_5p, -lastExDist_3p, -nextExDist_5p, -nextExDist_3p)


# Set working directory to extract the base directory name
setwd("../../../../.")

# Extract the folder name from the working directory
foldername <- basename(getwd())

# change working directory to where the final output file has to be stored
setwd("./BreakHunter/bams/tmp_bams/bp/")

# Create a new workbook
wb <- createWorkbook()

# Add a worksheet to the workbook
addWorksheet(wb, "Sheet1")

# copy final df name to a new variable to write into the files
df <- finale %>% type_convert(.)

# Get the column names
column_names <- colnames(df)

# Define the color to apply to the columns starting with "unique"
color <- "#ccccff"  # You can change this to any desired color

# hs <- createStyle(textDecoration = "BOLD", valign = "top", wrapText = T)
# writeData(wb, "Sheet1", finale, startCol = i, startRow = 1, colNames = T, headerStyle = hs)

# Iterate over the column names and add data to the worksheet
for (i in 1:length(column_names)) {
  column_name <- column_names[i]
  column_data <- df[[column_name]]
  
  # Write the column name and data to the worksheet
  hs <- createStyle(textDecoration = "BOLD", valign = "top", wrapText = T)
  writeData(wb, "Sheet1", column_name, startCol = i, startRow = 1, colNames = T, headerStyle = hs)
  writeData(wb, "Sheet1", column_data, startCol = i, startRow = 1, colNames = T, headerStyle = hs)
  
  if (grepl("^#Patients", column_name)) {
    cellStyle <- createStyle(fontColour = "blue")
    addStyle(wb, "Sheet1", cellStyle, rows = 1:nrow(df)+1, cols = i)
  }
  #Check if the column name starts with "unique" and apply the color
  if (grepl("^counts", column_name)) {
    cellStyle <- createStyle(fgFill =  color)
    addStyle(wb, "Sheet1", cellStyle, rows = 1:nrow(df)+1, cols = i)
  }
}

# Save the workbook to a file
writeDataTable(wb, "Sheet1", df, colNames = T, headerStyle = hs, rowNames = F)
saveWorkbook(wb, paste0(foldername, "_population_summary_breaks.xlsx"), overwrite = TRUE)

