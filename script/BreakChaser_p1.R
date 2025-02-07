rm(list=ls())

## load library
library(pacman)
p_load("stringr","tidyverse", "openxlsx", "stringr", "annotate", "GenomicAlignments",
       "BSgenome.Hsapiens.UCSC.hg19", "GenomicFeatures", "dplyr", "GenomicRanges", "readr", "purrr")
options(scipen = 999)
############################################################ Functions ###############################################################################

# Function to extract values of 'S'
extract_S_values <- function(cigar) {
  s_values <- str_extract_all(cigar, "[0-9]+S")[[1]]
  
  if (length(s_values) >= 3) {
    return(s_values[1:3])
  } else if (length(s_values) == 2) {
    return(c(s_values[1], s_values[2], 0))
  } else if (length(s_values) == 1) {
    return(c(s_values[1], 0, 0))
  } else {
    return(c(0, 0, 0))
  }
}

# Function to extract values of 'M'
extract_M_values <- function(cigar) {
  m_values <- str_extract_all(cigar, "[0-9]+M")[[1]]
  
  if (length(m_values) >= 3) {
    return(m_values[1:3])
  } else if (length(m_values) == 2) {
    return(c(m_values[1], m_values[2], 0))
  } else if (length(m_values) == 1) {
    return(c(m_values[1], 0, 0))
  } else {
    return(c(0, 0, 0))
  }
}

# Function to extract values of 'I'
extract_I_values <- function(cigar) {
  i_values <- str_extract_all(cigar, "[0-9]+I")[[1]]
  
  if (length(i_values) == 2) {
    return(c(i_values[1], i_values[2], 0))
  } else if (length(i_values) == 1) {
    return(c(i_values[1], 0))
  } else {
    return(c(0, 0))
  }
}

# Function to extract values of 'D'
extract_D_values <- function(cigar) {
  d_values <- str_extract_all(cigar, "[0-9]+D")[[1]]
  
  if (length(d_values) == 2) {
    return(c(d_values[1], d_values[2], 0))
  } else if (length(d_values) == 1) {
    return(c(d_values[1], 0))
  } else {
    return(c(0, 0))
  }
}

# Function to find the start positions of motif sequences in another sequence
find_motif_start_positions <- function(reference_sequence, motif_sequences) {
  start_positions_list <- lapply(motif_sequences, function(motif_sequence) {
    match_positions <- regexpr(motif_sequence, reference_sequence)
    start_positions <- ifelse(match_positions > -1, match_positions, NA)
    return(start_positions)
  })
  return(start_positions_list)
}

# Make the annotating function. It will annotate the intervals with gene_Ids
annotateIntervals <-function(intervals, txdb){
  stopifnot(is(intervals, "GRanges"), is(txdb, "TxDb"))
  anno = genes(txdb)
  olaps = findOverlaps(intervals, anno)
  mcols(olaps)$gene_id = genes$gene_id[subjectHits(olaps)]
  intervals_factor = factor(queryHits(olaps), levels=seq_len(queryLength(olaps)))
  intervals$gene_id = splitAsList(mcols(olaps)$gene_id, intervals_factor)
  intervals
}

######################################################################################################################################################

## Reading files with soft-clipped reads

files <- list.files(path="./bams/tmp_bams", 
                    pattern="_reads.tsv", full.names = T)

# Using the reference fasta and fai files (we use for mapping) I created a .bed file with Chr info and lengths and removed other
# columns except the 1st 3 and saved the file as .csv
### Creating txdb object to annotate the coordinates:
# Reading the bed file:
bed <- read.table("./database/WholeGenomeFasta/human_hg19.csv", header = F, sep = ",")

# Used the bed file info to generate chrominfo obj
chrominfo <- data.frame(chrom = bed[,1],
                        length = bed[, 3]) %>% mutate(is_circular = FALSE)
# Copied RefSeq_hg19 url from Google (NCBI)
metadata <- data.frame(name = "Resource URL",
                       value = paste0("https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.25/"))
# Putting the objects created above in the following function to create the txdb obj
txdb <- makeTxDbFromGFF(file = "./database/RefSeq.gtf",
                        format = "gtf",
                        chrominfo = chrominfo,
                        dataSource = "entrez",
                        organism = "Homo sapiens",
                        metadata = metadata)

# Make an header for the data.frame "myIntervals" containing
# the coordinates files
header <- c("Chromosome", "Start", "End")



for (k in 1:length(files)) {
  
  sc <- read.table(files[k], sep = "\t", header = F, stringsAsFactors =  F)
  
  ## extract id
  id1 <- str_extract(files[k], "[a-zA-Z0-9-_\\.]+_soft_clipped_reads.tsv")
  id <- str_trim(gsub("_soft_clipped_reads.tsv", "", id1))
  
  print(paste("processing ...  ", id1, sep=""))
  
  sc %>% 
    mutate(mateOr = ifelse(V4==V8, paste("same"), paste("opposite")),
           readStrand = ifelse(V9>0, paste("+"), "-")) %>%
    filter(V3 %in% c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X", "Y"))-> sc_0
  
  sc_0 <- sc_0[!grepl("N", sc_0$V6),] %>%
    mutate(s_count = str_count(V6, "S"),
           m_count = str_count(V6, "M"),
           d_count = str_count(V6, "D"),
           i_count = str_count(V6, "I"))
  
  ### pair-mates orientation is same and readstrand is +ve
  sc_0.1 <- sc_0 %>%
    filter(mateOr %in% "same" & readStrand %in% "+")
  ### pair-mates orientation is same and readstrand is -ve
  sc_0.2 <- sc_0 %>% 
    filter(mateOr %in% "same" & readStrand %in% "-")
  
  ### pair-mates in opp orientation and +ve read strand
  sc_0.3 <- sc_0 %>%
    filter(mateOr %in% "opposite" & readStrand %in% "+")
  
  ### pair-mates in opp orientation and +ve read strand
  sc_0.4 <- sc_0 %>%
    filter(mateOr %in% "opposite" & readStrand %in% "-")
  
  rm(sc)
  rm(sc_0)
  rm(bed)
  rm(chrominfo)
  rm(metadata)
  
  # Apply function to dataframe using purrr::map
  sc_0.1.2 <- sc_0.1 %>%
    mutate(all_s_values = purrr::map(V6, extract_S_values),
           first_s_value = as.numeric(sapply(all_s_values, function(x) gsub("S", "", x[[1]], ignore.case = T))),
           second_s_value = as.numeric(sapply(all_s_values, function(x) gsub("S", "", x[[2]], ignore.case = T)))) %>%
    mutate(all_m_values = purrr::map(V6, extract_M_values),
           first_m_value = as.numeric(sapply(all_m_values, function(x) gsub("M", "", x[[1]], ignore.case = T))),
           second_m_value = as.numeric(sapply(all_m_values, function(x) gsub("M", "", x[[2]], ignore.case = T))),
           third_m_value = as.numeric(sapply(all_m_values, function(x) gsub("M", "", x[[3]], ignore.case =T)))) %>%
    mutate(all_i_values = purrr::map(V6, extract_I_values),
           first_i_value = as.numeric(sapply(all_i_values, function(x) gsub("I", "", x[[1]], ignore.case = T))),
           second_i_value = as.numeric(sapply(all_i_values, function(x) gsub("I", "", x[[2]], ignore.case = T)))) %>%
    mutate(all_d_values = purrr::map(V6, extract_D_values),
           first_d_value = as.numeric(sapply(all_d_values, function(x) gsub("D", "", x[[1]], ignore.case = T))),
           second_d_value = as.numeric(sapply(all_d_values, function(x) gsub("D", "", x[[2]], ignore.case = T))))
  
  sc_0.2.2 <- sc_0.2 %>%
    mutate(all_s_values = purrr::map(V6, extract_S_values),
           first_s_value = as.numeric(sapply(all_s_values, function(x) gsub("S", "", x[[1]], ignore.case = T))),
           second_s_value = as.numeric(sapply(all_s_values, function(x) gsub("S", "", x[[2]], ignore.case = T)))) %>%
    mutate(all_m_values = purrr::map(V6, extract_M_values),
           first_m_value = as.numeric(sapply(all_m_values, function(x) gsub("M", "", x[[1]], ignore.case = T))),
           second_m_value = as.numeric(sapply(all_m_values, function(x) gsub("M", "", x[[2]], ignore.case = T))),
           third_m_value = as.numeric(sapply(all_m_values, function(x) gsub("M", "", x[[3]], ignore.case =T)))) %>%
    mutate(all_i_values = purrr::map(V6, extract_I_values),
           first_i_value = as.numeric(sapply(all_i_values, function(x) gsub("I", "", x[[1]], ignore.case = T))),
           second_i_value = as.numeric(sapply(all_i_values, function(x) gsub("I", "", x[[2]], ignore.case = T)))) %>%
    mutate(all_d_values = purrr::map(V6, extract_D_values),
           first_d_value = as.numeric(sapply(all_d_values, function(x) gsub("D", "", x[[1]], ignore.case = T))),
           second_d_value = as.numeric(sapply(all_d_values, function(x) gsub("D", "", x[[2]], ignore.case = T))))
  
  sc_0.3.2 <- sc_0.3 %>%
    mutate(all_s_values = purrr::map(V6, extract_S_values),
           first_s_value = as.numeric(sapply(all_s_values, function(x) gsub("S", "", x[[1]], ignore.case = T))),
           second_s_value = as.numeric(sapply(all_s_values, function(x) gsub("S", "", x[[2]], ignore.case = T))))%>%
    mutate(all_m_values = purrr::map(V6, extract_M_values),
           first_m_value = as.numeric(sapply(all_m_values, function(x) gsub("M", "", x[[1]], ignore.case = T))),
           second_m_value = as.numeric(sapply(all_m_values, function(x) gsub("M", "", x[[2]], ignore.case = T))),
           third_m_value = as.numeric(sapply(all_m_values, function(x) gsub("M", "", x[[3]], ignore.case =T)))) %>%
    mutate(all_i_values = purrr::map(V6, extract_I_values),
           first_i_value = as.numeric(sapply(all_i_values, function(x) gsub("I", "", x[[1]], ignore.case = T))),
           second_i_value = as.numeric(sapply(all_i_values, function(x) gsub("I", "", x[[2]], ignore.case = T)))) %>%
    mutate(all_d_values = purrr::map(V6, extract_D_values),
           first_d_value = as.numeric(sapply(all_d_values, function(x) gsub("D", "", x[[1]], ignore.case = T))),
           second_d_value = as.numeric(sapply(all_d_values, function(x) gsub("D", "", x[[2]], ignore.case = T))))
  
  sc_0.4.2 <- sc_0.4 %>%
    mutate(all_s_values = purrr::map(V6, extract_S_values),
           first_s_value = as.numeric(sapply(all_s_values, function(x) gsub("S", "", x[[1]], ignore.case = T))),
           second_s_value = as.numeric(sapply(all_s_values, function(x) gsub("S", "", x[[2]], ignore.case = T)))) %>%
    mutate(all_m_values = purrr::map(V6, extract_M_values),
           first_m_value = as.numeric(sapply(all_m_values, function(x) gsub("M", "", x[[1]], ignore.case = T))),
           second_m_value = as.numeric(sapply(all_m_values, function(x) gsub("M", "", x[[2]], ignore.case = T))),
           third_m_value = as.numeric(sapply(all_m_values, function(x) gsub("M", "", x[[3]], ignore.case =T)))) %>%
    mutate(all_i_values = purrr::map(V6, extract_I_values),
           first_i_value = as.numeric(sapply(all_i_values, function(x) gsub("I", "", x[[1]], ignore.case = T))),
           second_i_value = as.numeric(sapply(all_i_values, function(x) gsub("I", "", x[[2]], ignore.case = T)))) %>%
    mutate(all_d_values = purrr::map(V6, extract_D_values),
           first_d_value = as.numeric(sapply(all_d_values, function(x) gsub("D", "", x[[1]], ignore.case = T))),
           second_d_value = as.numeric(sapply(all_d_values, function(x) gsub("D", "", x[[2]], ignore.case = T))))
  
  # Removing the dataframes to clear up some memory
  rm(sc_0.1)
  rm(sc_0.2)
  rm(sc_0.3)
  rm(sc_0.4)
  
  # Extract first character of CIGAR string
  sc_0.1.3 <- sc_0.1.2 %>%
    mutate(first_alpha = substr(gsub("([^SMIDN]*)([SMIDN]).*", "\\2", V6), 1, 1))
  
  # Extract first character of CIGAR string
  sc_0.2.3 <- sc_0.2.2 %>%
    mutate(first_alpha = substr(gsub("([^SMIDN]*)([SMIDN]).*", "\\2", V6), 1, 1))
  
  # Extract first character of CIGAR string
  sc_0.3.3 <- sc_0.3.2 %>%
    mutate(first_alpha = substr(gsub("([^SMIDN]*)([SMIDN]).*", "\\2", V6), 1, 1))
  
  # Extract first character of CIGAR string
  sc_0.4.3 <- sc_0.4.2 %>%
    mutate(first_alpha = substr(gsub("([^SMIDN]*)([SMIDN]).*", "\\2", V6), 1, 1))
  
  # Removing the dataframes to clear up some memory
  rm(sc_0.1.2)
  rm(sc_0.2.2)
  rm(sc_0.3.2)
  rm(sc_0.4.2)
  
  # Extract last alphabet in CIGAR string
  sc_0.1.4 <- sc_0.1.3 %>%
    mutate(last_alpha = substr(gsub(".*([SMIDN])[^SMIDN]*$", "\\1", V6), nchar(gsub(".*([SMIDN])[^SMIDN]*$", "\\1", V6)), nchar(gsub(".*([SMIDN])[^SMIDN]*$", "\\1", V6))))
  
  # Extract last alphabet in CIGAR string
  sc_0.2.4 <- sc_0.2.3 %>%
    mutate(last_alpha = substr(gsub(".*([SMIDN])[^SMIDN]*$", "\\1", V6), nchar(gsub(".*([SMIDN])[^SMIDN]*$", "\\1", V6)), nchar(gsub(".*([SMIDN])[^SMIDN]*$", "\\1", V6))))
  
  # Extract last alphabet in CIGAR string
  sc_0.3.4 <- sc_0.3.3 %>%
    mutate(last_alpha = substr(gsub(".*([SMIDN])[^SMIDN]*$", "\\1", V6), nchar(gsub(".*([SMIDN])[^SMIDN]*$", "\\1", V6)), nchar(gsub(".*([SMIDN])[^SMIDN]*$", "\\1", V6))))
  
  # Extract last alphabet in CIGAR string
  sc_0.4.4 <- sc_0.4.3 %>%
    mutate(last_alpha = substr(gsub(".*([SMIDN])[^SMIDN]*$", "\\1", V6), nchar(gsub(".*([SMIDN])[^SMIDN]*$", "\\1", V6)), nchar(gsub(".*([SMIDN])[^SMIDN]*$", "\\1", V6))))
  
  # Removing the dataframes to clear up some memory
  rm(sc_0.1.3)
  rm(sc_0.2.3)
  rm(sc_0.3.3)
  rm(sc_0.4.3)
  # Calculating Breakpoints from CIGAR info
  break_p0.1 <- sc_0.1.4 %>% type.convert() %>%
    mutate(bp = ifelse(first_alpha == "S" & last_alpha == "M",
                       V4,
                       ifelse(first_alpha == "M" & last_alpha == "S",
                              V4 + first_m_value + second_m_value + third_m_value + first_i_value + second_i_value + first_d_value +
                                second_d_value, 
                              ifelse(first_alpha == "S" & last_alpha == "S" & first_s_value < second_s_value, 
                                     V4 + first_s_value + first_m_value + second_m_value + third_m_value + first_i_value + second_i_value + first_d_value +
                                       second_d_value, ifelse(first_alpha == "S" & last_alpha == "S" & first_s_value > second_s_value,
                                                              V4 + first_s_value, 0
                                       ))))) %>%
    filter(first_s_value != second_s_value)  %>%
    mutate(flag = ifelse(s_count == 1 & first_s_value < 5, "art", "")) %>%
    dplyr::filter(!flag %in% "art")
  
  break_p0.2 <- sc_0.2.4 %>% type.convert() %>%
    mutate(bp = ifelse(first_alpha == "S" & last_alpha == "M",
                       V4,
                       ifelse(first_alpha == "M" & last_alpha == "S",
                              V4 + first_m_value + second_m_value + third_m_value + first_i_value + second_i_value + first_d_value +
                                second_d_value, 
                              ifelse(first_alpha == "S" & last_alpha == "S" & first_s_value < second_s_value, 
                                     V4 + first_s_value + first_m_value + second_m_value + third_m_value + first_i_value + second_i_value + first_d_value +
                                       second_d_value, ifelse(first_alpha == "S" & last_alpha == "S" & first_s_value > second_s_value,
                                                              V4 + first_s_value, 0
                                       ))))) %>%
    filter(first_s_value != second_s_value)  %>%
    mutate(flag = ifelse(s_count == 1 & first_s_value < 5, "art", "")) %>%
    dplyr::filter(!flag %in% "art")
  
  break_p0.3 <- sc_0.3.4 %>% type.convert() %>%
    mutate(bp = ifelse(first_alpha == "S" & last_alpha == "M",
                       V4,
                       ifelse(first_alpha == "M" & last_alpha == "S",
                              V4 + first_m_value + second_m_value + third_m_value + first_i_value + second_i_value + first_d_value +
                                second_d_value, 
                              ifelse(first_alpha == "S" & last_alpha == "S" & first_s_value < second_s_value, 
                                     V4 + first_s_value + first_m_value + second_m_value + third_m_value + first_i_value + second_i_value + first_d_value +
                                       second_d_value, ifelse(first_alpha == "S" & last_alpha == "S" & first_s_value > second_s_value,
                                                              V4 + first_s_value, 0
                                       ))))) %>%
    filter(first_s_value != second_s_value)  %>%
    mutate(flag = ifelse(s_count == 1 & first_s_value < 5, "art", "")) %>%
    dplyr::filter(!flag %in% "art")
  
  break_p0.4 <- sc_0.4.4 %>% type.convert() %>%
    mutate(bp = ifelse(first_alpha == "S" & last_alpha == "M",
                       V4,
                       ifelse(first_alpha == "M" & last_alpha == "S",
                              V4 + first_m_value + second_m_value + third_m_value + first_i_value + second_i_value + first_d_value +
                                second_d_value, 
                              ifelse(first_alpha == "S" & last_alpha == "S" & first_s_value < second_s_value, 
                                     V4 + first_s_value + first_m_value + second_m_value + third_m_value + first_i_value + second_i_value + first_d_value +
                                       second_d_value, ifelse(first_alpha == "S" & last_alpha == "S" & first_s_value > second_s_value,
                                                              V4 + first_s_value, 0
                                       ))))) %>%
    filter(first_s_value != second_s_value)  %>%
    mutate(flag = ifelse(s_count == 1 & first_s_value < 5, "art", "")) %>%
    dplyr::filter(!flag %in% "art")
  
  # Removing the dataframes to clear up some memory
  rm(sc_0.1.4)
  rm(sc_0.2.4)
  rm(sc_0.3.4)
  rm(sc_0.4.4)
  
  break_p0.1.5p <- break_p0.1 %>%
    filter(first_alpha %in% "M")
  
  break_p0.1.3p <- break_p0.1 %>%
    filter(first_alpha %in% "S")
  
  break_p0.2.5p <- break_p0.2 %>%
    filter(first_alpha %in% "M")
  
  break_p0.2.3p <- break_p0.2 %>%
    filter(first_alpha %in% "S")
  
  break_p0.3.5p <- break_p0.3 %>%
    filter(first_alpha %in% "M")
  
  break_p0.3.3p <- break_p0.3 %>%
    filter(first_alpha %in% "S")
  
  break_p0.4.5p <- break_p0.4 %>%
    filter(first_alpha %in% "M")
  
  break_p0.4.3p <- break_p0.4 %>%
    filter(first_alpha %in% "S")
  
  # Removing the dataframes to clear up some memory
  rm(break_p0.1)
  rm(break_p0.2)
  rm(break_p0.3)
  rm(break_p0.4)
  
  bp01.5p <- break_p0.1.5p %>%
    dplyr::select(V1, V3, V4, V9, bp) %>%
    dplyr::rename(Len=V9)
  
  
  
  bp02.5p <- break_p0.2.5p %>%
    dplyr::select(V1, V3, V4, V9) %>%
    dplyr::rename(Len=V9) %>%
    mutate(Len = abs(Len))
  
  # Combining the dataframes based on readnames and insertsize 
  comb01.5p <- bp01.5p %>%
    full_join(bp02.5p, by=c("V1"="V1", "Len"="Len")) 
  
  comb02.5p <- comb01.5p %>%
    dplyr::select(V1, V3.x, bp, Len) %>%
    dplyr::mutate(count_5p = 2)%>%
    dplyr::rename(chr_5p = V3.x,
                  pos_5p = bp) %>%
    dplyr::filter(!(is.na(chr_5p) & is.na(pos_5p)))
  
  bp01.3p <- break_p0.1.3p %>%
    dplyr::select(V1, V3, V4, V9, bp) %>%
    dplyr::rename(Len=V9)
  
  
  
  bp02.3p <- break_p0.2.3p %>%
    dplyr::select(V1, V3, V4, V9) %>%
    dplyr::rename(Len=V9) %>%
    mutate(Len = abs(Len))
  
  
  comb01.3p <- bp01.3p %>%
    full_join(bp02.3p, by=c("V1"="V1", "Len"="Len")) 
  
  comb02.3p <- comb01.3p %>%
    dplyr::select(V1, V3.x, bp, Len) %>%
    dplyr::mutate(count_3p = 2) %>%
    dplyr::rename(chr_3p = V3.x,
                  pos_3p = bp) %>%
    dplyr::filter(!(is.na(chr_3p) & is.na(pos_3p)))
  
  bp03.5p <- break_p0.3.5p %>%
    dplyr::select(V1, V3, V4, V9, bp) %>%
    dplyr::rename(Len=V9)
  
  
  
  bp04.5p <- break_p0.4.5p %>%
    dplyr::select(V1, V3, V4, V9, bp) %>%
    dplyr::rename(Len=V9) %>%
    mutate(Len = abs(Len))
  
  
  comb03.5p <- bp03.5p %>%
    full_join(bp04.5p, by=c("V1"="V1", "Len"="Len")) %>%
    dplyr::mutate(chr_5p = ifelse(!is.na(V3.x), V3.x, V3.y ),
                  pos_5p = ifelse(!is.na(bp.x), bp.x, bp.y)) %>%
    dplyr::select(V1, chr_5p, pos_5p, Len)
  
  bp03.3p <- break_p0.3.3p %>%
    dplyr::select(V1, V3, V4, V9, bp) %>%
    dplyr::rename(Len=V9)
  
  
  
  bp04.3p <- break_p0.4.3p %>%
    dplyr::select(V1, V3, V4, V9) %>%
    dplyr::rename(Len=V9) %>%
    mutate(Len = abs(Len))
  
  
  comb03.3p <- bp03.3p %>%
    full_join(bp04.3p, by=c("V1"="V1", "Len"="Len")) %>%
    dplyr::mutate(chr_3p = ifelse(!is.na(V3.x), V3.x, V3.y ),
                  pos_3p = ifelse(!is.na(V4.x), V4.x, V4.y)) %>%
    dplyr::select(V1, chr_3p, pos_3p, Len, bp)
  
  comb04_all <- comb03.5p %>%
    dplyr::full_join(comb03.3p, by=c("V1"="V1")) %>%
    mutate(Len = round(rowMeans(dplyr::select(., starts_with("Len")), na.rm = TRUE), 0)) %>%
    dplyr::select(V1, chr_5p, pos_5p, chr_3p, pos_3p, Len) %>%
    dplyr::select(V1, chr_5p, pos_5p, chr_3p, pos_3p, Len) %>%
    mutate(counts = ifelse(!is.na(chr_5p) & !is.na(chr_3p), 1, NA),
           count_5p = ifelse(!is.na(chr_5p) & is.na(chr_3p), 1,
                             ifelse(counts == 1, 1, NA)),
           count_3p = ifelse(!is.na(chr_3p) & is.na(chr_5p), 1,
                             ifelse(counts == 1, 1, NA)))
  
  # combing all the dataframes together
  comb1 <- bind_rows(comb04_all, comb02.5p, comb02.3p) 
  
  comb2 <- comb1 %>%
    group_by(chr_5p, pos_5p, chr_3p, pos_3p) %>%
    dplyr::summarise(counts_wrt_5p_3p = sum(counts, na.rm = T),
              total_5p_counts = sum(count_5p, na.rm = T),
              total_3p_counts = sum(count_3p, na.rm = T),
              avg_insert = round(mean(Len), 0)) %>%
    ungroup() %>%
    group_by(chr_5p, pos_5p) %>%
    dplyr::mutate(total_5p_counts = sum(total_5p_counts, na.rm = T)) %>%
    ungroup() %>%
    group_by(chr_3p, pos_3p) %>%
    dplyr::mutate(total_3p_counts = sum(total_3p_counts, na.rm = T)) 
  
  # filtering out all the breaks with 1 read count at 3' and 5' ends
  comb3 <- comb2 %>%
    filter(!(counts_wrt_5p_3p < 1 & total_5p_counts <=1 & total_3p_counts <= 1))
  
  # Removing previous dataframes to clear up some memory
  df <- comb3
  rm(comb2)
  rm(comb1)
  rm(comb04_all)
  rm(comb03.3p)
  rm(comb03.5p)
  rm(comb02.3p)
  rm(comb02.5p)
  rm(comb01.3p)
  rm(comb01.5p)
  rm(break_p0.4.3p)
  rm(break_p0.4.5p)
  rm(break_p0.3.3p)
  rm(break_p0.3.5p)
  rm(break_p0.2.3p)
  rm(break_p0.2.5p)
  rm(break_p0.1.3p)
  rm(break_p0.1.5p)
  rm(bp01.3p)
  rm(bp01.5p)
  rm(bp02.3p)
  rm(bp02.5p)
  rm(bp03.3p)
  rm(bp03.5p)
  rm(bp04.3p)
  rm(bp04.5p)
  
  ## Annotating the breaks with hg19
  
  #getting the starting coordinate
  myIntervals_start_5p <- df[, c(1,2)] %>%
    dplyr::rename(chr = chr_5p,
                  start = pos_5p) %>%
    mutate(end = as.numeric(start))
  
  #getting the starting coordinate
  myIntervals_start_3p <- df[, c(3,4)] %>%
    dplyr::rename(chr = chr_3p,
                  start = pos_3p) %>%
    mutate(end = as.numeric(start))
  
  # changing colnames
  colnames(myIntervals_start_5p) <- header
  
  colnames(myIntervals_start_3p) <- header
  
  # The function "makeGRangesFromDataFrame" from the library 
  # GenomicRanges makes an object GRanges called "intervals" from "myIntervals"
  interval_start_5p <- GenomicRanges::makeGRangesFromDataFrame(myIntervals_start_5p, na.rm = T)
  interval_start_3p <- GenomicRanges::makeGRangesFromDataFrame(myIntervals_start_3p, na.rm = T)
  
  # extract the list of all gene_Ids from the txdb object
  genes = genes(txdb)
  
  ################################################ Annotations for 5' break ############################################################
  
  # Use the "annotateIntervals" funtion in order to annotate 
  #the intervals with gene_Ids and produce "myAnnotation" data.frame 
  myAnnotation_start_5p <- as.data.frame(annotateIntervals(interval_start_5p, txdb))
  
  
  myDf_master_start_5p <- data.frame()
  for (i in 1 :length(myAnnotation_start_5p$gene_id)){
    if(length(c(na.omit(myAnnotation_start_5p$gene_id[i])[[1]])) != 0) {
      # annotate the interval and copy into a myDf data.frame
      myDf_start_5p <- data.frame(chr = myAnnotation_start_5p$seqnames[i], start = myAnnotation_start_5p$start[i],
                                  end = myAnnotation_start_5p$end[i], startAnno = toString(c(na.omit(myAnnotation_start_5p$gene_id[i])[[1]])))
      myDf_master_start_5p <- rbind(myDf_master_start_5p, myDf_start_5p)
    }
    if (length(c(na.omit(myAnnotation_start_5p$gene_id[i])[[1]])) == 0){
      myDF_start_0_5p <- data.frame(chr = myAnnotation_start_5p$seqnames[i], start = myAnnotation_start_5p$start[i],
                                    end = myAnnotation_start_5p$end[i], startAnno = toString(paste("intergenic")))
      myDf_master_start_5p <- rbind(myDf_master_start_5p, myDF_start_0_5p)
    }
  }
  
  newDF_start_5p <- myDf_master_start_5p[, c(1:2,4)]
  
  # Assuming 'txdb' is your TxDb object and 'df' is your dataframe
  
  # Load required packages
  library(GenomicFeatures)
  
  # Filter gene names from the 'startAnno' column
  gene_names_5p <- newDF_start_5p$startAnno[newDF_start_5p$startAnno != "intergenic"]
  
  # Initialize an empty vector to store strand information
  strand_info_5p <- character(length(gene_names_5p))
  
  # Get gene information from TxDb for non-intergenic entries
  genes_info_5p <- genes(txdb, filter = list(gene_id = gene_names_5p))
  gene_strands_5p <- strand(genes_info_5p)
  
  # Match gene strands to their respective positions in 'startAnno'
  for (i in seq_along(gene_names_5p)) {
    gene_idx_5p <- match(gene_names_5p[i], genes_info_5p$gene_id)
    if (!is.na(gene_idx_5p)) {
      strand_info_5p[i] <- gene_strands_5p[gene_idx_5p]
    } else {
      strand_info_5p[i] <- "*"
    }
  }
  
  # Update 'strand' column in the dataframe for non-intergenic entries
  newDF_start_5p$strand <- ifelse(newDF_start_5p$startAnno != "intergenic", strand_info_5p, "*")
  
  newDF_start_5p <- newDF_start_5p %>%
    dplyr::rename(breakp1 = start,
                  geneName1 = startAnno)
  
  ################################################ Annotations for 3' break ############################################################
  
  # Use the "annotateIntervals" funtion in order to annotate 
  #the intervals with gene_Ids and produce "myAnnotation" data.frame 
  myAnnotation_start_3p <- as.data.frame(annotateIntervals(interval_start_3p, txdb))
  
  
  myDf_master_start_3p <- data.frame()
  for (i in 1 :length(myAnnotation_start_3p$gene_id)){
    if(length(c(na.omit(myAnnotation_start_3p$gene_id[i])[[1]])) != 0) {
      # annotate the interval and copy into a myDf data.frame
      myDf_start_3p <- data.frame(chr = myAnnotation_start_3p$seqnames[i], start = myAnnotation_start_3p$start[i],
                                  end = myAnnotation_start_3p$end[i], startAnno = toString(c(na.omit(myAnnotation_start_3p$gene_id[i])[[1]])))
      myDf_master_start_3p <- rbind(myDf_master_start_3p, myDf_start_3p)
    }
    if (length(c(na.omit(myAnnotation_start_3p$gene_id[i])[[1]])) == 0){
      myDF_start_0_3p <- data.frame(chr = myAnnotation_start_3p$seqnames[i], start = myAnnotation_start_3p$start[i],
                                    end = myAnnotation_start_3p$end[i], startAnno = toString(paste("intergenic")))
      myDf_master_start_3p <- rbind(myDf_master_start_3p, myDF_start_0_3p)
    }
  }
  
  newDF_start_3p <- myDf_master_start_3p[, c(1:2,4)]
  
  # Assuming 'txdb' is your TxDb object and 'df' is your dataframe
  
  # Load required packages
  library(GenomicFeatures)
  
  # Filter gene names from the 'startAnno' column
  gene_names_3p <- newDF_start_3p$startAnno[newDF_start_3p$startAnno != "intergenic"]
  
  # Initialize an empty vector to store strand information
  strand_info_3p <- character(length(gene_names_3p))
  
  # Get gene information from TxDb for non-intergenic entries
  genes_info_3p <- genes(txdb, filter = list(gene_id = gene_names_3p))
  gene_strands_3p <- strand(genes_info_3p)
  
  # Match gene strands to their respective positions in 'startAnno'
  for (i in seq_along(gene_names_3p)) {
    gene_idx_3p <- match(gene_names_3p[i], genes_info_3p$gene_id)
    if (!is.na(gene_idx_3p)) {
      strand_info_3p[i] <- gene_strands_3p[gene_idx_3p]
    } else {
      strand_info_3p[i] <- "*"
    }
  }
  
  # Update 'strand' column in the dataframe for non-intergenic entries
  newDF_start_3p$strand <- ifelse(newDF_start_3p$startAnno != "intergenic", strand_info_3p, "*")
  
  newDF_start_3p <- newDF_start_3p %>%
    dplyr::rename(breakp2 = start,
                  geneName2 = startAnno,
                  chr2 = chr,
                  strand_3p = strand)
  # linking 3' and 5' dataframes together before moving on to the next step
  fin <- dplyr::left_join(df, newDF_start_5p, by=c("chr_5p"="chr", "pos_5p"="breakp1")) %>%
    dplyr::distinct(chr_5p, pos_5p, .keep_all = T) %>%
    dplyr::rename(strand_5p = strand)
  fin <- left_join(fin, newDF_start_3p, by=c("chr_3p"="chr2", "pos_3p"="breakp2")) %>%
    dplyr::distinct(chr_5p, pos_5p, chr_3p,
                    pos_3p, .keep_all = T) %>%
    dplyr::rename(donor_chr = chr_5p,
                  donor_pos = pos_5p,
                  acceptor_chr = chr_3p,
                  acceptor_pos = pos_3p)
  
  
  
  fin1 <- fin %>%
    dplyr::select(donor_chr, donor_pos, acceptor_chr, acceptor_pos, counts_wrt_5p_3p, total_5p_counts, total_3p_counts, avg_insert, strand_5p,
                  geneName1, geneName2, strand_3p) %>%
    dplyr::rename(geneName_5p = geneName1,
                  geneName_3p = geneName2)
  
  fin1$donor_chr <- ifelse(!is.na(fin1$donor_chr), paste0("chr", fin1$donor_chr), fin1$donor_chr)
  
  fin1$acceptor_chr <- ifelse(!is.na(fin1$acceptor_chr), paste0("chr", fin1$acceptor_chr), fin1$acceptor_chr)
  
  df <- fin1
  
  ####  bed files with annotations for 5p break
  dat <- df %>%
    ungroup() %>%
    dplyr::select(donor_chr, donor_pos, strand_5p, geneName_5p, counts_wrt_5p_3p, total_5p_counts, avg_insert,
                  acceptor_chr, acceptor_pos, strand_3p, geneName_3p, total_3p_counts) %>% 
    mutate(donor_chr = ifelse(is.na(donor_chr), paste("chr1"), donor_chr), 
           donor_pos = ifelse(is.na(donor_pos), paste("1"), donor_pos)) %>%
    mutate(pos_end = as.numeric(donor_pos)+1, .after = donor_pos)
  
  write.table(dat, 
              file=paste0("./bams/tmp_bams/bp/tmp/",id,"_breaks_5p.bed"), 
              quote = F,sep = "\t", row.names = F, col.names = F)
  
  #### bed files with annotations for 3p break
  
  dat <- df %>%
    ungroup() %>%
    dplyr::select(acceptor_chr, acceptor_pos, strand_3p, geneName_3p,  total_3p_counts) %>% 
    dplyr::filter(!is.na(acceptor_chr)) %>%
    mutate(pos_end = as.numeric(acceptor_pos)+1, .after = acceptor_pos)
  
  write.table(dat, 
              file=paste0("./bams/tmp_bams/bp/tmp/",id,"_breaks_3p.bed"), 
              quote = F,sep = "\t", row.names = F, col.names = F)
  
}

