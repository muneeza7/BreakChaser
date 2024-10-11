rm(list=ls())

## load library
library(pacman)

p_load(
  "stringr","tidyverse", 
  "Biostrings",
  "BSgenome.Hsapiens.UCSC.hg19",
  "dplyr", 
  "stringr",
  "annotate",
  "org.Hs.eg.db",
  "TxDb.Hsapiens.UCSC.hg19.knownGene",
  "GenomicRanges",
  "Homo.sapiens",
  "GenomicFeatures",
  "GenomicAlignments",
  "SummarizedExperiment",
  "MatrixGenerics",
  "biomaRt", "memes")
pacman::p_load(readxl, tidyverse, xlsx, data.table, dplyr, purrr, matrixStats, ggplot2, openxlsx, magrittr, gapminder,GenomicRanges)
options(scipen = 999)
################################# Functions ##############################

## gene of interest
goi <- c( "ASXL1", "ATM","BCOR","BCORL1","BTG1","CALR","CDKN2A","CDKN2B","CUX1","DNMT3A","ETV6",
          "EZH2","GATA2","IKZF1","JAK2","KDM6A","KIT","KMT2D","PAX5","PHF6","PPM1D","RB1",
          "RUNX1","SETD1B","SETD2","SF3B1","SMC3","SRSF2","STAG2","TET2","TP53","UBE2A","WT1", "MTAP", "CDKN2A-AS1","CDKN2B-AS1", 
          "EP300", "NUP98", "CREBBP", "RAD21", "RCBTB2", "CBFB", "MYH11")


## Reading files with soft-clipped reads

files_5p <- list.files(path="./bams/tmp_bams/bp/tmp/", 
                       pattern="breaks_5p_mapped.bed", full.names = T)

files_3p <- list.files(path="./bams/tmp_bams/bp/tmp/", 
                       pattern="3p_mapped.bed", full.names = T)

### recovered intergenic breakpoints that got filtered out from mapped.bed files
files_inter_5p <- list.files(path="./bams/tmp_bams/bp/tmp/", 
                             pattern="intergenic_5p_mapped.bed", full.names = T)

for (k in 1:length(files_5p)) {
  
  # Annotations for intron, exon positions at 5' end
  dat <- read.table(files_5p[k], sep = "\t", header = F, stringsAsFactors =  F)
  
  ## extract id
  id1 <- str_extract(files_5p[k], "[a-zA-Z0-9-_\\.]+_breaks_5p_mapped.bed")
  id <- str_trim(gsub("_breaks_5p_mapped.bed|X__", "", id1))
  
  print(paste("processing ...  ", id1, sep=""))
  
  postdat.1 <- dat %>%
    mutate(V4 = V17) %>%
    mutate(V4=ifelse(V4=="+" | V4=="-", V4, V17)) %>%
    arrange(V1, V2, V3, V4, V22, V14, V15, V21) %>%
    mutate(overlap=V23-V24) %>%
    dplyr::filter(! (V19=="intron" & overlap==0)) %>% 
    group_by(V1, V2, V3, V18) %>%
    dplyr::mutate(
      first5 = ifelse(V17=="+", dplyr::first(V15), dplyr::last(V16) ), ## include strand ifelse
      last5 = ifelse(V17=="+", dplyr::first(V16), dplyr::last(V15)),
      first3 = ifelse(V17=="+", dplyr::last(V15), dplyr::first(V16)),
      last3 = ifelse(V17=="+", dplyr::last(V16), dplyr::first(V15)),
      event5 = ifelse(V17=="+", dplyr::first(V19), dplyr::last(V19)), 
      event3 = ifelse(V17=="+", dplyr::last(V19), dplyr::first(V19)),
      location5 = ifelse(V17=="+", dplyr::first(V21), dplyr::last(V21)),
      location3 = ifelse(V17=="+", dplyr::last(V21), dplyr::first(V21)),
      strandDirection= ifelse(V17==V4, "same", "opposite")) %>%
    mutate(event=paste(event5, event3, sep="-"),
           nExSkip=ifelse(abs(location5-location3)==0, abs(location5-location3), abs(abs(location5-location3)-1)),
           exon=paste(event5, location5, event3, location3, sep="")) %>% ## collapse exon number per transcript
    mutate(exon=gsub("exon", "E", exon),
           exon=gsub("intron", "I", exon),
           withinEx5=ifelse(V17=="+" & V2>first5 & V2 < last5, "Yes", 
                            ifelse(V17=="-" & V3<first5 & V3>last5, "Yes", "No")),
           withinEx3=ifelse(V17=="+" & V3>first3 & V3 < last3, "Yes", 
                            ifelse(V17=="-" & V2<first3 & V2>last3, "Yes", "No"))) %>%
    ungroup() %>%
    group_by(V1, V2, V3) %>%
    mutate(nGene = length(unique(V22)),
           nTranscripts=length(unique(V18)), ## have not take into account for nGene > 1
           distance_from_last5= ifelse(V4=="+", V2-last5, V3-last5),
           distance_from_first3= ifelse(V4=="+", V3-first3, V2-first3)) %>%
    mutate(total=paste(event, collapse = "|"),
           Gene=paste(unique(V18), collapse="|")) %>%
    mutate(abnormality=paste(ifelse(V4=="+" & str_detect(total, "exon-exon") & abs(distance_from_last5)==1 & 
                                      abs(distance_from_first3)==2 & nExSkip<=1 & nGene==1 & strandDirection=="same", "normal",
                                    ifelse(V4=="-" & str_detect(total, "exon-exon") & 
                                             abs(distance_from_last5)==2 & abs(distance_from_first3)==1 & 
                                             nExSkip<=1 & nGene==1 & strandDirection=="same", "normal", "abnormal")), 
                             collapse="|")) %>%
    group_by(V1, V2, V3) %>%
    arrange(desc(V20)) %>% 
    ungroup() %>%
    distinct(V1, V2, V9, V10, .keep_all = T) %>%  
    ungroup() %>%
    dplyr::select(-V3, -V14, -V15, -V16, -V17, -V18, -V20, -V21, -V22, -V23, -V24, -overlap,
                  -first5, -last5, -first3, -last3, -event5, -event3, -strandDirection, -nExSkip, -exon, -withinEx5,
                  -withinEx3, -nGene, -nTranscripts, -total, -Gene, -abnormality, -location3) %>%
    dplyr::rename(donor_chr = V1, 
                  donor_pos = V2, 
                  strand_5p = V4, 
                  geneName_5p = V5, 
                  counts_wrt_5p_3p = V6, 
                  total_5p_counts = V7, 
                  avg_insert = V8,
                  acceptor_chr = V9,
                  acceptor_pos = V10, 
                  strand_3p = V11,
                  geneName_3p = V12, 
                  total_3p_counts = V13, 
                  location_5p = V19,
                  locNum_5p = location5,
                  lastExDist_5p = distance_from_last5,
                  nextExDist_5p = distance_from_first3) %>%
    as.data.frame()
  
  # appending the list of recovered breaks to the list of 5' breaks
  dat0 <- read.table(files_inter_5p[k], sep = "\t", header = F, stringsAsFactors =  F) %>%
    dplyr::select(V1, V2, V4,V5, V6, V7, V8, V9, V10,  V11, V12, V13) %>%
    dplyr::rename(donor_chr = V1,
                  donor_pos = V2,
                  strand_5p = V4, 
                  geneName_5p = V5,
                  counts_wrt_5p_3p = V6,
                  total_5p_counts = V7,
                  avg_insert = V8,
                  acceptor_chr = V9,
                  acceptor_pos = V10,
                  strand_3p = V11,
                  geneName_3p = V12,
                  total_3p_counts =V13) %>%
    dplyr::mutate(donor_chr = ifelse((donor_chr == "chr1" & donor_pos == "1"), NA, donor_chr),
                  donor_pos = ifelse((donor_chr == "chr1" & donor_pos == "1"), NA, donor_pos))
  postdat.1.5p <- bind_rows(postdat.1, dat0)
  
  # Annotations for intron, exon positions at 3' end
  dat1 <- read.table(files_3p[k], sep = "\t", header = F, stringsAsFactors =  F)
  
  postdat.1.3p <- dat1 %>%
    mutate(V4.1 = V4) %>%
    mutate(V4 = V10) %>%
    mutate(V4=ifelse(V4=="+" | V4=="-", V4, V10)) %>%
    arrange(V1, V2, V3, V4, V15, V7, V8, V14) %>% #(V1, V2, V3, V4, V18, V11, V12, V17)
    mutate(overlap=V16-V17) %>%
    dplyr::filter(! (V12=="intron" & overlap==0)) %>%
    group_by(V1, V2, V3, V11) %>%
    dplyr::mutate(
      first5 = ifelse(V10=="+", dplyr::first(V8), dplyr::last(V9) ), ## include strand ifelse
      last5 = ifelse(V10=="+", dplyr::first(V9), dplyr::last(V8)),
      first3 = ifelse(V10=="+", dplyr::last(V8), dplyr::first(V9)),
      last3 = ifelse(V10=="+", dplyr::last(V9), dplyr::first(V8)),
      event5 = ifelse(V10=="+", dplyr::first(V12), dplyr::last(V12)),
      event3 = ifelse(V10=="+", dplyr::last(V12), dplyr::first(V12)),
      location5 = ifelse(V10=="+", dplyr::first(V14), dplyr::last(V14)),
      location3 = ifelse(V10=="+", dplyr::last(V14), dplyr::first(V14)),
      strandDirection= ifelse(V10==V4, "same", "opposite")) %>%
    mutate(event=paste(event5, event3, sep="-"),
           nExSkip=ifelse(abs(location5-location3)==0, abs(location5-location3), abs(abs(location5-location3)-1)),
           exon=paste(event5, location5, event3, location3, sep="")) %>% ## collapse exon number per transcript
    mutate(exon=gsub("exon", "E", exon),
           exon=gsub("intron", "I", exon),
           withinEx5=ifelse(V10=="+" & V2>first5 & V2 < last5, "Yes",
                            ifelse(V10=="-" & V3<first5 & V3>last5, "Yes", "No")),
           withinEx3=ifelse(V10=="+" & V3>first3 & V3 < last3, "Yes",
                            ifelse(V10=="-" & V2<first3 & V2>last3, "Yes", "No"))) %>%
    ungroup() %>%
    group_by(V1, V2, V3) %>%
    mutate(nGene = length(unique(V15)),
           nTranscripts=length(unique(V11)),
           distance_from_last5= ifelse(V4=="+", V2-last5, V3-last5),
           distance_from_first3= ifelse(V4=="+", V3-first3, V2-first3)) %>%
    mutate(total=paste(event, collapse = "|"),
           Gene=paste(unique(V11), collapse="|")) %>%
    mutate(abnormality=paste(ifelse(V4=="+" & str_detect(total, "exon-exon") & abs(distance_from_last5)==1 &
                                      abs(distance_from_first3)==2 & nExSkip<=1 & nGene==1 & strandDirection=="same", "normal",
                                    ifelse(V4=="-" & str_detect(total, "exon-exon") &
                                             abs(distance_from_last5)==2 & abs(distance_from_first3)==1 &
                                             nExSkip<=1 & nGene==1 & strandDirection=="same", "normal", "abnormal")),
                             collapse="|")) %>%
    group_by(V1, V2, V3) %>%
    dplyr::arrange(desc(V13)) %>%
    ungroup() %>%
    distinct(V1, V2, V4.1, .keep_all = T) %>%
    ungroup() %>%
    dplyr::select(-V3, -V7, -V8, -V9, -V10, -V11, -V13, -V14, -V15, -V16, -V17, -overlap,
                  -first5, -last5, -first3, -last3, -event5, -event3, -strandDirection, -nExSkip, -exon, -withinEx5, -withinEx3, -nGene, -nTranscripts, -total, -Gene, -abnormality, -location3) %>%
    dplyr::rename(acceptor_chr = V1,
                  acceptor_pos = V2,
                  strand_3p = V4,
                  geneName_3p = V5,
                  total_3p_counts = V6,
                  location_3p = V12,
                  locNum_3p = location5,
                  lastExDist_3p = distance_from_last5,
                  nextExDist_3p = distance_from_first3) %>%
    as.data.frame()
  
  # joing 5' and 3' end together
  postdat <- left_join(postdat.1.5p, postdat.1.3p, by=c("acceptor_chr"="acceptor_chr", "acceptor_pos"="acceptor_pos"
                                                        , "geneName_3p"="geneName_3p"))
  postdat <- postdat %>%
    dplyr::mutate(strand_3p = strand_3p.y, .after = acceptor_pos) %>%
    dplyr::rename(total_3p_counts = total_3p_counts.x,
                  event_5p = event.x,
                  event_3p = event.y) %>%
    dplyr::select(donor_chr, donor_pos, strand_5p, geneName_5p, 
                  location_5p, locNum_5p, event_5p, lastExDist_5p, nextExDist_5p, counts_wrt_5p_3p,
                  total_5p_counts, avg_insert, acceptor_chr, acceptor_pos, strand_3p, geneName_3p,
                  location_3p, locNum_3p, event_3p, lastExDist_3p, nextExDist_3p, total_3p_counts) %>%
    as.data.frame()
  
  # filtering for the events according to the genes of interest, region of interest, and location of the break and deletion size
  new <- postdat %>%
    dplyr::filter(geneName_5p %in% goi | geneName_5p %in% "intergenic" | geneName_3p %in% goi |
                    geneName_3p %in% "intergenic" | is.na(geneName_3p)) %>%
    dplyr::filter(!event_5p %in% c("intron-exon", "exon-intron") | event_5p %in% "intergenic") %>%
    dplyr::filter(abs(lastExDist_5p) > 10 & abs(nextExDist_5p) > 10 | is.na(lastExDist_5p) | is.na(nextExDist_5p)) %>%
    dplyr::filter(abs(lastExDist_3p) > 10 & abs(nextExDist_3p) > 10 |
                    is.na(lastExDist_3p) == T & is.na(nextExDist_3p) == T) %>%
    dplyr::distinct(donor_chr, donor_pos, acceptor_chr, acceptor_pos, .keep_all = T) %>%
    dplyr::filter(donor_chr %in% acceptor_chr | is.na(acceptor_chr)) %>%
    dplyr::mutate(DelSize = acceptor_pos - donor_pos,
                  chr = ifelse((donor_chr %in% acceptor_chr) & !(is.na(donor_chr)), donor_chr, acceptor_chr)) %>%
    dplyr::select(chr, donor_pos, acceptor_pos, geneName_5p, geneName_3p, DelSize, avg_insert, counts_wrt_5p_3p, total_5p_counts, total_3p_counts,
                  strand_5p, strand_3p, location_5p, locNum_5p, location_3p, locNum_3p, lastExDist_5p,
                  nextExDist_5p, lastExDist_3p, nextExDist_3p) %>%
    dplyr::filter(abs(DelSize) > 10000 | is.na(DelSize) == T) %>%
    as.data.frame()
  
  # Writing into csv files
  write.csv(new, 
            file=paste("./bams/tmp_bams/bp/",id, "_deletion_breaks.csv"), 
            row.names = F, na="")
}

