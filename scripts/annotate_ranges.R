#!/usr/bin/env Rscript
#
# Invoke with Rscript --vanilla script.R input output (hg_version)
# Example: Rscript --vanilla annotate_ranges.R SRX1760583.clusters.csv SRX1760583.clusters.annotated.csv hg19
#
# Script takes the output of PARalyzer showing read clusters and their genomic
# coordinates computed by bowtie. We then use AnnotationHub's databases for the human
# genome to assign Entrez Gene ID, Gene Symbols, the gene transcripts that they were
# matched to. Finally, The position in each gene transcript that the clusters mapped to
# are used to test if the gene coordinates are located in exon/intron/cds/5'UTR/3'UTR.
# The final output are returned as a ratio of hits to # of transcripts that a given
# cluster hit.
#

library(tibble)
library(tidyr)
library(dplyr)
library(readr)
library(purrr)
library(GenomicRanges)
library(GenomicFeatures)
library(AnnotationHub)

args = commandArgs(trailingOnly=TRUE)

input <- args[1]
output <- args[2]
genomeVersion <- args[3]

#input <- "SRR3502977.clusters.csv"
#output <- "SRR3502977.clusters.annotated.alt.csv"
#genomeVersion <- "hg19"


#Checking arguments
if (length(args) < 2) {
  stop("Input and/or output args missing", call.=FALSE)}
if (is.na(genomeVersion)) {
  genomeVersion <- "hg19"
} 
if (!genomeVersion %in% c("hg18", "hg19", "hg38")) {
  stop("genomeVersion is wrong, must be hg18, hg19 or hg38")
}

#Read file and make temp index.
df <- read_csv(file = input)
df <- rowid_to_column(df, "rowIndex")
# Detect columns:
r_chr <- base::colnames(df) %>% grep(pattern = "Chromosome")
r_start <- base::colnames(df) %>% grep(pattern = "Start")
r_end <- base::colnames(df) %>% grep(pattern = "End")
r_ID <- base::colnames(df) %>% grep(pattern = "ID")
r_str <- base::colnames(df) %>% grep(pattern = "Strand")
#Sanity check:
if (c(r_start, r_end, r_ID) %>% length != 3) stop("Problem with column names.")
#Construct GRanges object:
gr <- 
  GRanges(seqnames = dplyr::select(df, r_chr) %>% pull %>% Rle,
          ranges = IRanges(start = dplyr::select(df, r_start) %>% pull,
                           end = dplyr::select(df, r_end) %>% pull,
                           names = dplyr::select(df, r_ID) %>% pull),
          strand = dplyr::select(df, r_str) %>% pull %>% Rle)

# Load annotation.
hub <- AnnotationHub()
# Note: OrgDb is independent of human genome version.
OrgDb_An <- 
  with(hub, 
       subset(hub, species == "Homo sapiens" 
              & rdataclass == "OrgDb")) %>% .[[1]]
TxDb_An <-
  with(hub, 
       subset(hub, 
              hub$species == "Homo sapiens" 
              & rdataclass == "TxDb" 
              & paste0(genomeVersion, ".knownGene") %>% 
                grepl(., hub$title))
       ) %>% 
  .[[1]]

# Produce annotation by gene ID:
tx_genes <- genes(TxDb_An)
ov <- findOverlaps(gr, tx_genes)
#Get the gene IDs and symbols:
gn <- elementMetadata(tx_genes)[subjectHits(ov),"gene_id"]
sb <- AnnotationDbi::select(OrgDb_An, unlist(gn), "SYMBOL", "ENTREZID") %>% 
  as_tibble() %>%
  add_column("rowIndex" = queryHits(ov))
#Sanity check to ensure that the queryHits are unique
if (any(duplicated(queryHits(ov)))) {
  #print("WARNING - non-unique hits on gene")
  sb <- 
    sb %>%
    group_by(rowIndex) %>%
    summarize(ENTREZID = paste0(ENTREZID, collapse = ";"),
              SYMBOL = paste0(SYMBOL, collapse = ";"))
}
df_upd <- 
  left_join(df, sb, by = "rowIndex") %>%
  dplyr::rename(Symbol = SYMBOL, EntrezID = ENTREZID) %>%
  replace_na(list(Symbol = "Intergenic"))

#Find annotation by transcript:
tx <- transcripts(TxDb_An)
ov_tx <- findOverlaps(gr, tx)
#Make df to hold overlap hits, mapping between query index and database transcripts:
ovdf <- 
  tibble(queryHits = queryHits(ov_tx), subjectHits = subjectHits(ov_tx)) %>%
  tibble::add_column(tx_name = elementMetadata(tx)[subjectHits(ov_tx), "tx_name"])

# Now obtain the metadata based on transcript:
annotationHitsByTranscript <- function(type,
                                       transcript_overlap = ov_tx,
                                       database = TxDb_An,
                                       query_gr = gr) {
  #Transcript overlap computed outside of function to only calculate it once.
  #Function checks if there is an overlap in the range specified by the query (input)
  #and the transcript which it overlaps.
  #Essentially, we know that they overlap (due to findOverlaps), but we don't know
  #if the overlap is in a CDS, exon, intron, UTR etc.
  #Returns a vector of ratios of overlap as fractions: If 1 of 3 transcripts overlapped,
  #the float would be 0.33
  if (type == "cds") tx_target <- cdsBy(database, by = c("tx")) else 
    if (type == "exons") tx_target <- exonsBy(database, by = c("tx")) else
      if (type == "introns") tx_target <- intronsByTranscript(database) else
        if (type == "fiveUTRs") tx_target <- fiveUTRsByTranscript(database) else
          if (type == "threeUTRs") tx_target <- threeUTRsByTranscript(database)
  ov_target <- 
    map2_int(.x = queryHits(transcript_overlap),
       .y = subjectHits(transcript_overlap),
       .f = function(x, y) {
         if (as.character(y) %in% names(tx_target)) {
           findOverlaps(query_gr[x], tx_target[as.character(y)]) %>% 
             queryHits %>% 
             length %>% 
             as.integer() %>%
             return()
         } else {
           return(as.integer(0))}
       })
  ov_target %>% return()
}

ovdf_an <- add_column(ovdf, 
                      cds = annotationHitsByTranscript(type = "cds"),
                      exons = annotationHitsByTranscript(type = "exons"),
                      introns = annotationHitsByTranscript(type = "introns"),
                      fiveUTRs = annotationHitsByTranscript(type = "fiveUTRs"),
                      threeUTRs = annotationHitsByTranscript(type = "threeUTRs"))
ovdf_sum <- 
  ovdf_an %>%
  group_by(queryHits) %>%
  summarize(exons = mean(exons),
            introns = mean(introns),
            cds = mean(cds),
            fiveUTRs = mean(fiveUTRs),
            threeUTRs = mean(threeUTRs),
            transcript_num = n(),
            "transcripts" = paste0(tx_name, collapse=";"))
df_an <-  
  left_join(df_upd, ovdf_sum, by = c("rowIndex" = "queryHits")) %>% 
  dplyr::select(-rowIndex)

#Could rearrange dataframe here to bring interesting results more forward, but that also
#makes the script specific to only one output. Without it, it can work on any dataframe
#suitable for genomic ranges.

#Output the df.
write_csv(df_an, path = output, col_names = TRUE)


