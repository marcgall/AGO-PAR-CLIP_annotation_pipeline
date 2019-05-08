#!/usr/bin/env Rscript
#
# Invoke with Rscript --vanilla script.R input output (hg_version) (threads)
# Example: Rscript --vanilla annotate_ranges.R SRX1760583.clusters.csv SRX1760583.clusters.annotated.csv hg19 5
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
library(dplyr)
library(readr)
library(purrr)
library(GenomicRanges)
library(GenomicFeatures)
library(AnnotationHub)
library(foreach)
library(doParallel)

args = commandArgs(trailingOnly=TRUE)

input <- args[1]  #"SRX1760583.clusters.csv"
output <- args[2] # "SRX1760583.clusters.annotated.csv"
genomeVersion <- args[3] # "hg19"
threads <- args[4] %>% as.integer

#input <- "SRX1760583.clusters.csv"
#output <- "SRX1760583.clusters.annotated.csv"
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
if (is.na(threads)) {
  threads <- 1
} else if (threads > 5) { #can't use more than 5
  threads <- 5
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
  rowid_to_column("rowIndex")
df_upd <- 
  left_join(df, sb, by = "rowIndex") %>%
  dplyr::rename(Symbol = SYMBOL, EntrezID = ENTREZID)

#Find annotation by transcript:
tx <- transcripts(TxDb_An)
ov_tx <- findOverlaps(gr, tx)
#Make df to hold overlap hits, mapping between query index and database transcripts:
ovdf <- 
  tibble(queryHits = queryHits(ov_tx), subjectHits = subjectHits(ov_tx)) %>%
  add_column(tx_name = elementMetadata(tx)[subjectHits(ov_tx), "tx_name"])

# Now obtain the metadata based on transcript:
annotationHitsByTranscript_MC <- function(type, input_file) {
  # NOTE: For MC version, database and transcript overlap has to be recomputed
  # in this namespace.
  #
  #Transcript overlap computed outside of function to only calculate it once.
  #Function checks if there is an overlap in the range specified by the query (input)
  #and the transcript which it overlaps.
  #Essentially, we know that they overlap (due to findOverlaps), but we don't know
  #if the overlap is in a CDS, exon, intron, UTR etc.
  #Returns a vector of ratios of overlap as fractions: If 1 of 3 transcripts overlapped,
  #the float would be 0.33
  #
  # Load dataframe for GR construction
  df <- read_csv(file = input_file)
  df <- rowid_to_column(df, "rowIndex")
  # Detect columns:
  r_chr <- base::colnames(df) %>% grep(pattern = "Chromosome")
  r_start <- base::colnames(df) %>% grep(pattern = "Start")
  r_end <- base::colnames(df) %>% grep(pattern = "End")
  r_ID <- base::colnames(df) %>% grep(pattern = "ID")
  r_str <- base::colnames(df) %>% grep(pattern = "Strand")
  # Construct GR for comparions:
  gr <- 
    GRanges(seqnames = dplyr::select(df, r_chr) %>% pull %>% Rle,
            ranges = IRanges(start = dplyr::select(df, r_start) %>% pull,
                             end = dplyr::select(df, r_end) %>% pull,
                             names = dplyr::select(df, r_ID) %>% pull),
            strand = dplyr::select(df, r_str) %>% pull %>% Rle)
  #Loading database:
  hub <- AnnotationHub()
  database <-
    with(hub, 
         subset(hub, 
                hub$species == "Homo sapiens" 
                & rdataclass == "TxDb" 
                & paste0(genomeVersion, ".knownGene") %>% 
                  grepl(., hub$title))
    ) %>% 
    .[[1]]
  tx <- transcripts(database)
  transcript_overlap <- findOverlaps(gr, tx)
  #
  # Find annotation hits:
  print(type)
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
           findOverlaps(gr[x], tx_target[as.character(y)]) %>% 
             queryHits %>% 
             length %>% 
             as.integer() %>%
             return()
         } else {
           return(as.integer(0))}
       })
  ov_target %>% return()
}

#Multiprocessing
cl <- makeCluster(threads)
registerDoParallel(cl)

annotation_types <- c("cds", "exon", "introns", "fiveUTRs", "threeUTRs")

annotated_cols <- 
  foreach(n = 1:5, .packages = c("tibble", "dplyr", "readr", "purrr", "GenomicRanges", 
                                 "GenomicFeatures", "AnnotationHub")) %dopar% 
  annotationHitsByTranscript_MC(type = annotation_types[n], 
                                input_file = input)

ovdf_an <- add_column(ovdf, 
                      cds = annotated_cols[[1]],
                      exons = annotated_cols[[2]],
                      introns = annotated_cols[[3]],
                      fiveUTRs = annotated_cols[[4]],
                      threeUTRs = annotated_cols[[5]])
stopCluster(cl)

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


