#
# Add genomic information to PARalyzer output files.
#
library(tidyverse)
library(GenomicRanges)
library(AnnotationHub)
#library(Organism.dplyr)

setwd("Documents/PAR-CLIP/snakemake_pipeline/PARalyzer/")

df <- read_csv(file = "SRX1760583.clusters.csv")

# the info in csv closely matches the bed file format, which can be imported as a Granges object.
# Indeed, it pretty much IS the bed file with a few additions.


gr <- with(df, 
           GRanges(seqnames = Rle(Chromosome), #Rle converts the vector to a more efficient version
                   ranges = IRanges(start = ClusterStart, end = ClusterEnd, names = ClusterID),
                   strand = Rle(Strand), 
                   clusterSequence = ClusterSequence))

##################### Must be able to get this information from annotation:
##################### BUT is it necessary at this stage?
genome(gr) <- "hg19"
seqlengths(gr)

#check that all ClusterIDs are unique:
test <- df$ClusterID %>% duplicated
sum(test)
#They are unique.
  
genome <- BSgenome.Hsapiens.UCSC.hg19
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

transcriptsBy(txdb, by = "gene")

seqinfo(txdb)

#Finding AR by Entrez Gene ID
genes(txdb, filter = list(gene_id = 367))
exons(txdb, filter = list(gene_id = 367))
transcripts(txdb, filter = list(gene_id = 367))
genes(txdb, filter = list(gene_id = 367), columns = c("gene_id"))

####### Get annotation info from annotationhub:
#Get annotaiton for GRCh37 / hg19
hub <- AnnotationHub()
hs <- query(hub, "Hsapiens")
query(hub, "Hg19")

subset(hub, hub$species == "Homo sapiens")
subset(hub, hub$species == "Homo sapiens" & hub$rdataclass == "GRanges" & hub$dataprovider == "Ensembl")
subset(hub, hub$species == "Homo sapiens" & hub$rdataclass == "OrgDb")
OrgDb.An <- subset(hub, hub$species == "Homo sapiens" & hub$rdataclass == "OrgDb") %>% .[[1]]
#Here's how we can retrieve annotation from OrgDb:
select(OrgDb.An, keys = "AR", keytype = "SYMBOL", columns = c("GENENAME", "ENSEMBL", "ENTREZID"))


subset(hub, hub$species == "Homo sapiens")
# These do the same, but with grepl I can perform partial matches too.
subset(hub, grepl("Homo sapiens", hub$species))

#Get TxDb through annotationhub:
subset(hub, hub$species == "Homo sapiens" & grepl("TxDb", hub$title))
TxDb.An <- query(hub, "TxDb.Hsapiens.UCSC.hg19.knownGene.sqlite") %>% .[[1]]

### Now that I have the TxDb, I should be able to add Entrez Gene ID and genetic information to my PARalyzer results.
# Using e.g 'genes()' on txdb returns a GRanges object.
# Maybe this function? mapIdsToRanges
Txgenes <- genes(txdb)
findOverlaps(gr, Txgenes)
countOverlaps(gr, Txgenes) %>% table()

#This function could possibly be used to set metadata (like 'genes' if I can 
# extract the index of the output and apply it to the original gr dataframe.)
test <- subsetByOverlaps(gr, Txgenes)
# Could get the names as my index.
# Like so:
gr$type <- NA
gr[which(names(gr) %in% names(test)),] #This DIDNT work to set the values.

#move on and deal with it later.

#Use Organism.dplyr to combine TxDb and OrgDb.
# NOTE! This doesn't work with sequences retrieved from AnnotationHub. 
# It requires the package.
supportedOrganisms() %>% print(n=Inf)

library(org.Hs.eg.db)
src <- src_organism("TxDb.Hsapiens.UCSC.hg19.knownGene")

columns(src)
# Entered this info, but I'm unsure of how I can use it currently. I still
# need to key back to genomic range information on location.
columns(OrgDb.An)
select(src, keys="AR", keytype="symbol", 
       columns = c("cds_chrom", "cds_end", "cds_id", "ensembl", "entrez", "exon_start", "tx_name"))

#Use txdb to link entrez ID to gr object.
#
# My goal:
# Get Entrez ID by RANGE

#This works, BUT uses the Organism.dplyr which I might want to avoid.
genes(src, columns = "entrez")

#This returns the mapping that I want, BUT: it's not complete
gr.id <- mapRangesToIds(x = TxDb.An, ranges = gr, type = gene)

#What about the mappings that did not return an ID?
# I suppose that the region doen't correspond to a gene at all.
# This is similar to:
subsetByOverlaps(gr, Txgenes)
#but the latter does not add gene IDs at all. Method does not appear capable
# of conserving metadata.
# Findoverlaps instead gives me a dataframe (in another format) with the hits.
# Using this, I may be able to subset the ranges and thus preserve metadata.
overlaps <- findOverlaps(gr, Txgenes)
overlaps %>% subjectHits
overlaps %>% as_tibble

# This adds the metadata, but only to the overlap.
test <- subsetByOverlaps(gr, Txgenes)
idx <- overlaps %>% subjectHits
values <- DataFrame(gene_id=Txgenes$gene_id[idx])
mcols(test) <- c(mcols(test), values)
#Error here because there are duplicate values - but using 'unique' removes too much.

# Try locateVariants() from VariantAnnotation
library(VariantAnnotation)

locvar <- locateVariants(query = gr, subject = TxDb.An, AllVariants())
# I wonder if there's a bug in the program, since it returns that 1171 of the sequences
# are out of bound. On the other hand, this appears to fit the number that didn't map to
# a gene. But not quite.

# This locvar has more than twice as many rows as the starting gr object. Why?
# May be due to the range falling in more than one region.
# If I turn it to a tibble, the clusterId name is dropped.

locvar$GENEID %>% unique() %>% length()
#Something wrong with the conversion to tibble, filtering and selecting doesn't work.
locvar$LOCATION %>% unique

names(locvar)
#I don't understand this output, many of them here have NO name attribute. 
# Could these be the ranges that were out of bounds?
locvar[5:15,]
# Checking the first few ranges [5:15] reveals that the unnamed ranges are the
# duplicates of the next range, but labeled with a different location. So far so good.
#
# But further down the line, bigger gaps appear. Are these still all referring to the
# same range?
locvar[15:25,]
#Yes, they do refer to duplicate ranges. In this case, a promoter assigned to
# several TXIDs and GENEIDs.


table(locvar$LOCATION)

locateVariants(query = gr, subject = TxDb.An, CodingVariants())


#Perhaps I can aggregate somehow?

#Let's say I use locateVariants: How do I map these back to the original gr?
# I need to match it to either clusterID or to the ranges: These are the only
# unique identifiers.


#Try GRannot:
library(REMP)
# Note: This packages requires SO much loading, I don't think it's justifiable.
refseq <- 
  subset(hub, 
         hub$species == "Homo sapiens" & 
           grepl("RefSeq", hub$title) & 
           hub$genome == "hg19")[[1]]
GRannot(gr, refseq, symbol=TRUE)
# In addition to requiring WAY too many dependencies - the refseq data 
# does not provide interesting annotation to me.
# And in the end, it simply returned an nonsense error 
#(saying the input gr is not a GRange object)


##############################
############################## New day, new plan: Now use overlaps to produce cols
############################## of hit/no hit.
##############################


Txgenes <- genes(TxDb.An)

ov <- findOverlaps(gr, Txgenes)
#Can I use these two to filter the vectors?
queryHits(ov)
subjectHits(ov)

gr[queryHits(ov)]
# Yes, this works perfectly.
gr[queryHits(ov),]
# Doesn't appear to matter if I include the ',' or not to specify rows.
# Probably because slicing by index doesn't make as much sense on cols of the GR object

Txgenes[subjectHits(ov)]
###
subsetByOverlaps(gr, Txgenes)
subsetByOverlaps(Txgenes, gr)
### This produces overlap of 1713 or 1524 ranges. Why the difference?
### And the findOverlap hits gave 2099 ranges.
### Can this be influenced by the type of overlap? Probably not.
### It must have to do with single ranges mapping to more than one range
### in the other file.
###
### ov this by checking for duplicate entries in the findOverlaps hits.
queryHits(ov) %>% duplicated() %>% sum(.) #386
subjectHits(ov) %>% duplicated() %>% sum(.) #575
### Yes, this is the case. How do we deal with this?
### These add up to the same total overlap as the 2099 ranges.
###
### It is partly to be expected, because getting genes() doesn't take exons, etc into account
### And because we have duplicates in the input file.
genes(TxDb.An)
transcripts(TxDb.An)
transcriptsBy(TxDb.An, by = "gene")


ov <- findOverlaps(gr, Txgenes)
Txgenes <- genes(TxDb.An)
#Get the gene IDs
gn <- elementMetadata(Txgenes)[subjectHits(ov),"gene_id"]
sb <- select(OrgDb.An, unlist(gn), "SYMBOL", "ENTREZID")
#This works! Getting the symbol directly from the entrez id. More info from OrgDb can be queried.
#
# now, can I add it correctly back to the dataframe with correct spaces?
#
# Make a new column in original dataframe and access it by the index from findOverlaps()
df <- 
  df %>% add_column(Symbol = NA, Entrez_ID = NA, .before = "Chromosome")
df[queryHits(ov), "Symbol"] <- sb[,"SYMBOL"]
df[queryHits(ov), "Entrez_ID"] <- sb[,"ENTREZID"]
#
# Now add info on region:
# However, slight problem here: The info on UTRs is confined to extract by transcript.
# So there, I HAVE to get transcript information first. This means dealing with duplicate
# entries.
# How can I do this?

tx <- transcripts(TxDb.An)
ov <- findOverlaps(gr, tx)

#We get 5891 hits, twice the query-length.
queryHits(ov) %>% unique %>% length
#2029 of them are unique, leaving ~1000 which hits multiple transcripts accounting for
# the remaining ~4000.
#
#
subsetByOverlaps(gr, tx)

elementMetadata(tx)[subjectHits(ov), "tx_id"]

##### Can use 'select' to get index of the first overlapping element in each element in query.
findOverlaps(gr, tx, select = "first")
findOverlaps(gr, tx)

### Build dataframe out of the overlaps key, adding transcript ID per column and adding more
### information here.

ovdf <- tibble(queryHits = queryHits(ov), subjectHits = subjectHits(ov))
ovdf$tx_name <- elementMetadata(tx)[subjectHits(ov), "tx_name"]
#   (It turns out that the tx name is immaterial: Only the numeric ID is used internally)
# Now obtain the metadata based on transcript:

tx_cds <- cdsBy(TxDb.An, by = c("tx"))
# Here: Check if subjectHits(ov) are in tx_cds names.
# If yes, CDS is TRUE
ovdf$cds <- ovdf$subjectHits %in% names(tx_cds)

tx_exons <- exonsBy(TxDb.An, by = c("tx"))
ovdf$exons <- ovdf$subjectHits %in% names(tx_exons)

tx_introns <- intronsByTranscript(TxDb.An)
ovdf$introns <- ovdf$subjectHits %in% names(tx_introns)
#This seems wrong. They ALL register as introns? But also as exons?
# Ah, no. This is of course because the transcript CONTAINS both exons and introns.
# My approach is flawed. have to rethink.
#
# I CANNOT use the transcript number to fetch these, that is clear.
# The overlap function HAS to come into play. I can find overlaps for transcripts,
# cds, exons and promoters. (Also for miRNA and tRNA, but not necessary atm).
# However, I cannot get intron or UTR information from this.
#
# Can we use the 'nearest' methods in genomic ranges?
# Test
gr[1]
nearest(gr[1], genes(TxDb.An))
# No, nearest makes no sense. I want Overlap! (And this requires a Genomic Ranges,
# not a Genomic List object)
#
#
# I can generate the cds/exon tags, and from comparison with genes(), I can add an 
# intergenic tag as well.
# However, to get UTR I need to go by transcript.
tx_fiveUTR <- fiveUTRsByTranscript(TxDb.An)
# Is there a method to compare each of my ranges with the ranges in the transcript to find
# overlap? 
# Findoverlaps should be able to do this.
# 
# I can attach the IRanges to the df object and extract them to locate overlaps
# in each of the transcripts. If I receive 1 or more hits, it's TRUE, otherwise FALSE

findOverlaps(gr[1], tx_fiveUTR$`2`) %>% queryHits %>% length

# Rebuilding the df:
tx <- transcripts(TxDb.An)
ov_tx <- findOverlaps(gr, tx)
ovdf <- tibble(queryHits = queryHits(ov_tx), subjectHits = subjectHits(ov_tx))
ovdf$tx_name <- elementMetadata(tx)[subjectHits(ov_tx), "tx_name"]
### Adding this directly to the original df instead of the new one.
# new df might not serve a purpose if I can circumvent it.
# have to see what happens to the 'bytranscript' though.
ov_cds <- findOverlaps(gr, cds(TxDb.An))
df$cds <- FALSE
df$cds[queryHits(ov_cds)] <- TRUE

ov_exons <- findOverlaps(gr, exons(TxDb.An))
df$exons <- FALSE
df$exons[queryHits(ov_exons)] <- TRUE

ov_promoters <- findOverlaps(gr, promoters(TxDb.An))
df$promoters <- FALSE
df$promoters[queryHits(ov_promoters)] <- TRUE


tx_fiveUTR <- fiveUTRsByTranscript(TxDb.An)

# Now we get to 5UTR again.
#     make version of this:
findOverlaps(gr[1], tx_fiveUTR$`2`) %>% queryHits %>% length %>% as.logical()
#Index of transcript hits:
ov_tx
#For every pair of values here, use findOverlaps. Can use purrr's map2.

map2(.x = queryHits(ov_tx_sub),
     .y = subjectHits(ov_tx_sub),
     .f = function(x, y) findOverlaps(gr[x], 
                                      tx_fiveUTR[as.character(y)]) %>% return())
#We have a problem here, but it makes no sense to me:
# I assumed that the name of the fiveUTRByTranscript indicated the transcript ID
# but several transcript IDs that I have in ov_tx (overlap of gr and the transcript database)
# are NOT in the tx_fiveUTR object.
# The reason must be that there is no fiveUTR in this transcript?
# In that case, I need to account for this possibility. How many of the transcripts 
# are MISSING in the fiveUTRByTranscript list?
subjectHits(ov_tx) %in% names(tx_fiveUTR) %>% sum() #4165 out of 5891.
# So MOST have the fiveUTR.
#
# This seems plausible. Modify the map2 to account for this eventuality.
tx_fiveUTR[as.character(subjectHits(ov_tx)[7])]

findOverlaps(gr[queryHits(ov_tx)[8]], tx_fiveUTR[as.character(subjectHits(ov_tx)[8])])

ov_fiveUTR <-
  map2(.x = queryHits(ov_tx),
       .y = subjectHits(ov_tx),
       .f = function(x, y) {
         if (as.character(y) %in% names(tx_fiveUTR)) {
           findOverlaps(gr[x], tx_fiveUTR[as.character(y)]) %>% 
             queryHits %>% 
             length %>% 
             as.integer() %>% #if not, it converts to NULL
             #as.logical %>%
             return()
         } else {
             return(0)}
         })
# This works: Now to add it to the dataframe. Due to the length, it must first be added
# to the ovdf that was created to handle 'per transcript' data.

ov_fiveUTR %>% unlist() %>% sum() #only 80 hits in the end.

ovdf$fiveUTR <- ov_fiveUTR %>% unlist()

# This must be repeated for threeUTRByTranscript and intronsByTranscript.

# Finally, investigate if the result can be summarized with dplyr.
#These two functions do almost the same thing. Summarize is more appropriate, but
# mutate preserves the tx_name which may be interesting. Could perhaps keep it iwth summarize too.
ovdf %>% group_by(queryHits) %>% summarize(fiveUTR = paste0(fiveUTR, collapse=""))
ovdf %>% group_by(queryHits) %>% mutate(fiveUTR_paste = paste0(fiveUTR, collapse="")) %>% dplyr::slice(1)

# Next, have to make a more fitting method of fitting together the fiveUTR tags to count how many
# are present: This could simply be done using numbers of hits and summing. Then,
# it can be added to the dataframe by index.

#########################
# Decided to use numeric values and to summarize by the mean value.
#By also reporting the transcript names, it's easy to count how many there are.

ovdf_sum <- 
  ovdf %>% 
  group_by(queryHits) %>% 
  summarize(tx_name = paste0(tx_name, collapse=";"), 
            fiveUTR = mean(fiveUTR))

df <- tibble::rowid_to_column(df, "rowID")
df <- left_join(df, ovdf_sum, by = c("rowID" = "queryHits"))

### Now repeat for 3UTR

tx_threeUTR <- threeUTRsByTranscript(TxDb.An)
ov_threeUTR <-
  map2(.x = queryHits(ov_tx),
       .y = subjectHits(ov_tx),
       .f = function(x, y) {
         if (as.character(y) %in% names(tx_threeUTR)) {
           findOverlaps(gr[x], tx_threeUTR[as.character(y)]) %>% 
             queryHits %>% 
             length %>% 
             as.integer() %>% #if not, it converts to NULL
             #as.logical %>%
             return()
         } else {
           return(0)}
       })
ov_threeUTR %>% unlist() %>% sum() #1346 hits
ovdf$threeUTR <- ov_threeUTR %>% unlist()
ovdf_sum <- 
  ovdf %>% 
  group_by(queryHits) %>% 
  summarize(threeUTR = mean(threeUTR))
df <- left_join(df, ovdf_sum, by = c("rowID" = "queryHits"))

### Introns

tx_introns <-intronsByTranscript(TxDb.An)
ov_introns <-
  map2_int(.x = queryHits(ov_tx),
       .y = subjectHits(ov_tx),
       .f = function(x, y) {
         if (as.character(y) %in% names(tx_introns)) {
           findOverlaps(gr[x], tx_introns[as.character(y)]) %>% 
             queryHits %>% 
             length %>% 
             as.integer() %>% #if not, it converts to NULL
             #as.logical %>%
             return()
         } else {
           return(0)}
       })
ov_introns %>% unlist() %>% sum() #2863 hits
ovdf$introns <- ov_introns %>% unlist()
ovdf_sum <- 
  ovdf %>% 
  group_by(queryHits) %>% 
  summarize(introns = mean(introns))
df <- left_join(df, ovdf_sum, by = c("rowID" = "queryHits"))

#Finally: Find if any clusters are not labeled by any of these. If so, label as 'intergenic'
# In general, it should be any cluster that has no TRANSCRIPT associated.
# the rest are fine.

#df_an <- 
#  left_join(df, ovdf_sum, by = c("rowID" = "queryHits"))

df_an %>% filter(is.na(Symbol))
df_an %>% filter(exons.y == 0 & introns.y == 0) #gives NO outputs.

# Thus, there's no intergenic regions? How else could I detect them?
# Don't think it's interesting enough.

#### And this finishes the df.
#### Now to write this into a more concise version that can be called from shell.





########################
## 25/4
## Got erronious annotation output so have to find where it went wrong
## 
## Annotation output found AR ranges on chromosome 7 instead of chromosome X.
## The ranges were in the right area, but were extended beyond the gene as defined on UCSC.
#
#
library(tibble)
library(tidyr)
library(dplyr)
library(readr)
library(purrr)
library(GenomicRanges)
library(GenomicFeatures)
library(AnnotationHub)

setwd("/home/mlorentz/Documents/PAR-CLIP/snakemake_pipeline/PARalyzer/output")

#This file contains only a few ranges, most on AR.
#input <- "test2.csv"
input <- "SRR3502967.clusters.csv"
genomeVersion <- "hg19"

df <- read_csv(file = input)
df <- rowid_to_column(df, "rowIndex")



r_chr <- base::colnames(df) %>% grep(pattern = "Chromosome")
r_start <- base::colnames(df) %>% grep(pattern = "Start")
r_end <- base::colnames(df) %>% grep(pattern = "End")
r_ID <- base::colnames(df) %>% grep(pattern = "ID")
r_str <- base::colnames(df) %>% grep(pattern = "Strand")

gr <- 
  GRanges(seqnames = dplyr::select(df, r_chr) %>% pull %>% Rle,
          ranges = IRanges(start = dplyr::select(df, r_start) %>% pull,
                           end = dplyr::select(df, r_end) %>% pull,
                           names = dplyr::select(df, r_ID) %>% pull),
          strand = dplyr::select(df, r_str) %>% pull %>% Rle)

hub <- AnnotationHub()
OrgDb_An <- 
  with(hub, 
       subset(hub, species == "Homo sapiens" 
              & rdataclass == "OrgDb")) %>% .[[1]]
TxDb_An <-
  with(hub, 
       subset(hub, 
              species == "Homo sapiens" 
              & rdataclass == "TxDb" 
              & paste0(genomeVersion, ".knownGene") %>% 
                grepl(., title))
  ) %>% 
  .[[1]]

#Check that AR is location is defined correctly.

testg <- genes(TxDb_An)
testg[testg$gene_id == 367,]

#It is indeed. 
#So try to find overlaps with the AR gene.

findOverlaps(gr, testg[testg$gene_id == 367,])
#gives 147 overlaps. Again, we expect plenty of AR. But what happens when we annotate?

#Find annotation by transcript (OF AR)
tx <- transcripts(TxDb_An)
ov_tx <- findOverlaps(testg[testg$gene_id == 367,], tx)
#Make df to hold overlap hits, mapping between query index and database transcripts:
# HOLDING AR
ovdf <- 
  tibble(queryHits = queryHits(ov_tx), subjectHits = subjectHits(ov_tx)) %>%
  tibble::add_column(tx_name = elementMetadata(tx)[subjectHits(ov_tx), "tx_name"])


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
## Compute annotation
ovdf_an <- add_column(ovdf, 
                      cds = annotationHitsByTranscript(type = "cds"),
                      exons = annotationHitsByTranscript(type = "exons"),
                      introns = annotationHitsByTranscript(type = "introns"),
                      fiveUTRs = annotationHitsByTranscript(type = "fiveUTRs"),
                      threeUTRs = annotationHitsByTranscript(type = "threeUTRs"))
##### NOTE: Of 147 hits, only 7 hit a transcript. Can this be right????
#     I'm missing promoters. Is this really it?
#     The transcript names look right, though: They correspond to those on USCS when searched for.
# I think what I'm seeing that it here collapses to 7 transcripts because these are the 7 transcripts
# found in AR.

ovdf$tx_name %>% unique
#Yes, that's exactly why.
# However, none of the hits here return anything for genomic region type. This merits investigation.

###################
# Starting new attack: (with a trimmed down df consisting mainly of AR ranges)

# trim down to make it easily managable.
# Worked trimmed down, taking the whole of the test.csv
#df <- df[1:10,]

gr <- 
  GRanges(seqnames = dplyr::select(df, r_chr) %>% pull %>% Rle,
          ranges = IRanges(start = dplyr::select(df, r_start) %>% pull,
                           end = dplyr::select(df, r_end) %>% pull,
                           names = dplyr::select(df, r_ID) %>% pull),
          strand = dplyr::select(df, r_str) %>% pull %>% Rle)

#Get gene overlaps
tx_genes <- genes(TxDb_An)
ov <- findOverlaps(gr, tx_genes)
#Get the gene IDs and symbols:
gn <- elementMetadata(tx_genes)[subjectHits(ov),"gene_id"]
sb <- AnnotationDbi::select(OrgDb_An, unlist(gn), "SYMBOL", "ENTREZID") %>% 
  as_tibble() %>%
  add_column("rowIndex" = queryHits(ov))
  #rowid_to_column("rowIndex") #here is the problem. Can't give it a rowindex,must be given queryHits

#Sanity check to ensure that the queryHits are unique - ie, no double entry
#from multiple genes mapped to the same range
if (any(duplicated(queryHits(ov)))) {
  print("WARNING - non-unique hits on gene")
  sb <- 
    sb %>%
    group_by(rowIndex) %>%
    summarize(ENTREZID = paste0(ENTREZID, collapse = ";"),
              SYMBOL = paste0(SYMBOL, collapse = ";"))
}

df_upd <- 
  left_join(df, sb, by = "rowIndex") %>%
  dplyr::rename(Symbol = SYMBOL, EntrezID = ENTREZID) %>%
  replace_na(list(EntrezID = "-", Symbol = "Intergenic"))

######### LOCATED ERROR HERE:
df_upd %>% print(n=37)
# Two of the AR symbols have been shifted into the ch7 part. But this is because
# that there are two rows of NA at the very bottom - that should have been in chr7!
# Looking at the findOverlaps output, there is NO overlap for the two first queryHits
# But instead of returning NA here, their space remained empty.

#So the solution must be to make sure that the NA values are printed properly 
# (for intergenic regions). (Might still want to investigate why the order of the rows
# are changed for the annotated output, however.)



#Find annotation by transcript:
tx <- transcripts(TxDb_An)
ov_tx <- findOverlaps(gr, tx)
#Make df to hold overlap hits, mapping between query index and database transcripts:
ovdf <- 
  tibble(queryHits = queryHits(ov_tx), subjectHits = subjectHits(ov_tx)) %>%
  tibble::add_column(tx_name = elementMetadata(tx)[subjectHits(ov_tx), "tx_name"])


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











## Put annotation back into the original dataframe.
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







# Check that annotation of genes work:
tx_genes <- genes(TxDb_An)
ov <- findOverlaps(gr, tx_genes)
#Selecting first pair.
gr[48,] %>% ranges()
tx_genes[17162,] %>% ranges()
#This checks out, they definitely overlap.








