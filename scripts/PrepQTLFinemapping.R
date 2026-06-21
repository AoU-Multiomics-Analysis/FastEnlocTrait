library(readxl)
library(tidyverse)
library(data.table)
library(optparse)

FormatFastEnloc <- function(input_df) {
formatted <- input_df %>% 
      mutate(
    # clean variant string
    variant_clean = str_replace(variant, "^chrchr", "chr"),
    
    # parse chr, pos, ref, alt
    chrom = str_extract(variant_clean, "^chr[^_]+"),
    pos   = as.integer(str_extract(variant_clean, "(?<=_)\\d+(?=_)")),
    ref   = str_match(variant_clean, "^chr[^_]+_\\d+_([^_]+)_([^_]+)$")[,2],
    alt   = str_match(variant_clean, "^chr[^_]+_\\d+_([^_]+)_([^_]+)$")[,3],
    
    # variant ID
    id = paste(chrom, pos, ref, alt, sep = "_"),
    
    # convert L2 -> 2
    #cs_number = str_remove(cs_index, "^L"),
    cs_number = cs_id,
    # locus_id = gene:credible_set
    locus_id = paste0(molecular_trait_id, ":", cs_number),
    
    # no tissue, so use @=
    qtl_annot = paste0(
      locus_id,
      "@=",
      signif(pip, 6),
      "[",
      signif(cpip, 6),
      ":",
      number_variants,
      "]"
    )
  ) %>%
  group_by(chrom, pos, id, ref, alt) %>%
  summarise(
    qtl_annotation = paste(qtl_annot, collapse = "|"),
    .groups = "drop"
  ) %>%
  arrange(chrom, pos)
formatted  
}

####### PARSE ARGUMENTS #########
option_list <- list(
    optparse::make_option(c("--QTLData"), type = "character", default = NULL,
                          help = ""),
    optparse::make_option(c("--QTLType"), type = "character", default = NULL,
                          help = ""),
    optparse::make_option(c("--OutputFile"), type = "character", default = NULL,
                          help = "")
    )

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
QTLData <- opt$QTLData
QTLType <- opt$QTLType
OutputFile <- opt$OutputFile

if (is.null(QTLData)) stop("Must provide --QTLData")
if (is.null(QTLType)) stop("Must provide --QTLType")
if (is.null(OutputFile)) stop("Must provide --OutputFile")

######## LOAD DATA #############

if (QTLType == 'Expression') {
QTL <- fread(QTLData) %>% 
    select(molecular_trait_id,variant,pip,cs_id) %>% 
    mutate(molecular_trait_id = str_remove(molecular_trait_id,'\\..*')) %>%
    group_by(molecular_trait_id,cs_id) %>% 
    mutate(cpip = sum(pip),number_variants = dplyr::n())  %>% 
    ungroup() %>%
    FormatFastEnloc
} else if (QTLType == 'Splicing') {
QTL <- fread(QTLData) %>% 
    mutate(molecular_trait_id = stringr::str_extract(molecular_trait_id, "ENSG[0-9]+")) %>% 
    select(molecular_trait_id,variant,pip,cs_id) %>% 
    group_by(molecular_trait_id,cs_id) %>% 
    mutate(cpip = sum(pip),number_variants = dplyr::n()) %>%
    ungroup() %>% 
    FormatFastEnloc
} else if (QTLType == 'Protein') {
QTL <- fread(QTLData) %>% 
    filter(group == 'COMB') %>% 
    extract(
        col = variant,
        into = c("chrom", "pos", "ref", "alt"),
        regex = "^(chr[^:]+):(\\d+)([A-Z]+)_([A-Z]+)$",
        convert = TRUE
    ) %>% 
    mutate(variant = paste(chrom,pos,ref,alt,sep = '_')) %>% 
    mutate(molecular_trait_id = stringr::str_extract(molecular_trait_id, "ENSG[0-9]+")) %>% 
    select(molecular_trait_id,variant,pip,cs_id) %>% 
    group_by(molecular_trait_id,cs_id) %>% 
    mutate(cpip = sum(pip),number_variants = dplyr::n()) %>%
    ungroup() %>% 
    FormatFastEnloc
} else {
    stop("--QTLType must be one of Expression, Splicing, or Protein")
} 

# output extension should be .vcf.gz
QTL %>% write_tsv(OutputFile,col_names = FALSE)
