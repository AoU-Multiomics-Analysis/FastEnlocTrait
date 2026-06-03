library(readxl)
library(tidyverse)
library(dat.table)



CleanMVPData <- function(MVP_dat) {
MVPCleaned <- MVP_dat %>% 
    group_by(Locus,Population,`Population Signal`) %>% 
    mutate(cpip = sum(`CS-Level Pip`),
          alleles = str_extract(`MVP ID`, "[^:]+:[^:]+$"),
          CHR = paste0('chr',CHR)) %>%  
    mutate(variant = paste0(CHR,'_',BP38,'_',alleles)) %>%
    separate(alleles,into = c('ref','alt')) %>% 
    mutate(variant = str_replace(variant,':','_')) %>%
    mutate(cs_id = paste(Trait,Population,Locus,`Population Signal`,sep = '_')) %>% 
    ungroup() %>% 
    select(CHR,BP38,variant,ref,alt,cs_id,`CS-Level Pip`,cpip)   
 MVPCleaned   
}
make_fastenloc_input_MVP <- function(df,out_file = NULL,build = "b38") {
  out <- df %>%
    group_by(cs_id) %>%
    mutate(
      n_snps = n(),
      variant_id = paste(
        CHR,
        BP38,
        ref,
        alt,
        sep = "_"
      ),
      locus_string = paste0(
        cs_id,
        "=",
        format(`CS-Level Pip`,
               scientific = TRUE,
               digits = 4),
        "[",
        format(cpip,
               scientific = TRUE,
               digits = 4),
        ":",
        n_snps,
        "]"
      )
    ) %>%
    ungroup() %>%
    transmute(
      CHR,
      BP38,
      variant_id,
      ref,
      alt,
      locus_string
    )
  if (!is.null(out_file)) {
    write_tsv(
      out,
      out_file,
      col_names = FALSE
    )
  }
  out
}


