# Create supplementary tables
# By Ben Laufer

# Initialize --------------------------------------------------------------

.libPaths("/share/lasallelab/programs/DMRichR/R_4.1")
packages <- c("AnnotationHub", "ensembldb", "DMRichR", "org.Mmu.eg.db", "plyranges")
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))
AnnotationHub::setAnnotationHubOption("CACHE", "/share/lasallelab/programs/DMRichR/R_4.1")

# DMRs and Blocks ---------------------------------------------------------

updateAnno <- function(Source,
                       Sampling,
                       Contrast,
                       Type,
                       TxDb = TxDb,
                       annoDb = annoDb){
  
  print(glue::glue("Loading {Source} {Sampling} {Contrast} {Type}"))
  
  load(glue::glue("{Source}/DMRs/{Sampling}_{Contrast}/RData/{Type}.RData"))
  
  if(Type == "DMRs"){
    peaks <- sigRegions
  }else if(Type == "Blocks"){
    peaks <- blocks
  }
  
  # peaks <- dplyr::case_when(Type == "DMRs" ~ sigRegions,
  #                           Type == "Blocks" ~ blocks)
  
  peaks %>% 
    plyranges::mutate(direction = dplyr::case_when(stat > 0 ~ "Hypermethylated",
                                                   stat < 0 ~ "Hypomethylated"),
                      difference = round(beta/pi *100)) %>% 
    DMRichR::annotateRegions(TxDb = TxDb,
                             annoDb = annoDb) %>% 
    dplyr::mutate(Source = dplyr::case_when(Source == "brains" ~ "Brain",
                                            TRUE ~ Source),
                  Sampling = Sampling,
                  Contrast = Contrast) %>% 
    dplyr::select(Source, Sampling, Contrast, dplyr::everything())
}

## Run --------------------------------------------------------------------

TxDb <- AnnotationHub::AnnotationHub()[["AH83244"]]
annoDb <- "org.Mmu.eg.db"

setwd("/share/lasallelab/Ben/Obesity/")
dir.create("meta/supplementary_tables", recursive = TRUE)

contrasts <- tidyr::crossing("Source" = "cffDNA",
                             "Sampling" = forcats::as_factor(c("GD45", "GD90", "GD120", "GD150")),
                             "Contrast" = forcats::as_factor(c("OvC", "RvO", "PvO")),
                             "Type" = forcats::as_factor(c("DMRs", "Blocks"))) %>%
  dplyr::slice(-22) %>% # No background blocks for GD150 RvO, need to manually remove to prevent an error
  dplyr::bind_rows(tidyr::crossing("Source" = "brains",
                                   "Sampling" = forcats::as_factor(c("Hippocampus", "Hypothalamus", "PrefrontalCortex")),
                                   "Contrast" = forcats::as_factor(c("OvC", "RvO", "PvO")),
                                   "Type" = forcats::as_factor(c("DMRs", "Blocks"))))

purrr::walk(c("DMRs", "Blocks"), function(Type){

  contrasts %>% 
    dplyr::filter(Type == !!Type) %>% 
    purrr::pmap_dfr(updateAnno,
                    TxDb = TxDb,
                    annoDb = annoDb) %>%
    openxlsx::write.xlsx(glue::glue("meta/supplementary_tables/{Type}.xlsx"), overwrite = TRUE)
})

# Enrichments -------------------------------------------------------------

setwd("/Users/blaufer/Box Sync/Maternal Obesity/analyses/withGD45full")

## GO ---------------------------------------------------------------------

c("cffDNA", "brains") %>% 
  purrr::map_dfr(function(source){
    readxl::read_xlsx(glue::glue("consensusGO/GOfuncR_rrvgo_Ensembl_consenus_{source}.xlsx")) %>%
      dplyr::mutate(Source = dplyr::case_when(source == "brains" ~ "Brain",
                                              TRUE ~ source)) %>% 
      dplyr::select(Source, Term, `Gene Ontology`, "-log10.q-value" = `-log10.p-value`) %>% 
      dplyr::left_join(readxl::read_xlsx(glue::glue("consensusGO/GOfuncR_Ensembl_consensus_{source}.xlsx")) %>%
                         dplyr::select(Term = node_name, FWER = FWER_overrep)) %>%
      dplyr::arrange(`Gene Ontology`, dplyr::desc(`-log10.q-value`))
    }) %>%
  openxlsx::write.xlsx("Supplementary Tables/Supplementary Table 2.xlsx", overwrite = TRUE)

## PANTHER ----------------------------------------------------------------

c("cffDNA", "brains") %>% 
  purrr::map_dfr(function(source){
    readxl::read_xlsx(glue::glue("consensusGO/enrichr_Ensembl_consensus_{source}.xlsx"), sheet = "Panther_2016") %>%
      dplyr::mutate("Gene Ontology" = "Panther_2016",
                    "-log10.p-value" = -log10(Adjusted.P.value),
                    Source = dplyr::case_when(source == "brains" ~ "Brain",
                                              TRUE ~ source)) %>%
      dplyr::select(Source, "Gene Ontology", Term, "-log10.p-value", Combined.Score) %>%
      dplyr::mutate("Database" = dplyr::recode_factor(`Gene Ontology`,
                                                      "Panther_2016" = "Panther Pathways")) %>%
      dplyr::select(-`Gene Ontology`) %>% 
      dplyr::mutate(Term = Hmisc::capitalize(Term)) %>%
      dplyr::mutate(Term = stringr::str_remove(Term, "Homo sapiens.*$")) %>% 
      dplyr::mutate(Term = stringr::str_trim(Term)) 
    }) %>%
      openxlsx::write.xlsx("Supplementary Tables/Supplementary Table 3.xlsx", overwrite = TRUE)

## Motifs -----------------------------------------------------------------

c("cffDNA", "brains") %>% 
  purrr::map_dfr(function(source){
    readxl::read_xlsx(glue::glue("memes/{source}_meme_methylcytosine_results.xlsx")) %>%
      dplyr::mutate(Source = dplyr::case_when(source == "brains" ~ "Brain",
                                              TRUE ~ source),
                    pvalue = as.character(pvalue),
                    adj.pvalue = as.character(adj.pvalue),
                    evalue = as.character(evalue)) %>%
      dplyr::select(Source, motif_id, consensus, adj.pvalue, evalue)
  }) %>%
  openxlsx::write.xlsx("Supplementary Tables/Supplementary Table 4.xlsx", overwrite = TRUE)
