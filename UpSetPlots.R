# Upset Plots for cffDNA and brain gene symbol overlaps
# Ben Laufer

.libPaths("/share/lasallelab/programs/DMRichR/R_4.0")
packages <- c("AnnotationHub", "ensembldb", "DMRichR", "org.Mmu.eg.db", "UpSetR", "ComplexUpset", "tidyverse")
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))
AnnotationHub::setAnnotationHubOption("CACHE", "/share/lasallelab/programs/DMRichR/R_4.0")

setwd("/share/lasallelab/Ben/Obesity/meta")
dir.create("UpSets")

# Load DMRs ---------------------------------------------------------------

loadDMRs <- function(source, contrast){
  load(glue::glue("../{source}/DMRs/{contrast}/RData/DMRs.RData"))
  assign(contrast, sigRegions, envir = .GlobalEnv)
}

cffDNA <- c("GD45_OvC", "GD90_OvC", "GD120_OvC", "GD150_OvC",
            "GD45_RvO", "GD90_RvO", "GD120_RvO", "GD150_RvO",
            "GD45_PvO", "GD90_PvO", "GD120_PvO", "GD150_PvO")

purrr::walk(cffDNA,
            loadDMRs,
            source = "cffDNA")

brains <- c("Hippocampus_OvC", "Hypothalamus_OvC", "PrefrontalCortex_OvC",
            "Hippocampus_RvO", "Hypothalamus_RvO", "PrefrontalCortex_RvO",
            "Hippocampus_PvO", "Hypothalamus_PvO", "PrefrontalCortex_PvO")

purrr::walk(brains,
            loadDMRs,
            source = "brains")

# Get gene IDs ------------------------------------------------------------

extractGeneID <- function(regions = sigRegions,
                          TxDb = TxDb,
                          annoDb = annoDb){
  
  regions %>% 
    DMRichR::annotateRegions(TxDb = TxDb,
                             annoDb = annoDb) %>% 
    dplyr::select(geneId) %>%
    dplyr::distinct() %>% 
    purrr::pluck("geneId") %>%
    na.omit() %>% 
    return()
}

TxDb <- AnnotationHub::AnnotationHub()[["AH83244"]]
annoDb <- "org.Mmu.eg.db"

genes <- c(cffDNA, brains) %>%
  purrr::set_names() %>% 
  purrr::map(function(contrast){
  extractGeneID(get(contrast),
                TxDb = TxDb,
                annoDb = annoDb)
  })

# Plot --------------------------------------------------------------------

# https://github.com/krassowski/complex-upset/issues/51
# https://stackoverflow.com/questions/67094573/upsetr-sets-bar-interacting-different-color-and-its-sets-intersections

geneUpsetPlot <- function(list = list){
  list %>%
    UpSetR::fromList() %>%
    ComplexUpset::upset(.,
                        names(.),
                        n_intersections = 40, # Default from UpSetR
                        width_ratio = 0.2,
                        height_ratio = 0.95,
                        base_annotations = list(
                          'Intersection size'= intersection_size(
                            counts = FALSE
                          ) +
                            theme(
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              axis.line = element_line(colour = 'black'),
                              axis.text = element_text(size = 20),
                              axis.title = element_text(size = 20)
                            ) +
                            scale_y_continuous(expand = c(0, 0, 0.1, 0))
                        ),
                        sort_sets = FALSE,
                        queries = list(
                          upset_query(
                            intersect = c("Trimester 1", "Trimester 2", "Early Trimester 3", "Late Trimester 3",
                                          "Hippocampus", "Hypothalamus", "Prefrontal Cortex"),
                            color = '#E41A1C',
                            fill = '#E41A1C',
                            only_components = c('intersections_matrix', 'Intersection size')
                          ),
                          upset_query(
                            intersect = c("Hippocampus", "Hypothalamus", "Prefrontal Cortex"),
                            color = '#377EB8',
                            fill = '#377EB8',
                            only_components = c('intersections_matrix', 'Intersection size')
                          ),
                          upset_query(set = 'Trimester 1', fill = '#FF7F00'),
                          upset_query(set = 'Trimester 2', fill = '#FF7F00'),
                          upset_query(set = 'Early Trimester 3', fill = '#FF7F00'),
                          upset_query(set = 'Late Trimester 3', fill = '#FF7F00'),
                          upset_query(set = 'Hippocampus', fill = '#4DAF4A'),
                          upset_query(set = 'Hypothalamus', fill = '#4DAF4A'),
                          upset_query(set = 'Prefrontal Cortex', fill = '#4DAF4A')
                        ),
                        matrix = intersection_matrix(
                          geom = geom_point(
                            shape = 'circle filled',
                            size = 3.5,
                            stroke = 0.45
                          )
                        ),
                        set_sizes = (
                          upset_set_size(geom = geom_bar(color = 'black')
                          )
                        ),
                        themes = upset_modify_themes( # names(upset_themes) 
                          list(
                            'intersections_matrix' = theme(
                              axis.text = element_text(size = 20),
                              axis.title = element_blank()
                            ),
                            'overall_sizes' = theme(
                              axis.title = element_text(size = 14)
                            )
                          )
                        )
    )
}

combinedGenePlots <- function(contrast){
  
  pdf(glue::glue("UpSets/Combined_UpSet_Genes_{contrast}.pdf"),
      height = 6, width = 10, onefile = FALSE)
  
  print({
    genes %>%
      magrittr::extract(c(cffDNA, brains) %>% 
                          stringr::str_subset(contrast)) %>% 
      magrittr::set_names(c("Trimester 1", "Trimester 2", "Early Trimester 3", "Late Trimester 3",
                            "Hippocampus", "Hypothalamus", "Prefrontal Cortex")) %>% 
      rev() %>% 
      geneUpsetPlot()
  })
  
  dev.off()
}

purrr::walk(c("OvC", "RvO", "PvO"), combinedGenePlots)
