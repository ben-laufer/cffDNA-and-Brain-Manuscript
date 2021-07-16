# Consensus DMR functional enrichments 
# Gene Ontology, PANTHER Pathways, and Methylation Sensitive Transcription Factor Motifs
# Ben Laufer

# Create consensus --------------------------------------------------------

.libPaths("/share/lasallelab/programs/DMRichR/R_4.1")
AnnotationHub::setAnnotationHubOption("CACHE", "/share/lasallelab/programs/DMRichR/R_4.1")
packages <- c("org.Hs.eg.db", "org.Mmu.eg.db", "DMRichR", "magrittr", "ensembldb", "enrichR", "Brobdingnag",
              "BSgenome.Mmulatta.UCSC.rheMac10", "memes", "tidyverse", "cowplot")
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))

setwd("/share/lasallelab/Ben/Obesity/meta")

loadRegions <- function(source, contrast){
  load(glue::glue("../{source}/DMRs/{contrast}/RData/DMRs.RData"))
  assign(paste0(contrast,"_sigRegions"), sigRegions, envir = .GlobalEnv)
  assign(paste0(contrast,"_regions"), regions, envir = .GlobalEnv)
}

cffDNA <- c("GD45_OvC", "GD90_OvC", "GD120_OvC", "GD150_OvC",
            "GD45_RvO", "GD90_RvO", "GD120_RvO", "GD150_RvO",
            "GD45_PvO", "GD90_PvO", "GD120_PvO", "GD150_PvO")

purrr::walk(cffDNA,
            loadRegions,
            source = "cffDNA")

brains <- c("Hippocampus_OvC", "Hypothalamus_OvC", "PrefrontalCortex_OvC",
            "Hippocampus_RvO", "Hypothalamus_RvO", "PrefrontalCortex_RvO",
            "Hippocampus_PvO", "Hypothalamus_PvO", "PrefrontalCortex_PvO")

purrr::walk(brains,
            loadRegions,
            source = "brains")

getConsensus <- function(contrasts){
  contrasts %>% 
    purrr::map_dfr(function(regions){
      get(regions) %>%
        dplyr::as_tibble()
    }) %>%
    plyranges::as_granges() %>% 
    GenomicRanges::sort() %>%
    plyranges::reduce_ranges()
}

purrr::walk(c("cffDNA", "brains"), function(tissue){
  get(tissue) %>%
    paste0(.,"_sigRegions") %>% 
    getConsensus() %>% 
    assign(glue::glue("{tissue}_sigRegions_consensus"), ., envir = .GlobalEnv)
  
    get(tissue) %>%
      paste0(.,"_regions") %>% 
      getConsensus() %>% 
      assign(glue::glue("{tissue}_regions_consensus"), ., envir = .GlobalEnv)
})

# GOfuncR -----------------------------------------------------------------

setwd("/Users/blaufer/Desktop/obesity/consensusGO")

# This won't finish as a function, it needs to be split manually and takes a few days
# 
# dir.create("consensusGO")
# 
# purrr::walk(c("cffDNA", "brains"),
#             function(contrast){
#               print(glue::glue("Running GOfuncR for {contrast}"))
# 
#               get(glue::glue("{contrast}_sigRegions_consensus")) %>%
#                 DMRichR::GOfuncR(regions = get(glue::glue("{contrast}_regions_consensus")) ,
#                         n_randsets = 100,
#                         upstream = 5000,
#                         downstream = 1000,
#                         TxDb = AnnotationHub::AnnotationHub()[["AH83244"]],
#                         annoDb = "org.Hs.eg.db") %T>%
#                 openxlsx::write.xlsx(glue::glue("GOfuncR_Ensembl_consensus_{contrast}.xlsx")) %>%
#                 purrr::map( ~ dplyr::mutate_at(.x, dplyr::vars(dplyr::one_of("raw_p_overrep")), p.adjust)) %>%
#                 DMRichR::slimGO(tool = "GOfuncR",
#                                 annoDb = "org.Hs.eg.db",
#                                 plots = FALSE) %T>%
#                 openxlsx::write.xlsx(file = glue::glue("GOfuncR_rrvgo_Ensembl_consenus_{contrast}.xlsx")) %>%
#                 DMRichR::GOplot() %>%
#                 ggplot2::ggsave(glue::glue("GOfuncR_plot_Ensembl_conensus_{contrast}.pdf"),
#                                 plot = .,
#                                 device = NULL,
#                                 height = 8.5,
#                                 width = 10)
#             })

## Brain ------------------------------------------------------------------

brains_sigRegions_consensus %>% 
  DMRichR::GOfuncR(regions = brains_regions_consensus,
                   n_randsets = 100,
                   upstream = 5000,
                   downstream = 1000,
                   TxDb = AnnotationHub::AnnotationHub()[["AH83244"]],
                   annoDb = "org.Hs.eg.db") %T>%
  openxlsx::write.xlsx(glue::glue("GOfuncR_Ensembl_consensus_brains.xlsx")) %>% 
  purrr::map( ~ dplyr::mutate_at(.x, dplyr::vars(dplyr::one_of("raw_p_overrep")), p.adjust)) %>%
  DMRichR::slimGO(tool = "GOfuncR",
                  annoDb = "org.Hs.eg.db",
                  plots = FALSE) %T>%
  openxlsx::write.xlsx(file = glue::glue("GOfuncR_rrvgo_Ensembl_consenus_brains.xlsx")) %>% 
  DMRichR::GOplot() %>% 
  ggplot2::ggsave(glue::glue("GOfuncR_plot_Ensembl_conensus_brains.pdf"),
                  plot = .,
                  device = NULL,
                  height = 8.5,
                  width = 10)

## cffDNA ----------------------------------------------------------------

cffDNA_sigRegions_consensus %>% 
  DMRichR::GOfuncR(regions = cffDNA_regions_consensus,
                   n_randsets = 100,
                   upstream = 5000,
                   downstream = 1000,
                   TxDb = AnnotationHub::AnnotationHub()[["AH83244"]],
                   annoDb = "org.Hs.eg.db") %>%
  openxlsx::write.xlsx(glue::glue("GOfuncR_Ensembl_consensus_cffDNA.xlsx")) %>% 
  purrr::map( ~ dplyr::mutate_at(.x, dplyr::vars(dplyr::one_of("raw_p_overrep")), p.adjust)) %T>%
  DMRichR::slimGO(tool = "GOfuncR",
                  annoDb = "org.Hs.eg.db",
                  plots = FALSE) %T>%
  openxlsx::write.xlsx(file = glue::glue("GOfuncR_rrvgo_Ensembl_consenus_cffDNA.xlsx")) %>% 
  DMRichR::GOplot() %>% 
  ggplot2::ggsave(glue::glue("GOfuncR_plot_Ensembl_conensus_cffDNA.pdf"),
                  plot = .,
                  device = NULL,
                  height = 8.5,
                  width = 10)

# enrichR -----------------------------------------------------------------

enrichR:::.onAttach()

purrr::walk(c("cffDNA", "brains"), function(contrast){
  print(glue::glue("Running enrichR for {contrast}"))
  get(paste0(contrast,"_sigRegions_consensus")) %>%
    DMRichR::annotateRegions(TxDb = AnnotationHub::AnnotationHub()[["AH83244"]],
                             annoDb = "org.Mmu.eg.db") %>%  
    dplyr::select(geneSymbol) %>%
    purrr::flatten() %>%
    enrichR::enrichr(c("GO_Biological_Process_2018",
                       "GO_Cellular_Component_2018",
                       "GO_Molecular_Function_2018",
                       "KEGG_2019_Human",
                       "Panther_2016",
                       "Reactome_2016",
                       "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")) %>% 
    purrr::set_names(names(.) %>% stringr::str_trunc(31, ellipsis = "")) %T>%
    openxlsx::write.xlsx(file = glue::glue("enrichr_Ensembl_consensus_{contrast}.xlsx")) %>%
    DMRichR::slimGO(tool = "enrichR",
                    annoDb = "org.Mmu.eg.db",
                    plots = FALSE) %T>%
    openxlsx::write.xlsx(file = glue::glue("enrichr_rrvgo_Ensembl_consenus_{contrast}.xlsx")) %>% 
    DMRichR::GOplot() %>% 
    ggplot2::ggsave(glue::glue("enrichr_plot_Ensembl_conensus_{contrast}.pdf"),
                    plot = .,
                    device = NULL,
                    height = 8.5,
                    width = 10)
})

# MEME AME ----------------------------------------------------------------

# Run MEME AME (Analysis of Motif Enrichment) for methyl-senstive motifs
# https://meme-suite.org/meme/meme-software/5.3.3/meme-5.3.3.tar.gz
# METHYLCYTOSINE/yin2017.meme
# https://science.sciencemag.org/content/356/6337/eaaj2239

setwd("/share/lasallelab/Ben/Obesity/meta/memes")

#options(meme_bin = "/Users/blaufer/opt/miniconda3/pkgs/meme-5.3.0-py37pl5262he48d6b8_2/bin") 
options(meme_bin = "/software/meme/5.3.3/lssc0-linux/bin") 
memes::check_meme_install()
options(meme_db = "motif_databases/METHYLCYTOSINE/yin2017.meme")
goi <- BSgenome.Mmulatta.UCSC.rheMac10::BSgenome.Mmulatta.UCSC.rheMac10

methylMeme <- function(sigRegions = sigRegions,
                       regions = regions,
                       goi = goi){
  sigRegions %>%
    memes::get_sequence(goi) %>%
    memes::runAme(control = regions %>% 
                    memes::get_sequence(goi))
}

contrasts <- c("cffDNA", "brains")

parallel::mclapply(contrasts, 
                   function(contrast){
                     get(glue::glue("{contrast}_sigRegions_consensus")) %>% 
                       methylMeme(regions = get(glue::glue("{contrast}_regions_consensus")),
                                  goi = goi)  %>%
                       openxlsx::write.xlsx(file = glue::glue("{contrast}_meme_methylcytosine_results.xlsx"))
                   },
                   mc.cores = length(contrasts),
                   mc.silent = TRUE)

# Production plot ---------------------------------------------------------

GOplot <- function(contrast = contrast){
  readxl::read_xlsx(glue::glue("consensusGO/GOfuncR_rrvgo_Ensembl_consenus_{contrast}.xlsx")) %>% 
    dplyr::mutate("Database" = dplyr::recode_factor(`Gene Ontology`,
                                                    "Biological Process" = "GO Biological Process",
                                                    "Cellular Component" = "GO Cellular Component",
                                                    "Molecular Function" = "GO Molecular Function")) %>%
    dplyr::select(-`Gene Ontology`) %>% 
    dplyr::mutate(Term = stringr::str_trim(Term)) %>%
    dplyr::mutate(Term = Hmisc::capitalize(Term)) %>%
    dplyr::mutate(Term = stringr::str_trunc(Term, 50)) %>% 
    dplyr::group_by(Database) %>%
    dplyr::slice(1:7) %>%
    dplyr::ungroup() %>% 
    dplyr::mutate(Term = factor(.$Term, levels = unique(.$Term[order(forcats::fct_rev(.$Database), .$`-log10.p-value`)]))) %>%
    ggplot2::ggplot(aes(x = Term,
                        y = `-log10.p-value`,
                        fill = Database,
                        group = Database)) +
    ggplot2::geom_bar(stat = "identity",
                      position = position_dodge(),
                      color = "Black",
                      size = 0.25) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(expand = c(0, 0, 0.1, 0)) +
    ggsci::scale_fill_d3() +
    ggplot2::labs(y = expression("-log"[10](italic(q)))) +
    ggplot2::theme_classic() +
    ggplot2::theme(text = element_text(size = 12),
                   axis.title.y = element_blank(),
                   axis.title.x = element_text(size = 12),
                   plot.title = element_text(face = "bold", size = 12),
                   legend.position = "none") +
    ggtitle(dplyr::case_when(contrast == "brains" ~ "Brain DMRs",
                             contrast == "cffDNA" ~ "cffDNA DMRs")) %>%
    return()
}

pathwayPlot <- function(contrast = contrast){
  readxl::read_xlsx(glue::glue("consensusGO/enrichr_Ensembl_consensus_{contrast}.xlsx"), sheet = "Panther_2016") %>%
                       dplyr::mutate("Gene Ontology" = "Panther_2016",
                                     "-log10.p-value" = -log10(Adjusted.P.value)) %>%
                       dplyr::select("Gene Ontology", Term, "-log10.p-value", Combined.Score) %>%
    dplyr::mutate("Database" = dplyr::recode_factor(`Gene Ontology`,
                                                    "Panther_2016" = "Panther Pathways")) %>%
    dplyr::select(-`Gene Ontology`) %>% 
    dplyr::mutate(Term = Hmisc::capitalize(Term)) %>%
    dplyr::mutate(Term = stringr::str_remove(Term, "Homo sapiens.*$")) %>%
    dplyr::mutate(Term = stringr::str_trim(Term)) %>%
    dplyr::mutate(Term = stringr::str_trunc(Term, 50)) %>% 
    dplyr::slice(1:10) %>%
    dplyr::mutate(Term = factor(.$Term, levels = rev(.$Term))) %>%
    ggplot2::ggplot(aes(x = Term,
                        y = `-log10.p-value`,
                        fill = Database,
                        group = Database)) +
    ggplot2::geom_bar(stat = "identity",
                      position = position_dodge(),
                      color = "Black",
                      size = 0.25) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(expand = c(0, 0, 0.1, 0)) +
    ggplot2::scale_fill_manual(values = ggsci::pal_d3()(4)[4]) + 
    ggplot2::labs(y = expression("-log"[10](italic(q)))) +
    ggplot2::theme_classic() +
    ggplot2::theme(text = element_text(size = 12),
                   axis.title.y = element_blank(),
                   axis.title.x = element_text(size = 12),
                   plot.title = element_text(face = "bold", size = 12),
                   legend.position = "none") +
    ggtitle("PANTHER Pathways") %>%
    return()
}

amePlot <- function(contrast = contrast){
  readxl::read_xlsx(glue::glue("memes/{contrast}_meme_methylcytosine_results.xlsx")) %>%
    dplyr::slice(1:10) %>% 
    tidyr::separate(evalue, into = c("base", "exponent"), "e-") %>% 
    dplyr::rename("Motif ID" = motif_id) %>% 
    dplyr::mutate("Transcription Factor" = stringr::word(`Motif ID`, sep = "-"),
                  Domain = dplyr::case_when(stringr::str_detect(`Motif ID`, "eDBD") ~ "Extended DNA-binding Domain",
                                            stringr::str_detect(`Motif ID`, "FL") ~ "Full Length"),
                  Methylated = dplyr::case_when(stringr::str_detect(`Motif ID`, "methyl") ~ "Methylated",
                                                TRUE ~ "Not Methylated"),
                  "-log10.E-value" = -log10(as.numeric(.$base)*Brobdingnag::as.brob(10)^(-as.numeric(.$exponent))),
                  Contrast = !!contrast) %>% 
    ggplot2::ggplot(ggplot2::aes(x = reorder(`Motif ID`,`-log10.E-value`),
                                 y = `-log10.E-value`,
                                 fill = Methylated)) +
    ggplot2::geom_bar(stat = "identity",
                      position = ggplot2::position_dodge(),
                      color = "Black",
                      size = 0.25) +
    ggplot2::scale_fill_manual(values = c("Methylated" = "black",
                                          "Not Methylated" = "white")) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(expand = c(0, 0, 0.1, 0)) +
    ggplot2::labs(y = expression("-log"[10]("E-value"))) +
    ggplot2::theme_classic() +
    ggplot2::theme(text = element_text(size = 12),
                   axis.title.y = element_blank(),
                   axis.title.x = element_text(size = 12),
                   plot.title = element_text(face = "bold", size = 12),
                   legend.position = "none") +
    ggtitle("Transcription Factor Motifs") %>% 
    return()
}

## Run --------------------------------------------------------------------

setwd("..")

cowplot::plot_grid(GOplot("cffDNA"),
                   GOplot("brains"),
                   pathwayPlot("cffDNA"),
                   pathwayPlot("brains"),
                   amePlot("cffDNA"),
                   amePlot("brains"),
                   ncol = 2,
                   align = "v",
                   labels = c("A)", "B)", "C)", "D)", "E)", "F)"),
                   label_size = 12,
                   rel_heights = c(1.6,1,1),
                   scale = 0.95) %>%
  ggplot2::ggsave("consensus DMR functional enrichments.pdf",
                  plot = ., 
                  width = 12,
                  height = 10)
