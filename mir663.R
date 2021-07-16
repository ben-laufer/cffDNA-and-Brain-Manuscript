# Values and plots for the mir-663 block
# Ben Laufer

.libPaths("/share/lasallelab/programs/DMRichR/R_4.0")
packages <- c("AnnotationHub", "ensembldb", "Gviz", "DMRichR", "plyranges", "magrittr", "tidyverse", "WGCNA")
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))
AnnotationHub::setAnnotationHubOption("CACHE", "/share/lasallelab/programs/DMRichR/R_4.0")

setwd("/share/lasallelab/Ben/Obesity/meta")
dir.create("mir663")

# Values for statistics ---------------------------------------------------
# For objects with all samples from a source

purrr::walk(c("cffDNA", "brains"), function(source){
  load(glue::glue("../{source}/DMRs/all/{source}_bsseq.RData"))
  assign(glue::glue("{source}.bs.filtered.bsseq"), bs.filtered.bsseq, envir = .GlobalEnv)
})

block <- data.frame("seqnames" = "chr20",
                    "start" = (29790471 + 17000),
                    "end" = 29823571 - 1500) %>%
  plyranges::as_granges()

purrr::walk(c("cffDNA", "brains"), function(source){
  regionMeth <- get(glue::glue("{source}.bs.filtered.bsseq")) %>%
    bsseq::getMeth(BSseq = .,
                   regions = block,
                   type = "smooth",
                   what = "perRegion") %>% 
    t() %>%
    dplyr::as_tibble() %>%
    cbind(., pData(get(glue::glue("{source}.bs.filtered.bsseq")))) %T>%
    write.csv(glue::glue("mir663/chr20_block_{source}_individual_smoothed_methylation.csv"))
  
  saveRDS(regionMeth, file = glue::glue("mir663/chr20_block_{source}_individual_smoothed_methylation.rds"))
})

# Plot --------------------------------------------------------------------
# For objects broken up by timepoint and region

cffDNA <- c("GD45",
            "GD90",
            "GD120",
            "GD150")

brains <- c("Hippocampus",
            "Hypothalamus",
            "PrefrontalCortex")

loadBsseq <- function(source, contrast){
  load(glue::glue("../{source}/DMRs/all/{contrast}_bsseq.RData"))
  assign(glue::glue("{contrast}.bs.filtered.bsseq"), bs.filtered.bsseq, envir = .GlobalEnv)
}

purrr::walk(cffDNA,
            loadBsseq,
            source = "cffDNA")

purrr::walk(brains,
            loadBsseq,
            source = "brains")

glue::glue("Assigning colors for plotting...")

purrr::walk(c(cffDNA, brains), function(source){
  bs.filtered.bsseq <- get(glue::glue("{source}.bs.filtered.bsseq"))
  pData(bs.filtered.bsseq) <- pData(bs.filtered.bsseq) %>%
    dplyr::as_tibble() %>% 
    dplyr::mutate("Group" = dplyr::case_when(stringr::str_detect(Group, "Restriction") ~ "Caloric Restriction",
                                             TRUE ~ as.character(Group))) %>% 
    dplyr::mutate("Group" = dplyr::recode_factor(Group,
                                                 "Obese" = "Obese",
                                                 "Caloric Restriction" = "Intervention",
                                                 "Pravastatin" = "Intervention", 
                                                 "Control" = "Control")) %>% 
    dplyr::mutate("col" = dplyr::case_when(Group == "Control" ~ "mediumblue",
                                           Group == "Obese" ~ "firebrick3",
                                           Group == "Intervention" ~ "darkgoldenrod1", #darkgreen
                                           Group ==  "Intervetnion" ~ "darkgoldenrod1"))
  assign(glue::glue("{source}.bs.filtered.bsseq"), bs.filtered.bsseq, envir = .GlobalEnv)
})

glue::glue("Creating plotting annotations...")
genome <- "rheMac10"
# AnnotationHub::display(AnnotationHub::query(AnnotationHub::AnnotationHub(), c("EnsDb", "Macaca mulatta")))
TxDb <- AnnotationHub::AnnotationHub()[["AH83244"]]
seqlevelsStyle(TxDb) <- "UCSC"
annotations <- GenomicRanges::GRangesList(CpGs = DMRichR::getCpGs(genome),
                                          Exons = DMRichR::getExons(TxDb),
                                          compress = FALSE)

testCovariate <- "Group"

purrr::walk(c(cffDNA, brains), function(source){
  
  print(glue::glue("Now plotting {source} chr20 block"))
  
  #Main region
  sigRegions <- data.frame("seqnames" = "chr20",
                           "start" = (29790471 + 17000),
                           "end" = 29823571 - 1500) %>%
    plyranges::as_granges()

  pdf(glue::glue("mir663/chr20_block_Ensembl_{source}_mainregion.pdf"), height = 4, width = 8)
  
  DMRichR::plotDMRs2(get(glue::glue("{source}.bs.filtered.bsseq")),
                     regions = sigRegions,
                     testCovariate = testCovariate,
                     extend = (end(sigRegions) - start(sigRegions) + 1)*0.1,
                     addRegions = sigRegions,
                     annoTrack = annotations,
                     regionCol = "#FF00001A",
                     lwd = 2,
                     qval = FALSE,
                     stat = FALSE,
                     main = NULL,
                     horizLegend = FALSE,
                     addLines = FALSE,
                     addTicks = FALSE)
  
  dev.off()
})

## Gene track -------------------------------------------------------------
# https://bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/ensembldb.html#7_Plotting_genetranscript_features_using_ensembldb_and_Gviz_and_ggbio

## Retrieving a Gviz compatible GRanges object with all genes
gr <- ensembldb::getGeneRegionTrackForGviz(TxDb,
                                           chromosome = "20",
                                           start = 29807471,
                                           end = 29822071)

mcols(gr) <- gr %>% 
  mcols() %>% 
  dplyr::as_tibble() %>% 
  dplyr::mutate(symbol = dplyr::case_when(symbol == "mml-mir-663-9" ~ "mml-mir-663-9",
                                          symbol == "5_8S_rRNA" ~ "5_8S_rRNA",
                                          TRUE ~ gene))

## Define a genome axis track
gat <- GenomeAxisTrack()

ideo_track <- IdeogramTrack(genome = "rheMac10", chromosome = "chr20")

pdf(file = "mir663/chr20_track.pdf", height = 4, width = 8)
plotTracks(list(ideo_track,
                gat,
                GeneRegionTrack(gr)),
           transcriptAnnotation = "symbol")
dev.off()

# Correlations ------------------------------------------------------------

load("datTraits.RData") # From WGCNA script

MEs <- readr::read_rds("mir663/chr20_block_cffDNA_individual_smoothed_methylation.rds") %>%
  dplyr::as_tibble() %>% 
  dplyr::select("mir663" = V1, Infant, Timepoint) %>%
  dplyr::filter(Timepoint %in% c("GD45", "GD90", "GD120", "GD150")) %>%
  dplyr::mutate(Infant = as.character(Infant)) %>% 
  tidyr::pivot_wider(names_from = Timepoint, values_from = mir663) %>%
  dplyr::left_join(readr::read_rds("chr20_Block_main_regionbrains_individual_smoothed_methylation.rds") %>%
                     dplyr::as_tibble() %>% 
                     dplyr::select("mir663" = V1, "Infant" = InfantID, Region) %>%
                     dplyr::mutate(Infant = as.character(Infant)) %>% 
                     tidyr::pivot_wider(names_from = Region, values_from = mir663) %>%
                     dplyr::select(Infant, Hippocampus, Hypothalamus, `Prefrontal Cortex`)) %>%
  tibble::column_to_rownames("Infant")

MEs <- MEs[order(match(rownames(MEs),rownames(datTraits))),]

## Maternal Serum ---------------------------------------------------------

getModule <- function(GD,MEs){
  
  print(glue::glue("Loading {GD}"))
  
  tidyGD <- dplyr::case_when(GD == "GD040" ~ "GD45",
                             GD == "GD090" ~ "GD90",
                             TRUE ~ as.character(GD))
  
  MEs <- MEs %>%
    dplyr::select(!!tidyGD)
  
  datTraits %>% 
    tibble::rownames_to_column(var = "InfantID") %>% 
    dplyr::left_join(readxl::read_xlsx("IDs.xlsx") %>%
                       dplyr::mutate(InfantID = as.character(InfantID))) %>% 
    dplyr::left_join(readr::read_csv("Maternal_plasma_meta.csv") %>%
                       dplyr::select(Exp,hsCRP,"IFN-g" = IFN_g,"IL-1RA" = IL_ra,"IL-6" = IL_6,
                                     "IL-8" = IL_8,"IL-10" = IL_10,"IL-12/IL-23 p40" = "IL_12/23_p40",
                                     "IL-13" = IL_13,"IL-15" = IL_15,
                                     "MIP-1a (CCL3)" = MIP_1a, sCD40L, Plasma_color, GD_targeted, Animal_ID) %>%
                       dplyr::left_join(readr::read_csv("Maternal_plasma_conc_adj.csv") %>%
                                          dplyr::select("2-Hydroxyisovalerate" = "2_Hydroxyisovalerate", 
                                                        Acetate, Formate, Glucose, Glutamate,
                                                        "N,N-Dimethylglycine" = N.N_Dimethylglycine, Proline, Urea,
                                                        Exp), by = "Exp") %>%
                       dplyr::filter(Plasma_color == "normal") %>%
                       dplyr::filter(GD_targeted == !!GD) %>%
                       dplyr::rename(MotherID = Animal_ID)) %>%
    dplyr::select(-Exp, -MotherID, -Plasma_color, -GD_targeted,) %>%
    dplyr::mutate(InfantID = as.character(InfantID)) %>%
    tibble::column_to_rownames(var = "InfantID") %>%
    as.matrix() %>%
    WGCNA::cor(MEs,
               .,
               use = "p")
}

moduleTraitCor <- purrr::map(c("GD040", "GD090", "GD120", "GD150"),
                             getModule,
                             MEs = MEs) %>%
  do.call(rbind, .) %>%
  magrittr::set_rownames(c("Trimester 1",
                           "Trimester 2",
                           "Early Trimester 3",
                           "Late Trimester 3"))

moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, 25)

#### Heatmap --------------------------------------------------------------

pdf("mir663/mir663_maternal_plasma.pdf", height = 4, width = 18)

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(moduleTraitCor) %>%
                 stringr::str_wrap(25),
               yLabels = rownames(moduleTraitCor),
               ySymbols = rownames(moduleTraitCor),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               main = paste("mir663-trait Relationships with Maternal Blood"),
               verticalSeparator.x = c(4,6,17))

dev.off()

## Brain ------------------------------------------------------------------

getModule <- function(region,MEs){
  
  print(glue::glue("Loading {region}"))
  
  tidyMeta <- dplyr::case_when(region == "Hippocampus" ~ "hippocampus",
                               region == "Prefrontal Cortex" ~ "p.cortex",
                               region == "Hypothalamus" ~ "hypothalamus")
  
  datTraits %>% 
    tibble::rownames_to_column(var = "InfantID") %>% 
    dplyr::left_join(readxl::read_xlsx("IDs.xlsx") %>%
                       dplyr::mutate(InfantID = as.character(InfantID))) %>%
    dplyr::left_join(readxl::read_xlsx(glue::glue("{tidyMeta}Lipids.xlsx")) %>%
                       dplyr::mutate(InfantID = as.character(InfantID))) %>%
    dplyr::left_join(readr::read_csv("infant_brain_meta.csv") %>%
                       dplyr::select(Exp_ID,InfantID = Animal_ID) %>% 
                       dplyr::mutate(InfantID = as.character(InfantID)) %>% 
                       dplyr::right_join(readr::read_csv(glue::glue("Infant_{tidyMeta}_conc_adj.csv")) %>%
                                           dplyr::select(Exp_ID,
                                                         Adenosine, Citrate, GTP, Inosine, NAD,
                                                         "O-Phosphocholine" = O_Phosphocholine),
                                         by = "Exp_ID")) %>%
    dplyr::select(-Exp_ID, -MotherID) %>%
    dplyr::mutate(InfantID = as.character(InfantID)) %>%
    tibble::column_to_rownames(var = "InfantID") %>%
    as.matrix() %>%
    WGCNA::cor(MEs %>%
                 dplyr::select(!!region),
               .,
               use = "p")
}

moduleTraitCor <- purrr::map(c("Hippocampus", "Hypothalamus", "Prefrontal Cortex"),
                             getModule,
                             MEs = MEs) %>%
  do.call(rbind, .) 

moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, 25) 

### Heatmap ---------------------------------------------------------------

pdf("mir663/mir663_brain.pdf", height = 3.5, width = 18)

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(moduleTraitCor) %>%
                 stringr::str_wrap(25),
               yLabels = rownames(moduleTraitCor),
               ySymbols = rownames(moduleTraitCor),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               main = paste("mir663-trait Relationships with Infant Brain"),
               verticalSeparator.x = c(4,6,17))

dev.off()
