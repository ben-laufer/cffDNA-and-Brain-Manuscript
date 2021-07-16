# Consensus DMR and human NDD/ASD overlaps
# Ben Laufer

.libPaths("/share/lasallelab/programs/DMRichR/R_4.0")
AnnotationHub::setAnnotationHubOption("CACHE", "/share/lasallelab/programs/DMRichR/R_4.0")

# Use custom built masked BSgenome for rheMac10
packages <- c("DMRichR", "regioneR", "BSgenome.Mmulatta.UCSC.rheMac10.masked","magrittr","Biostrings",
              "BSgenome.Hsapiens.UCSC.hg38.masked", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db")
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))

setwd("/share/lasallelab/Ben/Obesity/meta")
dir.create("overlaps")

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

# Make consensus regions --------------------------------------------------

getConsensus <- function(contrasts){
  contrasts %>% 
    purrr::map_dfr(function(dmrs){
      get(dmrs) %>%
        dplyr::as_tibble()
    }) %>%
    plyranges::as_granges() %>% 
    GenomicRanges::sort() %>%
    plyranges::reduce_ranges()
}

purrr::walk(c("cffDNA", "brains"), function(tissue){
  get(tissue) %>%
    getConsensus() %>% 
    assign(glue::glue("{tissue}_consensus"), ., envir = .GlobalEnv)
})

# cffDNA and brain overlaps -----------------------------------------------

# Remove contigs (Solution for the Rhesus pieces)
# https://support.bioconductor.org/p/83588/

keepBSgenomeSequences <- function(genome, seqnames)
{
  stopifnot(all(seqnames %in% seqnames(genome)))
  genome@user_seqnames <- setNames(seqnames, seqnames)
  genome@seqinfo <- genome@seqinfo[seqnames]
  genome
}

BSgenome.Mmulatta.UCSC.rheMac10.masked
genome <- BSgenome.Mmulatta.UCSC.rheMac10.masked
sequences_to_keep <- paste0("chr", c(1:20, "X", "Y"))
genome <- keepBSgenomeSequences(BSgenome.Mmulatta.UCSC.rheMac10.masked, sequences_to_keep)
genome

# Count overlaps
regioneR::numOverlaps(A = cffDNA_consensus, 
                      B = brains_consensus,
                      count.once = TRUE)

# Test overlaps
pt <- regioneR::overlapPermTest(A = cffDNA_consensus, 
                                B = brains_consensus,
                                alternative = "greater",
                                genome = genome,
                                ntimes = 10000,
                                count.once = TRUE)

pdf("overlaps/cffDNA_brain_Overlap.pdf", height = 7.50, width = 11.50)
plot(pt)
dev.off()

# Local Z-score analysis
lz <- regioneR::localZScore(A = cffDNA_consensus, 
                            B = brains_consensus,
                            pt = pt, 
                            window = 10*mean(width(brains_consensus)), 
                            step = mean(width(brains_consensus))/2,
                            count.once = TRUE)

pdf("overlaps/cffDNA_Brain_Overlap_Zscore.pdf", height = 7.50, width = 11.50)
plot(lz)
dev.off()

save(pt, lz, file = "overlaps/regioneR_cffDNA_Brain_overlap.RData")

# Annotate overlaps
c(cffDNA_consensus, brains_consensus) %>%
  plyranges::reduce_ranges() %>%
  plyranges::filter_by_overlaps(cffDNA_consensus) %>%
  plyranges::filter_by_overlaps(brains_consensus) %>%
  DMRichR::annotateRegions(TxDb = AnnotationHub::AnnotationHub()[["AH83244"]],
                           annoDb = "org.Mmu.eg.db") %>%
  openxlsx::write.xlsx(file = "overlaps/brain_cffDNA_coordinate_overlaps_annotated.xlsx")

# liftOver to hg38 --------------------------------------------------------
# R version of liftOver isn't the same algorithm and gives worse results for regions

system("rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/liftOver/rheMac10ToHg38.over.chain.gz .")
system("gunzip rheMac10ToHg38.over.chain.gz")

regions <- c("brains_regions_consensus","brains_sigRegions_consensus",
             "cffDNA_regions_consensus","cffDNA_sigRegions_consensus")

genomes <- c("hg38")

# Make bed files
purrr::walk(regions,
            function(regions){
              
              print(glue::glue("Processing {regions}"))
              
              get(regions) %>% 
                DMRichR::gr2bed(glue::glue("{regions}_rheMac10.bed"))
            })

# liftOver
tidyr::crossing(regions,
                genomes) %>% 
  purrr::pwalk(function(regions, genomes){
    
    print(glue::glue("Lifting over {regions} to {genomes}"))
    
    system(glue::glue("liftOver \\
                       -minMatch=0.1 \\
                       {regions}_rheMac10.bed \\
                       rheMac10To{Hmisc::capitalize(genomes)}.over.chain \\
                       {regions}_{genomes}.bed \\
                       {regions}_unmapped_{genomes}.txt"))
    
    print(glue::glue("Finished lifting over {regions} to {genomes}"))
    
  })

# Load liftOver
loadBED <- function(name){
  sigRegions <- readr::read_tsv(paste0(name,"_hg38.bed"),
                                col_names = c("chr", "start", "end"))  %>% 
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  assign(name, sigRegions, envir = .GlobalEnv)
  rm(sigRegions)
}

contrasts <- c("cffDNA_sigRegions_consensus", "brains_sigRegions_consensus")
purrr::walk(contrasts, loadBED)

# Load human DMRs ---------------------------------------------------------

loadDMRs <- function(name){
  load(glue::glue("humans/{name}/RData/DMRs.RData"))
  assign(name, sigRegions, envir = .GlobalEnv)
  rm(sigRegions, regions)
}

contrasts <- c("ASD_placenta_male", "ASD_placenta_female", "ASD_brain_male",
               "Rett_brain_female", "Dup15_brain_male", "DS_brain_male")

purrr::walk(contrasts, loadDMRs)

# Test for human NDD/ASD overlaps -----------------------------------------

ASDoverlaps <- function(contrast = contrast,
                        dataset = dataset){
  
  print(glue::glue("Running permutation test for {contrast} within {dataset}"))
  
  pt <- regioneR::overlapPermTest(A = get(contrast), 
                                  B = get(dataset), 
                                  alternative = "greater",
                                  genome = "hg38",
                                  ntimes = 10000,
                                  count.once = TRUE)
  
  tibble::tibble(contrast = !!contrast,
                 dataset = !!dataset,
                 overlaps = pt$numOverlaps$observed,
                 zscore = pt$numOverlaps$zscore,
                 pvalue = pt$numOverlaps$pval)
}

## Run --------------------------------------------------------------------

results <- tidyr::crossing(contrast = c("cffDNA_sigRegions_consensus", "brains_sigRegions_consensus"),
                           dataset = c("ASD_placenta_male", "ASD_placenta_female", "ASD_brain_male",
                                       "Rett_brain_female", "Dup15_brain_male", "DS_brain_male")) %>% 
  purrr::pmap_dfr(ASDoverlaps) %>% 
  dplyr::mutate(p.adjust = p.adjust(pvalue)) 

save(results, file = "overlaps/NDDoverlaps.RData")

## Annotate overlaps ------------------------------------------------------

annotateOverlaps <- function(contrast = contrast,
                             dataset = dataset){
  
  print(glue::glue("Annotating {contrast} within {dataset}"))
  
  get(contrast) %>%
    plyranges::filter_by_overlaps(get(dataset)) %>%
    DMRichR::annotateRegions(TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, 
                             annoDb = "org.Hs.eg.db") %>%
    dplyr::as_tibble() %>%
    dplyr::select(-geneId) %>%
    tibble::add_column(contrast = contrast, dataset = dataset, .before = 1) %>%
    dplyr::mutate(contrast = dplyr::case_when(contrast == "cffDNA_sigRegions_consensus" ~ "cffDNA",
                                              contrast == "brains_sigRegions_consensus" ~ " Brain"),
                  dataset = dplyr::case_when(dataset == "ASD_placenta_male" ~ "ASD Placenta (Male)",
                                             dataset == "ASD_placenta_female" ~ "ASD Placenta (Female)"
                                             dataset == "ASD_brain_male" ~ "ASD Brain (Male)",
                                             dataset == "Rett_brain_female" ~ "Rett Syndrome Brain (Female)",
                                             dataset == "Dup15_brain_male" ~ "Dup15q Syndrome Brain (Male)",
                                             dataset == "DS_brain_male" ~ "Down Syndrome Brain (Male)"))
}

### Run -------------------------------------------------------------------

tidyr::crossing(contrast = c("cffDNA_consensus", "brains_consensus"),
                dataset = c("ASD_placenta_male", "ASD_placenta_female", "ASD_brain_male",
                            "Rett_brain_female", "Dup15_brain_male", "DS_brain_male")) %>% 
  purrr::pmap_dfr(annotateOverlaps) %>%
  openxlsx::write.xlsx("overlaps/human_overlaps_annotated.xlsx")
