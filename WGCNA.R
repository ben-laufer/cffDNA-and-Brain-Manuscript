# WGCNA for DMRichR results
# Modified from: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html
# Ben Laufer

rm(list=ls())
options(scipen=999)

cat("\n[DM.R] Loading packages \t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

.libPaths("/share/lasallelab/programs/DMRichR/R_4.0")

packages <- c("DMRichR", "sva", "org.Mmu.eg.db",
              "matrixStats","Hmisc","splines","foreach","doParallel","fastcluster",
              "AnnotationHub", "ensembldb", "org.Mmu.eg.db", "plyranges",
              "dynamicTreeCut","survival","GO.db","preprocessCore","impute", "WGCNA")


stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))
AnnotationHub::setAnnotationHubOption("CACHE", "/share/lasallelab/programs/DMRichR/R_4.0")
rm(packages)

# Prepare data ------------------------------------------------------------

setwd("/share/lasallelab/Ben/Obesity/meta")

## Load DMRs --------------------------------------------------------------

loadDMRs <- function(source, contrast){
  load(glue::glue("../{source}/DMRs/{contrast}/RData/DMRs.RData"))
  assign(contrast, sigRegions, envir = .GlobalEnv)
}

cffDNA <- c("GD45_OvC", "GD90_OvC", "GD120_OvC", "GD150_OvC",
            "GD45_RvC", "GD90_RvO", "GD120_RvO", "GD150_RvO",
            "GD45_PvO", "GD90_PvO", "GD120_PvO", "GD150_PvO",
            "GD45_RvC", "GD90_RvC", "GD120_RvC", "GD150_RvC")

purrr::walk(cffDNA,
            loadDMRs,
            source = "cffDNA")

brains <- c("Hippocampus_OvC", "Hypothalamus_OvC", "PrefrontalCortex_OvC",
            "Hippocampus_RvO", "Hypothalamus_RvO", "PrefrontalCortex_RvO",
            "Hippocampus_PvO", "Hypothalamus_PvO", "PrefrontalCortex_PvO",
            "Hippocampus_RvC", "Hypothalamus_RvC", "PrefrontalCortex_RvC")

purrr::walk(brains,
            loadDMRs,
            source = "brains")

## Make consensus regions -------------------------------------------------

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

## Smooth hippocampal methylomes ------------------------------------------

system("cp ../cytosine_reports/*.gz ./")

purrr::walk(c("Hippocampus"), function(Region,
                                       testCovariate = "Group",
                                       adjustCovariate = NULL,
                                       matchCovariate = NULL,
                                       coverage = 1,
                                       cores = 10,
                                       perGroup = 0.75,
                                       genome = "rheMac10"){
  cat("\n[DMRichR] Loading Bismark genome-wide cytosine reports \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  start_time <- Sys.time()
  
  meta <- openxlsx::read.xlsx("sample_info_master.xlsx", colNames = TRUE) %>%
    dplyr::mutate_if(is.character, as.factor) %>%
    dplyr::filter(Region == !!Region)
  
  bs.filtered <- DMRichR::processBismark(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                                         meta = meta,
                                         testCovariate = testCovariate,
                                         adjustCovariate = adjustCovariate,
                                         matchCovariate = matchCovariate,
                                         coverage = coverage,
                                         cores = cores,
                                         perGroup = perGroup
  ) 
  
  glue::glue("Saving Rdata...")
  save(bs.filtered, file = glue::glue("{Region}_bismark.RData"))
  
  glue::glue("Filtering timing...")
  end_time <- Sys.time()
  end_time - start_time
  
  cat("\n[DMRichR] Smoothing individual methylation values \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
  start_time <- Sys.time()
  
  bs.filtered.bsseq <- bsseq::BSmooth(bs.filtered,
                                      BPPARAM = MulticoreParam(workers = cores,
                                                               progressbar = TRUE)
  )
  
  bs.filtered.bsseq
  
  glue::glue("Saving Rdata...")
  save(bs.filtered.bsseq, file = glue::glue("{Region}_bsseq.RData"))
  
  glue::glue("Individual smoothing timing...")
  end_time <- Sys.time()
  end_time - start_time
})

system("rm *.txt.gz")

## Extract methylation values ---------------------------------------------

extractMethyl <- function(dataset, consensus){
  load(glue::glue("{dataset}_bsseq.RData"))
  
  smoothed <- bs.filtered.bsseq %>%
      bsseq::getMeth(BSseq = .,
                     regions = consensus,
                     type = "smooth",
                     what = "perRegion") %>% 
    as.data.frame(check.names = FALSE) %>% 
    dplyr::bind_cols(as.data.frame(consensus), .)

  assign(glue::glue("{dataset}_meth"), smoothed, envir = .GlobalEnv)
  save(smoothed, file = glue::glue("WGCNA/{dataset}/meth_consensus.RData"))
}

dir.create("/share/lasallelab/Ben/Obesity/meta/WGCNA/Hippocampus", recursive = TRUE)

consensus <- c(brains_consensus, cffDNA_consensus) %>%
  plyranges::reduce_ranges()

extractMethyl(dataset = "Hippocampus",
              consensus = consensus)

# Start WGCNA -------------------------------------------------------------

options(stringsAsFactors = FALSE)
enableWGCNAThreads(60)
setwd("/share/lasallelab/Ben/Obesity/meta/WGCNA/Hippocampus")
load("meth_consensus.RData")

## WGCNA 1: Data input and cleaning ----

smoothed <- Hippocampus_meth

names(smoothed) <- smoothed %>%
  names() %>%
  gsub("\\-.*","",.)
  
WGCNA_data0 <- as.data.frame(t(smoothed[, c(6:length(smoothed))]))
names(WGCNA_data0) <- paste(smoothed$seqnames,":", smoothed$start, "-", smoothed$end, sep = "")
rownames(WGCNA_data0) <- names(smoothed[, c(6:length(smoothed))])

# Check for missing values and outliers
WGCNA_gsg <- goodSamplesGenes(WGCNA_data0, verbose = 3)
WGCNA_gsg$allOK

# Remove any offending regions and/or samples
if (!WGCNA_gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!WGCNA_gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(WGCNA_data0)[!WGCNA_gsg$goodGenes], collapse = ", ")));
  if (sum(!WGCNA_gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(WGCNA_data0)[!WGCNA_gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  WGCNA_data0 = WGCNA_data0[WGCNA_gsg$goodSamples, WGCNA_gsg$goodGenes]
}

# Cluster samples to check for outliers
sampleTree = hclust(dist(WGCNA_data0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
#sizeGrWindow(12,9)
pdf(file = "sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

## Manually remove outliers based on visualization (Change for each analysis)
# Plot a line to show what the cut would exclude
abline(h = 34, col = "red");
# Determine cluster under the line (use same value as line to cut)
clust = cutreeStatic(sampleTree, cutHeight = 34, minSize = 10)
table(clust)

# clust 1 contains the samples we want to keep, so remove outliers here (clust problem! minsize?)
keepSamples = (clust==1)
WGCNA_data = WGCNA_data0[keepSamples, ]
nGenes = ncol(WGCNA_data)
nSamples = nrow(WGCNA_data)

# Load and clean traits

sampleInfo <- readxl::read_xlsx("../sample_info_brain_master.xlsx") %>%
  dplyr::filter(Region == "Hippocampus") %>%
  dplyr::select(InfantID, Group, Obese, C_section, Foster, Cohort) %>%
  dplyr::mutate("Obese Only" = dplyr::case_when(Group == "Obese" ~ 1,
                                                TRUE  ~ 0)) %>% 
  dplyr::mutate(Control = dplyr::case_when(Group == "Control" ~ 1,
                                           TRUE  ~ 0)) %>% 
  dplyr::mutate(Obese = dplyr::case_when(Obese == "Obese" ~ 1,
                                         Obese == "Control" ~ 0)) %>% 
  dplyr::mutate(C_section = dplyr::case_when(C_section == "Yes" ~ 1,
                                             C_section == "No" ~ 0)) %>% 
  dplyr::mutate(Foster = dplyr::case_when(Foster== "Yes" ~ 1,
                                             Foster == "No" ~ 0)) %>% 
  dplyr::mutate(Intervention = dplyr::case_when(Group == "Restriction" ~ 1,
                                                Group == "Pravastatin" ~ 1,
                                                TRUE ~ 0)) %>% 
  dplyr::mutate(Restriction = dplyr::case_when(Group == "Restriction" ~ 1,
                                               TRUE ~ 0)) %>% 
  dplyr::mutate(Pravastatin = dplyr::case_when(Group == "Pravastatin" ~ 1,
                                               TRUE ~ 0)) %>% 
  dplyr::mutate(Year_1 = dplyr::case_when(Cohort == "YEAR 1" ~ 1,
                                               TRUE ~ 0)) %>% 
  dplyr::mutate(Year_2 = dplyr::case_when(Cohort == "YEAR 2" ~ 1,
                                          TRUE ~ 0)) %>% 
  #dplyr::select(-Group, -Cohort)  %>%
  dplyr::select(InfantID, Obese = "Obese Only", Control, Restriction, Pravastatin)

traits <- readxl::read_xlsx("../OBESITY Nov Pref.xlsx") %>% 
  dplyr::select(-COHORT,-"Group Assignment", -"contol obese/normal") %>%
  dplyr::rename(InfantID = INFANT_ID) %>%
  dplyr::mutate(InfantID = as.character(InfantID)) %>%
  dplyr::select(InfantID,
                "Recognition Memory: Abstract Stimuli" = "no. looks N/N+F",
                "Recognition Memory: Social Stimuli" = NOVPA)

meta <- sampleInfo %>%
  dplyr::inner_join(traits, by = "InfantID") %>% 
  as.data.frame() 

str(meta)

WGCNA_samples <- rownames(WGCNA_data)
traitRows <- match(WGCNA_samples, meta$InfantID)
datTraits <- meta[traitRows, -1]

rownames(datTraits) <- meta[traitRows, 1]
collectGarbage()

# Re-cluster samples
sampleTree2 = hclust(dist(WGCNA_data), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = TRUE);

# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")

dev.off()

save(WGCNA_data, datTraits, file = "WGCNA_1.RData") 
#load("WGCNA_1.RData")

## WGCNA 2: Automatic, one-step network construction and module detection ----

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 30, by = 2))
# Call the network topology analysis function
sft = pickSoftThreshold(WGCNA_data,
                        powerVector = powers,
                        networkType = "signed",
                        corFnc = "bicor", # "pearson"
                        corOptions = list(maxPOutliers = 0.1), # list(use = 'p')
                        verbose = 5) 

# Text way to assign (modified from SVA network paper)
if(!is.na(sft$powerEstimate)){
  print(paste("Soft power should be", sft$powerEstimate))
  wgcna_power <- sft$powerEstimate
}else if(is.na(sft$powerEstimate)){
  print(paste("no power reached r-squared cut-off, assing power based on number of samples"))
  wgcna_power <- 16
  print(paste("Soft power should be", wgcna_power))
}

# Plot the results:
pdf("soft_thresholding_power.pdf", height = 5, width =9)
#sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()

save(sft, wgcna_power, file = "WGCNA_sft.RData")
# load("WGCNA_sft.RData")

# Construct network and detect modules 
# Use identified soft threshold (power) above for power below or if no obvious plateau use Horvath's pre-caluclated (topic 6 from link below)
# https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/faq.html
net = blockwiseModules(WGCNA_data,
                       power = wgcna_power, # Change based on results of wgcna_power, 24
                       networkType = "signed", # signed is recommended for methylation
                       TOMType = "signed", 
                       corType = "bicor", # "bicor" is more powerful than "pearson", but also needs to be selected in soft thresholding  (https://www.ncbi.nlm.nih.gov/pubmed/23217028)
                       maxPOutliers = 0.1, # Forces bicor to never regard more than the specified proportion of samples as outliers (https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html)
                       minModuleSize = min(20, ncol(WGCNA_data)/2), # Functions default is used here
                       numericLabels = TRUE,
                       pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "regionTOM", 
                       loadTOM = TRUE,
                       blocks = NULL,
                       maxBlockSize = length(WGCNA_data), # split calculations into blocks or don't using = length(WGCNA_data), limit is = sqrt(2^31)
                       nThreads = 60,
                       verbose = 3)

# open a graphics window
#sizeGrWindow(12, 9)
pdf("module_dendogram.pdf", height = 5, width = 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

dev.off()

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "WGCNA-networkConstruction-auto.RData")

## WGCNA 3: Relating modules to phenotypes and identifying important regions ----

# load("WGCNA_1.RData")
# load("WGCNA-networkConstruction-auto.RData")

# Quantify module-trait associations 
# Define numbers of genes and samples
nGenes = ncol(WGCNA_data)
nSamples = nrow(WGCNA_data)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(WGCNA_data, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

# Remove Grey
MEs <- MEs %>%
  dplyr::select(-MEgrey)

save(datTraits, file = "daTraits.RData") # for mir663 analysis
save(MEs, datTraits, file = "hippocampus_MEs.RData")

### Infant hippocampus ----------------------------------------------------

load("hippocampus_MEs.RData")

datTraits <- datTraits %>% 
  tibble::rownames_to_column(var = "InfantID") %>%
  dplyr::left_join(readxl::read_xlsx("hippocampusLipids.xlsx") %>%
                     dplyr::mutate(InfantID = as.character(InfantID)) %>%
                     dplyr::select(InfantID,
                                   "C16:0 Palmitic Acid/Hexadecanoic Acid",
                                   "C18:0 Stearic Acid/Octadecanoic Acid",
                                   "C18:2 n-6 Linoleic Acid",
                                   "C20:3 Homo-y-Linolenic Acid/8,11,14-Eicosatrienoic Acid" = "C20:3 Homo-gamma-linolenic acid/8,11,14-eicosatrienoic acid",
                                   "C20:4 n-6 Arachidonic Acid" = "C20:4 n-6 Arachidonic acid" ,
                                   "C22:6/C24:1 Docosahexaenoic Acid (DHA)" = "C22:6 (DHA)/C24:1 Docosahexaenoic Acid (DHA)"
                     )) %>%
  dplyr::left_join(readr::read_csv("Infant_hippocampus_conc_adj.csv") %>%
                     dplyr::left_join(readr::read_csv("Infant_brain_meta.csv")) %>%
                     dplyr::select(Exp_ID,	InfantID = Animal_ID,
                                   "b-Hydroxybutyric Acid" = `3_Hydroxybutyrate`,
                                   Asparagine, Citrate, Glutathione, Glycerol, Guanosine, UMP) %>%
                     dplyr::mutate(InfantID = as.character(InfantID)) %>%
                     dplyr::select(-Exp_ID)) %>% 
  tibble::column_to_rownames(var = "InfantID")

moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, 25)

#### Heatmap --------------------------------------------------------------

pdf("hippocampus_tidy_module_trait_correlations.pdf", height = 4, width = 15)

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits) %>%
                 stringr::str_wrap(25),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               main = paste("Infant Hippocampus Module-trait Relationships"),
               verticalSeparator.x = c(4,6, 12),
               horizontalSeparator.y = c(2,3))

dev.off()

### Blue Module Plots ----------------------------------------------------

data <- MEs %>%
  tibble::rownames_to_column("ID") %>%
  dplyr::left_join(datTraits %>%
                     tibble::rownames_to_column("ID")) %>%
  dplyr::mutate(Group = dplyr::case_when(Obese == 1 ~ "Obese",
                                         Control == 1 ~ "Control",
                                         Restriction == 1 ~ "Restriction",
                                         Pravastatin == 1 ~ "Pravastatin")) %>%
  dplyr::mutate(Group = factor(Group, levels = c("Obese", "Pravastatin", "Restriction", "Control"))) %>%
  dplyr::as_tibble()

#### Dot plot -------------------------------------------------------------
# https://github.com/cemordaunt/comethyl/blob/0b2b0a124d5e6c118aa25431e3d7bf58ccfcea9a/R/Explore_Module_Trait_Correlations.R#L109

# blueModule <- data %>% 
#   ggplot(aes(x = Group, y = MEblue, group = Group, color = Group)) +
#   geom_boxplot() +
#   geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = -.75, aes(fill = Group)) +
#   ggplot2::theme_classic() +
#   ggplot2::theme(text = element_text(size = 10),
#                  axis.title = element_text(size = 10),
#                  axis.title.x = element_blank(),
#                  legend.position = "none") +
#   ylab("Eigengene") +
#   scale_color_manual(values=c("firebrick3", "darkgoldenrod1", "forestgreen", "mediumblue")) + 
#   scale_fill_manual(values=c("firebrick3", "darkgoldenrod1", "forestgreen", "mediumblue"))

blueModule <- data %>% 
  ggplot2::ggplot(aes(x = Group,
                      y = MEblue,
                      fill = Group,
                      group = Group,)) +
  ggplot2::geom_bar(stat = "summary",
                    position = position_dodge(),
                    color = "Black",
                    size = 0.25) +
  ggplot2::theme_classic() +
  ggplot2::theme(text = element_text(size = 8),
                 axis.title = element_text(size = 8),
                 axis.title.x = element_blank(),
                 axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
                 legend.position = "none") +
  scale_fill_manual(values = c("firebrick3", "darkgoldenrod1", "forestgreen", "mediumblue")) +
  ylab("Eigengene")

#### Scatter plot ---------------------------------------------------------
# Modified from: https://github.com/cemordaunt/comethyl/blob/0b2b0a124d5e6c118aa25431e3d7bf58ccfcea9a/R/Explore_Module_Trait_Correlations.R#L151

smoothWGCNA <- function(data = data,
                        trait = trait){
  data %>% 
    ggplot(ggplot2::aes_string(x = rlang::as_name(trait), y = "MEblue")) + # color = "Group"
    geom_smooth(method = MASS::rlm,
                formula = y ~ x,
                color = "#56B1F7",
                fill = "#336A98",
                se = FALSE) +
    geom_point(size = 1.5) +
    ggplot2::theme_classic() +
    ggplot2::theme(text = element_text(size = 8),
                   axis.title = element_text(size = 8),
                   axis.text.x = element_text(angle = 45, hjust = 1),
                   legend.position = "none") +
    #scale_color_manual(values=c("Obese" = "firebrick3", "Pravastatin" = "darkgoldenrod1", "Restriction" = "forestgreen", "Control" ="mediumblue")) +
    ylab("Eigengene") %>% 
    return()
}

#### Combine --------------------------------------------------------------

cowplot::plot_grid(blueModule,
                   smoothWGCNA(data, "`Recognition Memory: Abstract Stimuli`"),
                   smoothWGCNA(data, "`Recognition Memory: Social Stimuli`"),
                   smoothWGCNA(data, "`C18:2 n-6 Linoleic Acid`"),
                   smoothWGCNA(data, "Asparagine"),
                   smoothWGCNA(data, "Citrate"),
                   nrow = 1,
                   align = "h",
                   labels = c("B)", "C)", "D)", "E)", "F)", "G)"),
                   label_size = 12,
                   scale = 0.9) %>%
  ggplot2::ggsave("blueModule.pdf",
                  plot = ., 
                  width = 12,
                  height = 2)

### Maternal Serum --------------------------------------------------------

getModule <- function(GD,MEs){
  
  print(glue::glue("Loading {GD}"))
  
  tibble::tibble(InfantID = rownames(MEs)) %>%
    dplyr::left_join(readRDS("chr9_all_DUX4_Block_cffDNA_individual_smoothed_methylation.rds") %>% # "chr9_blue_dmrs_whole_regioncffDNA_individual_smoothed_methylation.rds"
                       dplyr::mutate(Timepoint = dplyr::case_when(Timepoint == "GD45" ~ "GD040",
                                                                  Timepoint == "GD90" ~ "GD090",
                                                                  TRUE ~ as.character(Timepoint))) %>% 
                       dplyr::filter(Timepoint == !!GD) %>% 
                       dplyr::select(InfantID = Infant, Timepoint, "DUX4 Block", "DUX4 Region 1",
                                     "DUX4 Region 3", "DUX4 Region 9", "DUX4 Region 19", "DUX4 Region 20") %>%
                       dplyr::mutate(InfantID = as.character(InfantID)) %>%
                       tidyr::pivot_wider(names_from = Timepoint)) %>%
    dplyr::left_join(readxl::read_xlsx("IDs.xlsx") %>%
                       dplyr::mutate(InfantID = as.character(InfantID))) %>% 
    dplyr::left_join(readr::read_csv("Maternal_plasma_meta.csv") %>%
                       dplyr::left_join(readr::read_csv("Maternal_plasma_conc_adj.csv"), by = "Exp") %>%
                       dplyr::filter(Plasma_color == "normal") %>% 
                       dplyr::filter(GD_targeted == !!GD) %>% 
                       dplyr::select(Exp, MotherID = Animal_ID,
                                     "Alkaline Phosphatase" = Alk_Phos,
                                     Triglyceride,
                                     Hematocrit,
                                     "Lymphocytes" = Lymphocytes_ul,
                                     "White Blood Cells" = WBC,
                                     "MCP-1 (CCL2)" = MCP_1,
                                     "IFN-g" = IFN_g,
                                     "IL-1RA" = IL_ra,
                                     "IL-2" = IL_2,
                                     "IL-8" = IL_8,
                                     "IL-10" = IL_10,
                                     "IL-13" = IL_13,
                                     "sCD40L" = sCD40L,
                                     "TGF-a" = TGFa,
                                     "hs-CRP" = hsCRP,
                                     "a-Ketoglutaric acid" = `2_Oxoglutarate`,
                                     "b-Hydroxybutyric acid" = `3_Hydroxybutyrate`,
                                     Arginine,
                                     "Acetoacetate" = Acetoacetate,
                                     "Choline" = Choline,
                                     Creatine,
                                     Glycine,
                                     Glutamine,
                                     "Myo-inositol" = myo_Inositol,
                                     Pyruvate,
                                     Succinate,
                                     Uridine)) %>%
    dplyr::select(-Exp) %>%
    dplyr::select(-MotherID) %>%
    dplyr::mutate(InfantID = as.character(InfantID)) %>%
    tibble::column_to_rownames(var = "InfantID") %>% 
    as.matrix() %>%
    WGCNA::cor(MEs %>%
                 dplyr::select(!!GD := MEblue),
               .,
               use = "p")
}

load("hippocampus_MEs.RData")

moduleTraitCor <- purrr::map(c("GD040", "GD090", "GD120", "GD150"),
                             getModule,
                             MEs = MEs) %>%
  do.call(rbind, .) %>%
  magrittr::set_rownames(c("Trimester 1",
                           "Trimester 2",
                           "Early Trimester 3",
                           "Late Trimester 3"))

moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, 25) %>%
  t()

moduleTraitCor <- moduleTraitCor %>%
  t()

#### Heatmap --------------------------------------------------------------

pdf("hippocampus_tidy_module_trait_correlations_maternal_plasma_select_blue.pdf", height = 15.5, width = 6) #50

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
               main = paste("Hippocampus Blue Module-trait Relationships with Maternal Blood") %>%
                 stringr::str_wrap(40),
               horizontalSeparator.y = c(6,11,21))

dev.off()

### Supplementary ---------------------------------------------------------

#### Brain ----------------------------------------------------------------

load("hippocampus_MEs.RData")

datTraits <- datTraits %>% 
  tibble::rownames_to_column(var = "InfantID") %>%
  dplyr::left_join(readxl::read_xlsx("hippocampusLipids.xlsx") %>%
                     dplyr::mutate(InfantID = as.character(InfantID))) %>%
  dplyr::left_join(readr::read_csv("Infant_hippocampus_conc_adj.csv") %>%
                     dplyr::left_join(readr::read_csv("Infant_brain_meta.csv") %>%
                                        dplyr::select(Exp_ID,	InfantID = Animal_ID) %>%
                                        dplyr::mutate(InfantID = as.character(InfantID))) %>%
                     dplyr::select(-Exp_ID)) %>% 
  tibble::column_to_rownames(var = "InfantID")

moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, 25)

##### Heatmap -------------------------------------------------------------

pdf("hippocampus_tidy_module_trait_correlations_long.pdf", height = 4, width = 50)

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits) %>%
                 stringr::str_wrap(25),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               main = paste("Infant Hippocampus Module-trait Relationships"))

dev.off()

#### Maternal Serum -------------------------------------------------------

getModule <- function(GD,MEs){
  
  print(glue::glue("Loading {GD}"))
  
  tibble::tibble(InfantID = rownames(MEs)) %>%
    dplyr::left_join(readRDS("chr9_all_DUX4_Block_cffDNA_individual_smoothed_methylation.rds") %>%
                       dplyr::mutate(Timepoint = dplyr::case_when(Timepoint == "GD45" ~ "GD040",
                                                                  Timepoint == "GD90" ~ "GD090",
                                                                  TRUE ~ as.character(Timepoint))) %>%
                       dplyr::filter(Timepoint == !!GD) %>%
                       dplyr::select(InfantID = Infant, Timepoint, contains("DUX4")) %>%
                       dplyr::mutate(InfantID = as.character(InfantID)) %>%
                       tidyr::pivot_wider(names_from = Timepoint)) %>%
    dplyr::left_join(readxl::read_xlsx("IDs.xlsx") %>%
                       dplyr::mutate(InfantID = as.character(InfantID))) %>%
    dplyr::left_join(readr::read_csv("Maternal_plasma_meta.csv") %>%
                       dplyr::select(Exp, EPV_b,EPV,WBC,RBC,Hemoglobin,Hematocrit,MCV,MCH,MCHC,Platelets,
                                     Seg_Neutrophils_percent,Seg_Neutrophils_per.ul,Lymphocytes_percent,
                                     Lymphocytes_ul,Monocytes_percent,Monocytes_per.ul,Eosinophils_percent,
                                     Eosinophils_per.ul,Plasma_Protein,Fibrinogen,Sodium,Potassium,
                                     Chloride,TCO2,Anion_Gap,Phosphorous,Calcium,BUN,Total_Protein,
                                     Albumin,ALT,AST,CPK,Alk_Phos,GGT,LDH,Cholesterol,Triglyceride,Bili_Total,
                                     Direct,hsCRP,GM_CSF,IFN_g,IL_1b,IL_ra,IL_2,IL_6,IL_8,IL_10,"IL_12/23_p40",IL_13,IL_15,
                                     IL_17a,MCP_1,MIP_1b,MIP_1a,sCD40L,TGFa,TNFa,VEGF,C_Peptide,GIP,Insulin,
                                     Insulin_uU.mL,Leptin,PP_53,PYY_54,Sx, Plasma_color, GD_targeted, Animal_ID) %>%
                       dplyr::left_join(readr::read_csv("Maternal_plasma_conc_adj.csv"), by = "Exp") %>%
                       dplyr::filter(Plasma_color == "normal") %>%
                       dplyr::filter(GD_targeted == !!GD) %>%
                       dplyr::rename(MotherID = Animal_ID)) %>%
    dplyr::select(-Exp, -Plasma_color, -GD_targeted, -MotherID) %>%
    dplyr::mutate(InfantID = as.character(InfantID)) %>%
    tibble::column_to_rownames(var = "InfantID") %>%
    as.matrix() %>%
    WGCNA::cor(MEs %>%
                 dplyr::select(!!GD := MEblue),
               .,
               use = "p")
}

load("hippocampus_MEs.RData")

moduleTraitCor <- purrr::map(c("GD040", "GD090", "GD120", "GD150"),
                             getModule,
                             MEs = MEs) %>%
  do.call(rbind, .) %>%
  magrittr::set_rownames(c("Trimester 1",
                           "Trimester 2",
                           "Early Trimester 3",
                           "Late Trimester 3"))

moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, 25) %>%
  t()

moduleTraitCor <- moduleTraitCor %>%
  t()

##### Heatmap -------------------------------------------------------------

pdf("hippocampus_tidy_module_trait_correlations_maternal_plasma_long_blue.pdf", height = 55, width = 6) #50

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
               main = paste("Hippocampus Blue Module-trait Relationships with Maternal Blood") %>%
                 stringr::str_wrap(40))

dev.off()

### DMR significance and module membership ----------------------------------

# DMR significance and module membership
trait <- "Obese"

# Define variable diagnosis containing the diagnosis column of datTrait
trait <- datTraits %>%
  dplyr::select(!!trait)

# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(WGCNA_data, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(WGCNA_data, trait, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(trait), sep="");
names(GSPvalue) = paste("p.GS.", names(trait), sep="");

# Intramodular analysis
# Choose module from correlation heatmap
module = "blue"
column = match(module, modNames);
moduleGenes = moduleColors==module;

# sizeGrWindow(7, 7);
# par(mfrow = c(1,1));
pdf(paste0("hippocampus_blue_module_membership_", names(trait), ".pdf"), height = 7, width = 7)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", names(trait)),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

dev.off()

# Create the starting data frame for saving
geneInfo0 = data.frame(moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for diagnosis
modOrder = order(-abs(cor(MEs, trait, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Obese))
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "WGCNA_DMR_Info_Hippocampus_Obese.csv")

#' tidyModules
#' @description Extract and tidy regions of methylation from significant WGCNA modules
#' @param sigModules Character vector with names of signficant modules 
#' @param geneInfo geneInfo dataframe from WGCNA analysis (end of tutorial 3)
#' @return Genomic ranges object with coordinates for regions of methylation and module names as meta information
#' @export tidyModules
tidyModules <- function(sigModules = sigModules,
                        geneInfo = geneInfo){
  geneInfo %>%
    rownames_to_column() %>%
    dplyr::as_tibble() %>%
    dplyr::filter(moduleColor %in% sigModules) %>%
    dplyr::select(rowname, moduleColor) %>%
    tidyr::separate(rowname, into = c('seqnames', 'start', 'end')) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE, ignore.strand = TRUE) %>%
    return()
}

sigModuleRanges <- tidyModules(sigModules = "blue", geneInfo = geneInfo)

#sigModuleRanges <- split(sigModuleRanges, sigModuleRanges$moduleColor)

seqlevelsStyle(sigModuleRanges) <- "UCSC" 

backgroundRanges <- makeGRangesFromDataFrame(smoothed)
seqlevelsStyle(backgroundRanges) <- "UCSC" 

save(sigModuleRanges,backgroundRanges, file = "Hippocampus_blue_module_ranges.RData")

#### Annotate ---------------------------------------------------------------

annotateRegions <- function(regions = sigRegions,
                            TxDb = TxDb,
                            annoDb = annoDb){
  
  print(glue::glue("Annotating {tidyRegions} regions from {tidyGenome} with gene symbols",
                   tidyRegions = length(regions),
                   tidyGenome = TxDb %>%
                     GenomeInfoDb::genome() %>%
                     unique()))
  
  if(class(TxDb) == "EnsDb"){
    
    seqlevelsStyle(regions) <- "Ensembl" # Work around for organism not supported
    
    TxDb <- TxDb %>%
      ensembldb::filter(GeneBiotypeFilter("protein_coding")) 
    #ensembldb::filter(~ symbol != "NA" & gene_name != "NA" & entrez != "NA" & tx_biotype == "protein_coding") 
  }
  
  regions %>% 
    ChIPseeker::annotatePeak(TxDb = TxDb,
                             annoDb = annoDb,
                             overlap = "all",
                             verbose = FALSE
    ) %>%
    dplyr::as_tibble() %>%
    dplyr::select("seqnames",
                  "start",
                  "end",
                  "width",
                  "annotation",
                  "distanceToTSS",
                  "SYMBOL",
                  "GENENAME",
                  "geneId") %>%
    dplyr::rename(geneSymbol = SYMBOL,
                  gene = GENENAME
    ) %>%
    dplyr::mutate(annotation = gsub(" \\(.*","", annotation)) %>% 
    return()
}


DMReport <- function(sigRegions = sigRegions,
                     regions = regions,
                     bs.filtered = bs.filtered,
                     coverage = coverage,
                     name = "DMReport"){
  cat("\n","Preparing HTML report...")
  
  stopifnot(class(sigRegions) == c("tbl_df", "tbl", "data.frame"))
  
  sigRegions %>%
    dplyr::select(seqnames,
                  start,
                  end,
                  width,
                  annotation,
                  distanceToTSS,
                  geneSymbol,
                  gene,
                  geneId,
                  Name) %>%
    gt::gt() %>%
    gt::tab_header(
      title = name,
      subtitle = glue::glue("There are {tidySigRegions} regions \\
             from {tidyRegions} background regions \\
             assayed at {coverage}x coverage.
             On average, the DMRs are {avgLength} bp long.", 
                            tidySigRegions = nrow(sigRegions),
                            tidyRegions = length(regions),
                            avgLength = mean(sigRegions$width) %>% round()
      )
    ) %>% 
    gt::fmt_number(
      columns = gt::vars("width"),
      decimals = 0) %>% 
    gt::as_raw_html(inline_css = FALSE) %>%
    write(glue::glue("{name}.html"))
  cat("Done", "\n")
}

##### Run -----------------------------------------------------------------

setwd("/share/lasallelab/Ben/Obesity/meta/WGCNA/Hippocampus")
load("Hippocampus_blue_module_ranges.RData")

.libPaths("/share/lasallelab/programs/DMRichR/R_4.0") 
#BiocManager::install("jorainer/ChIPseeker")
packages <- c("AnnotationHub", "ensembldb", "DMRichR", "org.Mmu.eg.db", "plyranges")
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))
AnnotationHub::setAnnotationHubOption("CACHE", "/share/lasallelab/programs/DMRichR/R_4.0")

TxDb <- AnnotationHub::AnnotationHub()[["AH83244"]]
annoDb <- "org.Mmu.eg.db"
coverage <-  1
genome <- "rheMac10"
load("../../../brains/DMRs/all/Hippocampus_bsseq.RData")

print("Annotating module ranges")
sigModuleRanges %>%
  annotateRegions(TxDb = TxDb,
                  annoDb = annoDb) %>%
  dplyr::mutate(Name = dplyr::case_when(geneSymbol != "NA" ~ geneSymbol,
                                        is.na(geneSymbol) ~ geneId)) %>% 
  dplyr::mutate(Name = dplyr::case_when(Name == "ENSMMUG00000060367" ~ "DUX4",
                                        Name == "LOC715898" ~ "GOLGA6C",
                                        Name == "LOC699789" ~ "DPY19L2",
                                        Name == "LOC715936" ~ "ZNF721",
                                        Name == "C14H11orf74" ~ "IFTAP",
                                        Name == "ENSMMUG00000001136" ~ "KLK13",
                                        Name == "ENSMMUG00000012704" ~ "ISOC2",
                                        TRUE ~ as.character(Name))) %T>%
  DMReport(regions = backgroundRanges,
           bs.filtered = bs.filtered,
           coverage = coverage,
           name = "WGCNA_hippocampus_blue_module") %>% 
  openxlsx::write.xlsx("DMRs_annotated_Ensembl_WGCNA_hippocampus_blue_module.xlsx")

bs.filtered.bsseq %>% 
  DMRichR::getSmooth(sigModuleRanges) %>%
  write.csv("blue_module_individual_smoothed_methylation.csv")

### Hub genes -------------------------------------------------------------

#setwd("/share/lasallelab/Ben/Obesity/meta/WGCNA/Hippocampus")
#load("WGCNA_1.RData")
#load("Hippocampus_tidy_datTraits.RData")
#load("WGCNA-networkConstruction-auto.RData")

hubGenes <- chooseTopHubInEachModule(
  WGCNA_data, 
  moduleColors, 
  omitColors = "grey", 
  power = 4, # https://support.bioconductor.org/p/46342/
  type = "signed") %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>% 
  purrr::set_names("module", "coordinate") %>% 
  tidyr::separate(coordinate, into = c('seqnames', 'start', 'end')) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE, ignore.strand = TRUE) 

hubGenes %>%
  annotateRegions(TxDb = TxDb,
                  annoDb = annoDb) %>%
  tibble::add_column(hubGenes$module) %>% 
  openxlsx::write.xlsx("Hubgenes_WGCNA_hippocampus.xlsx")

### Export to cytoscape ---------------------------------------------------

setwd("/share/lasallelab/Ben/Obesity/meta")
setwd("WGCNA/Hippocampus/")
load("regionTOM-block.1.RData") # Takes a long time to load
load("WGCNA_1.RData")
#load("Hippocampus_tidy_datTraits.RData")
load("WGCNA-networkConstruction-auto.RData")

# Select modules
modules = "blue"

# Select module probes
probes = names(WGCNA_data)
inModule = is.finite(match(moduleColors, modules))
modProbes = probes[inModule]

# Annotate module probes
TxDb <- AnnotationHub::AnnotationHub()[["AH83244"]]
annoDb <- "org.Mmu.eg.db"

# https://stackoverflow.com/questions/7659891/r-make-unique-starting-in-1
make.unique.2 = function(x, sep='.'){
  ave(x, x, FUN=function(a){if(length(a) > 1){paste(a, 1:length(a), sep=sep)} else {a}})
}

modGenes <- modProbes %>%
  dplyr::as_tibble() %>% 
  tidyr::separate(value, into = c('seqnames', 'start', 'end')) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE, ignore.strand = TRUE) %>%
  annotateRegions(TxDb = TxDb,
                  annoDb = annoDb) %>% 
  dplyr::mutate(nodeName = dplyr::case_when(geneSymbol != "NA" ~ geneSymbol,
                                            is.na(geneSymbol) ~ geneId)) %>% 
  dplyr::mutate(nodeName = dplyr::case_when(nodeName == "ENSMMUG00000060367" ~ "DUX4",
                                            nodeName == "LOC715898" ~ "GOLGA6C",
                                            nodeName == "LOC699789" ~ "DPY19L2",
                                            nodeName == "LOC715936" ~ "ZNF721",
                                            nodeName == "C14H11orf74" ~ "IFTAP",
                                            nodeName == "ENSMMUG00000001136" ~ "KLK13",
                                            nodeName == "ENSMMUG00000012704" ~ "ISOC2",
                                            nodeName == "ENSMMUG00000051744" ~ "Novel Gene 1",
                                            nodeName == "ENSMMUG00000052059" ~ "Novel Gene 2",
                                            nodeName == "ENSMMUG00000062616" ~ "Novel Gene 3",
                                            nodeName == "ENSMMUG00000031128" ~ "Novel Gene 4",
                                            nodeName == "ENSMMUG00000056847" ~ "Novel Gene 5",
                                            nodeName == "ENSMMUG00000060266" ~ "Novel Gene 6",
                                            nodeName == "ENSMMUG00000055454" ~ "Novel Gene 7",
                                            nodeName == "ENSMMUG00000055981" ~ "Novel Gene 8",
                                            nodeName == "ENSMMUG00000049528" ~ "Novel Gene 9",
                                            TRUE ~ as.character(nodeName))) %>% 
  dplyr::mutate(nodeName = make.unique.2(.$nodeName, sep = " ")) %>%
  purrr::pluck("nodeName")

# Select the corresponding Topological Overlap (This takes a long time)
# https://support.bioconductor.org/p/69715/
modTOM = as.matrix(TOM)[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
save(modTOM, file = "modTOM.RData")

# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = "Hippocampus_blue_CytoscapeInput-edges.txt",
                               nodeFile = "Hippocampus_blue_CytoscapeInput-nodes.txt",
                               weighted = TRUE,
                               threshold = 0.005, # 0.5 is default for function, 0.02 from tutorial
                               nodeNames = modGenes,
                               altNodeNames = modProbes,
                               nodeAttr = moduleColors[inModule])

# https://www.biostars.org/p/60896/
# On Cytoscape 3.8.0, you can import WGCNA network by File -> Import -> Network from File
# and selecting the module edge file. On the import dialogue box, you will typically 
# select the fromNode column as the Source Node and the toNode column as the Target Node. 
# The weight column should be left as an Edge Attribute. 
# The direction column should be changed to interaction type. 
# fromAltName is a Source Node Attribute while the toAltName is a Target Node Attribute.
