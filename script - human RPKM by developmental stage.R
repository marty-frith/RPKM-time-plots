setwd("C:/R projects/Gene expression over development plots")

### load packages
library(tidyverse)
library(dplyr)
library(reshape2) # reshape2 is needed for the melt function
library(plotrix)
library(stringr)
library(biomaRt)

select <- dplyr::select # dplyr and biomaRt both have a select function, so defining the select function in this script as the dplyr version

### use biomaRt to add HGNC symbol column to data

Raw <- read_tsv("Data/brainseq_phase1_RPKM_all.txt") %>%   # read in gene expression data
  rename("ensembl_ID" = "...1") # rename column 1 to ensembl_ID

ensembl <- useEnsembl(biomart =  "genes", dataset = "hsapiens_gene_ensembl") #define data to be used 

genemap <-  getBM(attributes = c("ensembl_gene_id_version", "hgnc_symbol"), # biomaRt query, switch to gene id no version
                  filters = "ensembl_gene_id_version",
                  values = Raw$ensembl_ID,
                  mart = ensembl)

idx <- match(Raw$ensembl_ID, genemap$ensembl_gene_id_version) # assign matches between our ensembl IDs and hsapiens ensembl IDs to idx 
Raw$hgnc_symbol <- genemap$hgnc_symbol[ idx ]  #match ensembl IDs found in idx to their equivalent hgnc symbols in genemap and assign to column in Raw 




### select genes 
desiredGenes <- c("GRIN2A", "TRIO", "CACNA1C", "NADK") # user defines gene names here desiredStages 
desiredGenes <- sort(desiredGenes)
desiredGenesLength <- c(length(desiredGenes))

### data manipulation  
Meta <- read_tsv("Data/brainseq_phase1_sample-meta.txt")   #read in meta data

RawFiltered <- Raw %>%
  group_by(hgnc_symbol) %>%
  filter(hgnc_symbol %in% desiredGenes, .preserve = TRUE) %>%
  arrange(hgnc_symbol)# filter ensembl_ID column based on user input of genes desired 

RawFiltered <- data.frame(t(RawFiltered[-1]))  #transform dataframe so that ensembl_ID rows become columns 
colnames(RawFiltered) <- desiredGenes   #give new ensembl_ID columns the names of the user defined genes 
RawFiltered <-tibble::rownames_to_column(RawFiltered, "Sample_ID")#set row names to be a column called sample_ID

RawClean <- RawFiltered %>%
  select(-Sample_ID) %>%
  mutate_if(is.character, as.numeric, na.rm = TRUE) %>%
  mutate(Sample_ID = RawFiltered$Sample_ID)

Joined <- Meta %>%
  left_join(RawClean, by = c("Sample_ID")) %>%  #join raw data table and meta data table together
  mutate(across(Age, round, 3)) %>% #round numbers in age column to 3dp to make compatible with age bins
  mutate(age_bins = cut(Age, breaks = c(-0.864, -0.810, -0.648, -0.621, -0.432, -0.081, -0.001, 0, 0.499, 0.999, 5, 12, 19, 29, 59, 100))) %>%
  filter(age_bins != "(-0.081,-0.001]")


  
# age_bins --> Early Fetal = 0.864:-0.811, Early Midfetal = -0.81: -0.649, Midfetal = -0.648:-0.622, 
#            Late Midfetal = -0.621:-0.433, Late Fetal = -0.432:-0.081, Early Infancy = 0:0.499, Late Infancy = 0.5:0.999
#            Early Childhood = 1:5, Late Childhood = 6:12, Adolescence = 13:19, 
#            Early Adulthood = 20:29, Mid Adulthood = 30:59, Late Adulthood = 60:100)

#fetal stages calculated as 1 week back from birth = -0.027, calculate bottom of each bin and bring top of bin up to meet next bin
#did not include 38 weeks post conception to birth as was not included in example paper, also only one sample in this data set 
#e.g late midfetal = 17-23 weeks post conception, gestation = 40 weeks --> bottom = ((40-17) * 0.027) = -0.0621
#                                                                          top = bottom of late fetal bin - 0.001 (= 23 weeks, 6 days, 23 hours, etc.)


  
Means <- Joined %>%
  # filter(genotype == "WT") %>%   can choose to filter based on a condition here  
  select(age_bins, (desiredGenes)) %>%   #select age column and any column containing an ensemble ID
  group_by(age_bins) %>%
  summarise_each(funs(mean)) %>% 
  melt(id.vars = "age_bins", variable.name = "hgnc_symbol", value.name = "RPKM") 


StdError <- Joined %>%
  # filter(genotype == "WT") %>%   can choose to filter based on a condition here  
  select(age_bins, (desiredGenes)) %>%   #select age column and any column containing an ensemble ID
  group_by(age_bins) %>%
  summarise_each(funs(std.error)) %>%  
  melt(id.vars = "age_bins", variable.name = "hgnc_symbol", value.name = "se") 


PlotData <- left_join(Means, StdError) %>%
  filter(age_bins != "(0.499,0.999]") %>%
  mutate(se_min = RPKM - (se/2),
         se_max = RPKM + (se/2))
  
  

  
  ###plot data
  #create vectors for plot features
breakLabs <- c(0, 0.00098, 0.00195, 0.00390, 0.00781, 0.01563, 0.03125, 0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096) #breaks in log2 for yscale

##values plotted on log2 scale

geneExprPlot <- ggplot(PlotData, aes(x = age_bins, y=RPKM, color = hgnc_symbol, group = hgnc_symbol)) + #group and color by gene selected 
  geom_line(linewidth = 1.1) + #set line geometry
  geom_ribbon(aes(x=age_bins, ymin= se_min, ymax = se_max, color = hgnc_symbol, fill = hgnc_symbol), alpha = 0.25) +
  scale_y_continuous(trans = "log2") + #log transform y scale to improve readability
  coord_cartesian(ylim = c(min(PlotData$se_min), max(PlotData$se_max))) + #set limits for y scale, can change params to 'zoom' in wihtout data loss
  scale_x_discrete(labels = c("Early Midfetal", "Midfetal",  "Late Midfetal", "Late Fetal", "Early Infancy", 
                              "Early Childhood", "Late Childhood", "Adolescence", "Early Adulthood", "Late Adulthood")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), #adjusting text on x axis 
        panel.background = element_rect(fill = "white"), #background parameters
        panel.grid.major = element_line(color = "grey 95"), #major gridline parameters
        panel.grid.minor.y = element_line(color = "grey 96"),
        axis.line = element_line((color = "light grey")))

print(geneExprPlot)











# ## Values log2 transformed before plotting, <1 log2 transformation makes negativevalue, unsure if this is correct / how to plot
# 
#log2 of value +1 y: log2(rpkm+1)

# PlotDataLog2 <- bind_cols(SE, Means) %>%
#   select(!age_bins...4) %>%
#   select(!Standard_error) %>%   #remove columns not needed from the table
#   rename("age_bins" = "age_bins...1", "se" = "value...3", "RPKM" = "value...6")%>%   #rename columns for clarity
#   filter(age_bins != "(0.499,0.999]") %>%
#   arrange(ensembl_ID)%>%
#   mutate(RPKM_log2 = log2(RPKM),
#          se_log2 = log2(se),
#          se_min_log2 = RPKM_log2 - (se_log2/2),   #create min and max values for standard error ribbon
#          se_max_log2 = RPKM_log2 + (se_log2/2))
# 
# 
# 
# geneExprPlotLog2 <- ggplot(PlotDataLog2, aes(x = age_bins, y=RPKM_log2, color = ensembl_ID, group = ensembl_ID)) + #group and color by gene selected
#   geom_line(linewidth = 1.1) + #set line geometry
#   geom_ribbon(aes(x=age_bins, ymin = se_min_log2, ymax = se_max_log2, color = ensembl_ID, fill = ensembl_ID), alpha = 0.25) +
#   scale_x_discrete(labels = c("Early Midfetal", "Midfetal",  "Late Midfetal", "Late Fetal", "Early Infancy",
#                               "Early Childhood", "Late Childhood", "Adolescence", "Early Adulthood", "Late Adulthood")) +
#   coord_cartesian(ylim = c(max(PlotDataLog2$se_min_log2), min(PlotDataLog2$se_max_log2))) + #set limits for y scale, can change params to 'zoom' in wihtout data loss
#   #annotation_logticks(base = 2, sides = "l", scaled = FALSE) + #add log scale ticks
#   labs(x="Developmental Stage", y="log2(RPKM)", title = "log2(RPKM) by devlopmental stage" ) +
#   theme(axis.text.x = element_text(vjust = 1, hjust=1, angle = 45), #adjusting text on x axis
#         panel.background = element_rect(fill = "white"), #background parameters
#         panel.grid.major = element_line(color = "grey 95"), #major gridline parameters
#         panel.grid.minor.y = element_line(color = "grey 96"),
#         axis.line = element_line((color = "light grey"))) +
#   facet_wrap(vars(ensembl_ID), scales="free")   #facet these plots be gene because very hard to seperate plots w/ se ribbons when log2
# 
# 
# print(geneExprPlotLog2)