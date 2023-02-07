# user options 

WorkingDirectory <- ("C:/R projects/Gene expression over development plots")
rawData <- ("Data/Setd1a-Het-WT_E14-P70_AllGenes_rpkm.txt") # user gives path to data file
metaData <- ("Data/sample_meta.csv")
desiredGenes <- c("Trio", "Grin2a") # user defines gene names here desired
Log2TransRPKM <- TRUE


### load packages
setwd(WorkingDirectory)
library(tidyverse)
library(dplyr)
library(reshape2) # reshape2 is needed for the melt function
library(plotrix)
library(biomaRt)

select <- dplyr::select # dplyr and biomaRt both have a select function, so defining the select function in this script as the dplyr version

### use biomaRt to add gene symbol column to data
Meta <- read_csv(metaData)   #read in meta data

Raw <- read_tsv(rawData) %>%   # read in gene expression data
  rename("ensembl_ID" = "...1")# rename column 1 to ensembl_ID
Raw$ensembl_ID <- gsub("\\..*", "", Raw$ensembl_ID)


ensembl = useEnsembl(biomart =  "ensembl", dataset = "mmusculus_gene_ensembl") #define data to be used 

genemap <-  getBM(attributes = c("ensembl_gene_id", "external_gene_name"), # biomaRt query 
                  filters = "ensembl_gene_id",
                  values = Raw$ensembl_ID,
                  mart = ensembl)

idx <- match(Raw$ensembl_ID, genemap$ensembl_gene_id) # assign matches between our ensembl IDs and hsapiens ensembl IDs to idx 
Raw$gene_symbol <- genemap$external_gene_name[ idx ]  #match ensembl IDs found in idx to their equivalent hgnc symbols in genemap and assign to column in Raw 



### data manipulation  

RawFiltered <- Raw %>%
  group_by(gene_symbol) %>%
  filter(gene_symbol %in% desiredGenes, .preserve = TRUE) %>%
  arrange(gene_symbol)# filter ensembl_ID column based on user input of genes desired 

RawFiltered <- data.frame(t(RawFiltered[-1]))  #transform dataframe so that ensembl_ID rows become columns 
colnames(RawFiltered) <- desiredGenes   #give new ensembl_ID columns the names of the user defined genes 
RawFiltered <- rownames_to_column(RawFiltered, "Sample_ID")#set row names to be a column called sample_ID


Joined <- Meta %>%
  left_join(RawFiltered, by = c("sampleID" = "Sample_ID")) %>% #join raw data table and meta data table together
  rename("Sample_ID" = "sampleID") %>%
  mutate(across(desiredGenes, as.numeric))

Joined$age <- ordered(Joined$age, levels = c("E14", "E18", "P7", "P35", "P70"))

Means <- Joined %>%
  # filter(genotype == "WT") %>%   can choose to filter based on a condition here  
  select((desiredGenes), age) %>%   #select age column and any column containing an ensemble ID
  group_by(age) %>%
  summarise_each(funs(mean)) %>% 
  melt(id.vars = "age", variable.name = "gene_symbol", value.name = "RPKM") 


StdError <- Joined %>%
  # filter(genotype == "WT") %>%   can choose to filter based on a condition here  
  select((desiredGenes), age) %>%   #select age column and any column containing an ensemble ID
  group_by(age) %>%
  summarise_each(funs(std.error)) %>%  
  melt(id.vars = "age", variable.name = "gene_symbol", value.name = "se") 


PlotData <- left_join(Means, StdError) %>%
  mutate(se_min = RPKM - (se/2),
         se_max = RPKM + (se/2))

PlotDataLog2 <- left_join(StdError, Means) %>%
  mutate(RPKM_log2 = log2(RPKM +1),
         se_log2 = log2(se +1),
         se_min_log2 = RPKM_log2 - (se_log2/2),   #create min and max values for standard error ribbon
         se_max_log2 = RPKM_log2 + (se_log2/2))


###plot data
#create vectors for plot features
breakLabs <- c(0, 0.00098, 0.00195, 0.00390, 0.00781, 0.01563, 0.03125, 0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096) #breaks in log2 for yscale

if (Log2TransRPKM == TRUE) {
  
  geneExprPlotLog2 <- ggplot(PlotDataLog2, aes(x = age, y=RPKM_log2, color = gene_symbol, group = gene_symbol)) + #group and color by gene selected
    geom_line(linewidth = 1.1) + #set line geometry
    geom_ribbon(aes(x=age, ymax = se_min_log2, ymin = se_max_log2, color = gene_symbol, fill = gene_symbol), alpha = 0.25) +
    scale_x_discrete(labels = c("E14", "E18", "P7", "P35", "P70")) +
    labs(x="Developmental Stage", y="log2(RPKM + 1)", title = "log2(RPKM + 1) by Developmental Stage" ) +
    theme(axis.text.x = element_text(vjust = 1, hjust=1, angle = 45), #adjusting text on x axis
          panel.background = element_rect(fill = "white"), #background parameters
          panel.grid.major = element_line(color = "grey 95"), #major gridline parameters
          panel.grid.minor.y = element_line(color = "grey 96"),
          axis.line = element_line((color = "light grey")))
  print(geneExprPlotLog2)
  
}else{
  
  geneExprPlot <- ggplot(PlotData, aes(x = age, y=RPKM, color = gene_symbol, group = gene_symbol)) + #group and color by gene selected
    geom_line(linewidth = 1.1) + #set line geometry
    geom_ribbon(aes(x=age, ymin= se_min, ymax = se_max, color = gene_symbol, fill = gene_symbol), alpha = 0.25) +
    scale_y_continuous(trans = "log2", breaks = breakLabs) + #log transform y scale to improve readability
    coord_cartesian(ylim = c(min(PlotData$se_min), max(PlotData$se_max))) + #set limits for y scale, can change params to 'zoom' in wihtout data loss
    scale_x_discrete(labels = c("E14", "E18", "P7", "P35", "P70")) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), #adjusting text on x axis
          panel.background = element_rect(fill = "white"), #background parameters
          panel.grid.major = element_line(color = "grey 95"), #major gridline parameters
          panel.grid.minor.y = element_line(color = "grey 96"),
          axis.line = element_line(color = "light grey"),
          axis.title.x = element_blank())
  print(geneExprPlot)
}

