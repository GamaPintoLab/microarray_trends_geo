## Analysis of GEO database


BiocManager::install("rentrez")
library(rentrez)
library(tidyverse)
library(scales) # For formatting
library(ggbreak)

entrez_db_summary("gds")

#DbName: gds
#MenuName: GEO DataSets
#Description: GEO DataSets
#DbBuild: Build241112-1852.1
#Count: 7763677
#LastUpdate: 2024/11/12 19:45 

entrez_db_searchable("gds")

#Searchable fields for database 'gds'
#ALL 	 All terms from all searchable fields 
#UID 	 Unique number assigned to publication 
#FILT 	 Limits the records 
#ORGN 	 exploded organism names 
#ACCN 	 accession for GDS (DataSet), GPL (Platform), GSM (Sample), GSE (Series) 
#TITL 	 Words in title of record 
#DESC 	 Text from description, summary and other similar fields 
#SFIL 	 Supplementary Files 
#ETYP 	 Entry type (DataSet or Series)         <<<<<<<<<<<<<<
#STYP 	 Sample type 
#VTYP 	 type of values, e.g. log ratio, count 
#PTYP 	 Platform technology type               <<<<<<<<<<<<<<
#GTYP 	 type of dataset                        <<<<<<<<<<<<<<
#NSAM 	 Number of samples 
#SRC 	 sample source 
#AUTH 	 author of the GEO Sample, Platform or Series 
#INST 	 institute, or organization affiliatedd with contributers 
#NPRO 	 number of platform probes 
#SSTP 	 subset variable type 
#SSDE 	 subset description 
#GEID 	 name or identifier for the spot, e.g. GenBank, UniGene ID, Locus Link ID etc. 
#PDAT 	 publication date from the GEO related entities 
#UDAT 	 date                                   <<<<<<<<<<<<<<
#TAGL 	 Tag/Signature length for SAGE/MPSS 
#RGSE 	 Related Series 
#RGPL 	 Related Platform 
#MESH 	 Medical Subject Headings 
#PROJ 	 Project 
#ATNM 	 Attribute Name 
#ATTR 	 Attribute 
#PROP 	 Properties


### My option searh
##PTYP 	 Platform technology type 
# in situ oligonucleotide (1510611)
# mixed spotted oligonucleotide/cdna (352)
# oligonucleotide beads (638844)
# spotted dna/cdna (131981)
# spotted oligonucleotide (221705)
# high throughput sequencing (5072307)

##GTYP 	 type of dataset     
# expression profiling by array (72791)
# expression profiling by genome tiling array (763)
# expression profiling by high throughput sequencing (106298)
# genome binding/occupancy profiling by array (238)
# genome binding/occupancy profiling by genome tiling array (2391)
# genome binding/occupancy profiling by high throughput sequencing
# genome binding/occupancy profiling by snp array (23)
# genome variation profiling by array (884)
# genome variation profiling by genome tiling array (1650)
# genome variation profiling by high throughput sequencing (335)
# genome variation profiling by snp array (1561)
# methylation profiling by array (1697)
# methylation profiling by genome tiling array (2573)
# methylation profiling by high throughput sequencing (5866)
# methylation profiling by snp array (5708)
# non coding rna profiling by array (1697)
# non coding rna profiling by genome tiling array (112)
# non coding rna profiling by high throughput sequencing (8236)
# snp genotyping by snp array (900)


r_search <- entrez_search(db="gds", term="array [DTYP]")
r_search

#Entrez search result with 174932 hits (object contains 20 IDs and no web_history object)
#Search term (as translated):  array[All Fields]


search_year <- function(year, term){
  query <- paste(term, "AND (", year, "[PDAT])")
  entrez_search(db="gds", term=query, retmax=0)$count
}

year <- 2000:2024
datasets <- sapply(year, search_year, term="array", USE.NAMES=FALSE)

plot(year, datasets, type='b', main="The Rise and Fall of the Microarray")

class(datasets)
#[1] "integer"

datasets
#[1]     1    53   336  1872  5236  4396  5988  8908  8830  8892 10270  9652 12212 13139
#[15] 10214 14197 10111  9283  5667  8489  7162  6166  5696  4958  32

sum(datasets)
#[1] 174932



### Data per year per array type + NGS

# Define the search function
search_tech <- function(year, term) {
  # Format the query using the term for PTYP and the year for PDAT
  query <- paste(term, "[PTYP] AND (", year, "[PDAT])")
  entrez_search(db="gds", term=query, retmax=0)$count
}

# Define the year range
year <- 2000:2024

# Specify the platform technology types to search
technology_terms <- c("in situ oligonucleotide", 
                      "mixed spotted oligonucleotide/cdna", 
                      "oligonucleotide beads", 
                      "spotted dna/cdna", 
                      "spotted oligonucleotide",
                      "high throughput sequencing")

# Run the search for each year and technology term, storing results in a list
results <- lapply(technology_terms, function(term) {
  sapply(year, search_tech, term=term, USE.NAMES=FALSE)
})

# Convert the results list to a matrix for easier plotting
results_matrix <- do.call(cbind, results)
colnames(results_matrix) <- technology_terms

## editing for plotting

class(results_matrix)
class(year)

data <- as.data.frame(results_matrix)
data$year <- year

colnames(data) <- c("in situ oligo", "mixed spotted","oligo beads","spotted cDNA","spotted oligo", "HTS", "year")

write.csv(data, file="technology_type")

colSums(data)
#in situ oligo mixed spotted   oligo beads  spotted cDNA spotted oligo           HTS 
#   1510647           352        639186        131981        221705             5078363 

summary(data)

#in situ oligo    mixed spotted    oligo beads     spotted cDNA   spotted oligo  
#Min.   :     0   Min.   : 0.00   Min.   :    0   Min.   :    9   Min.   :    0  
#1st Qu.: 28024   1st Qu.: 0.00   1st Qu.:  759   1st Qu.:  585   1st Qu.: 4363  
#Median : 55311   Median : 0.00   Median :24150   Median : 2512   Median : 8312  
#Mean   : 60426   Mean   :14.08   Mean   :25567   Mean   : 5279   Mean   : 8868  
#3rd Qu.: 98408   3rd Qu.:10.00   3rd Qu.:44689   3rd Qu.: 9297   3rd Qu.:15189  
#Max.   :121065   Max.   :84.00   Max.   :59296   Max.   :18831   Max.   :23721  

#HTS               year     
#Min.   :      0   Min.   :2000  
#1st Qu.:      2   1st Qu.:2006  
#Median :  14750   Median :2012  
#Mean   : 203134   Mean   :2012  
#3rd Qu.: 362204   3rd Qu.:2018  
#Max.   :1418371   Max.   :2024

# Combine all data into a long format, including NGS
long_data <- data %>%
  pivot_longer(
    cols = 1:6, # Include HTS in the reshaping
    names_to = "Technology",
    values_to = "Value"
  )

long_data_filtered <- long_data %>% filter(year <= 2023)

# Define desired order for the technologies 
desired_order <- c("in situ oligo", "spotted oligo", "spotted cDNA", "mixed spotted", "oligo beads", "HTS")

long_data_filtered$Technology <- factor(long_data_filtered$Technology, levels = desired_order)

# Plot with scale break

ggplot(long_data_filtered, aes(x = year, y = Value, color = Technology, linetype = Technology)) +
  geom_line(linewidth = 0.75) +                               # Add lines
  geom_point(size = 2, shape = 18) +                        # Add open circles
  scale_y_break(
    c(125000, 207500),                                     # Break between 150,000 and 300,000
    scales = 1                                              # Expand higher range scale
  ) +
  scale_linetype_manual(
    values = c("HTS" = "dotted", "spotted cDNA" = "twodash", "spotted oligo" = "solid", "in situ oligo" = "solid", "mixed spotted" = "twodash", 
               "oligo beads" = "longdash")) +
  scale_color_manual(
    values = c("spotted cDNA" = "#0000CC", "spotted oligo" = "#1976D2", 
               "in situ oligo" = "#D32F2F", "mixed spotted" = "#90CAF9", 
               "oligo beads" = "#F1C40F", "HTS" = "#707B7C")) +
  labs(
    x = "Year",
    y = "Number of GEO Datasets",
    color = "Technology",
    linetype = "Technology",
    title = "GEO Datasets by Technology"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5))      




### Data per year per array assay type
##GTYP 	 type of dataset     

# expression profiling by array (72791)
# expression profiling by genome tiling array (763)
# genome binding/occupancy profiling by array (238)
# genome binding/occupancy profiling by genome tiling array (2391)
# genome binding/occupancy profiling by snp array (23)
# genome variation profiling by array (884)
# genome variation profiling by genome tiling array (1650)
# genome variation profiling by snp array (1561)
# methylation profiling by array (1697)
# methylation profiling by genome tiling array (2573)
# methylation profiling by snp array (5708)
# non coding rna profiling by array (1697)
# non coding rna profiling by genome tiling array (112)
# snp genotyping by snp array (900)

total <- 72791+763+238+2391+23+884+1650+1561+1697+2573+5708+112+900
total
#[1] 91291

dataset_classes <- c("expression profiling", "non-coding RNA profiling", "genome occupancy", "genome variation", "snp genotyping", "methylation profiling")

### Data per year per array type + NGS

# Define the search function
search_assay <- function(year, term) {
  # Format the query using the term for PTYP and the year for PDAT
  query <- paste(term, "[GTYP] AND (", year, "[PDAT])")
  entrez_search(db="gds", term=query, retmax=0)$count
}

# Define the year range
year <- 2000:2024

# Specify the assay types to search

assay_types <- c("expression profiling by array", 
"expression profiling by genome tiling array",
"genome binding/occupancy profiling by array", 
"genome binding/occupancy profiling by genome tiling array", 
"genome binding/occupancy profiling by snp array",
"genome variation profiling by array",
"genome variation profiling by genome tiling array",
"genome variation profiling by snp array",
"methylation profiling by array",
"methylation profiling by genome tiling array",
"methylation profiling by snp array",
"non coding rna profiling by array",
"non coding rna profiling by genome tiling array",
"snp genotyping by snp array")

short_assay_types <- c("expression_array", 
                 "expression_tiling",
                 "occupancy_array", 
                 "occupancy_tiling", 
                 "occupancy_snp_array",
                 "variation_array",
                 "variation_tiling",
                 "variation_snp_array",
                 "methylation_array",
                 "methylation_tiling",
                 "methylation_snp_array",
                 "ncrna_array",
                 "ncrna_tiling",
                 "snp_array")

assay_types <- unique(assay_types)

# Run the search for each year and technology term, storing results in a list
results_assay <- lapply(assay_types, function(term) {
  sapply(year, search_assay, term=term, USE.NAMES=FALSE)
})

# Convert the results list to a matrix for easier plotting
results_assay_matrix <- do.call(cbind, results_assay)
colnames(results_assay_matrix) <- short_assay_types
assay_data <- as.data.frame(results_assay_matrix)
assay_data$year <- year

d <- assay_data[,1:14]
sums <- colSums(d)
sum(d)
#[1] 91317

short_assay_data <- dplyr::mutate(assay_data, `expression profiling` = expression_array+ expression_tiling, 
                                  `ncRNA profiling` = ncrna_array+ ncrna_tiling,
                                  `genome occupancy` = occupancy_array + occupancy_tiling + occupancy_snp_array, 
                                  `genome variation`= variation_array + variation_tiling + variation_snp_array, 
                                  `snp genotyping` = snp_array, 
                                  `methylation profiling` = methylation_array + methylation_tiling + methylation_snp_array)

plot_data <- dplyr::select(short_assay_data, year, `expression profiling`, 
                           `ncRNA profiling`,
                           `genome occupancy`, 
                           `genome variation`, 
                           `snp genotyping`, 
                           `methylation profiling`)

write.csv(plot_data, file="assay_type")

# Combine all data into a long format
long_data2 <- plot_data %>%
  pivot_longer(
    cols = 2:7, 
    names_to = "Assay type",
    values_to = "Value"
  )

long_data2_filtered <- long_data2 %>% filter(year <= 2023)

# Define desired order for the assays 
desired_assay_order <- c("expression profiling", 
                         "ncRNA profiling",
                         "genome variation", 
                         "snp genotyping", 
                         "genome occupancy",
                         "methylation profiling")

long_data2_filtered$`Assay type` <- factor(long_data2_filtered$`Assay type`, levels = desired_assay_order)

# Plot with scale break

ggplot(long_data2_filtered, aes(x = year, y = Value, color = `Assay type`, linetype = `Assay type`)) +
  geom_line(linewidth = 0.75) +                               # Add lines
  geom_point(size = 2, shape = 18) +                        # Add open circles
  scale_y_break(
    c(750, 1250),                                     # Break between 150,000 and 300,000
    scales = 1                                              # Expand higher range scale
  ) +
  scale_linetype_manual(
    values = c("genome variation" = "solid", "snp genotyping" = "dotted", 
               "expression profiling" = "solid", "ncRNA profiling" = "dotted", 
               "genome occupancy" = "solid", "methylation profiling" ="dotted")) +
  scale_color_manual(
    values = c("genome variation" = "#1565C0", "snp genotyping" = "#64B5F6", 
               "expression profiling" = "#D32F2F", "ncRNA profiling" = "#E57373", 
               "genome occupancy" = "#F9A825", "methylation profiling" ="#FDD835" )) +
   labs(
    x = "Year",
    y = "Number of Datasets/Series",
    color = "Assay type",
    title = "GEO Datasets by Assay Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(vjust = +30),            # Push the title downward
    legend.position = "bottom",                                  # Place legend at the bottom
    legend.margin = margin(t = 15, unit = "pt"),        # Add space above the legend
    legend.title = element_text(size = 8, face = "bold"),                   # Bold the legend title
    legend.text = element_text(size = 8))


